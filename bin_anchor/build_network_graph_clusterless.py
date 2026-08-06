#!/usr/bin/env python3
"""
Build the feature-correspondence graph for one cycle of the CLUSTERLESS workflow.

WHERE THIS SITS
    pair_features_1_clusterless.py   -> paired_features_edges.pkl   (edges)
    build_network_graph_clusterless.py  (THIS SCRIPT)
        -> accepted_components.json  (feature groups that survived this cycle)
        -> deleted_vertices.json     (vertices sent back for re-pairing next cycle)
    orchestrate_clusterless_cycles.py loops the two above, then
    clusters_json_to_consensusxml_clusterless.py emits the final consensusXML.

WHAT IT DOES
    1. Load the edges the pairing step emitted. One edge = "these two features, in
       two different files, are the same analyte".
    2. Build an igraph graph: one vertex per (filename, feature_idx).
    3. Take the CONNECTED COMPONENTS as the feature groups. That is what makes this
       "clusterless" - there is no Leiden/Louvain sub-clustering inside a component.
    4. Enforce the one-feature-per-file rule: a component may contain at most one
       vertex from any given input file. Where a file contributed several, all but
       one are deleted and handed to the next cycle.
    5. Export the surviving components and the deleted vertices.

DETERMINISM
    Output must be reproducible run to run. Two places make it so, and both matter:
      - export_accepted_components() sorts components by id (filtering results
        arrive in worker-completion order).
      - filter_one_component() sorts the kept vertex ids (set iteration order
        depends on how the set was built).
    Worker count, batching and load balancing do NOT affect the result.

ONE REPRESENTATION OF THE PARTITION
    graph.components() is the only place component membership is decided, and
    components[i] is already a plain list of vertex ids, so it is passed around
    directly rather than copied into another structure. Filtering returns
    {comp_idx: {'vertices': [...]}} covering every component - passed through
    untouched or with duplicates removed - which is exactly what the export
    consumes.
"""

import argparse
import pickle
import json
import math
from itertools import product
import numpy as np
import os
import sys
import gc
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
import time
from multiprocessing import Pool, Value, Lock

try:
    import igraph as ig
except ImportError:
    print("igraph not found. Install with: pip install python-igraph")
    exit(1)


# ============================================================================
# SECTION 1 - TUNING CONSTANTS
#
# None of these change the result. They trade memory against speed, or control
# how chatty the run is. The one knob that DOES change results is
# --mega_complexity_threshold, which is a CLI argument, not a constant.
# ============================================================================

# Worker processes used by the duplicate-vertex filtering stage when --num_workers
# is not given. 1 means "run in-process", i.e. the same path as
# --disable-multiprocessing: parallelism is opt-in, not the default.
DEFAULT_NUM_WORKERS = 1

# Each worker walks its own chunk of components and accumulates them until the
# summed load weight would exceed this limit, then builds the work items for that
# batch, processes them, and forces a gc.collect() before starting the next.
#
# This is purely a MEMORY-PACING knob: it does NOT affect IPC. The metadata has
# already been shipped to the worker by the time any of this runs - what crosses
# the process boundary is fixed earlier, where the metadata is split into one
# chunk per worker. Raising the limit gives fewer, larger batches (higher peak
# RAM, fewer collections); lowering it gives more, smaller batches.
#
# It also does NOT affect how components are distributed across workers - that is
# separate greedy load balancing over every component, independent of this value.
#
# Tune to the machine: larger with more RAM, smaller with more input files (more
# files -> bigger per-component vertex_attrs slices held for the whole batch).
#
# NOTE: one component's weight is bounded by mega_complexity_threshold, so once
# this limit exceeds that threshold every batch holds several components; and if
# a worker's entire chunk fits under the limit it forms a single batch.
COMBINATIONS_BATCH_LIMIT = 500000

# Per-component progress is logged only every Nth component (per worker process).
# Filtering routinely handles thousands of components per cycle, and logging each
# one buried the run in ~5 lines x N components of near-identical output.
# Anomalies (lost CPU affinity, missing graph) are ALWAYS logged, regardless.
LOG_EVERY_N_COMPONENTS = 250


# ============================================================================
# SECTION 2 - SHARED STATE FOR WORKER PROCESSES
#
# The filtering stage forks a Pool. Rather than pickling the graph and the
# vertex attributes into every worker via initargs (measured at ~2GB of copies),
# they are published here as module globals BEFORE the Pool is created, so
# workers inherit them through fork's copy-on-write. Set in
# filter_duplicate_file_vertices(); read in filter_one_component() and
# process_component_chunk().
#
# Consequence: this stage is fork-only. It will not work under spawn (macOS
# default, Windows). The pipeline runs in a Linux container, so that is fine.
# ============================================================================

graph = None
vertex_attrs = None

# Anchor-filter configuration, published the same way and for the same reason:
# every worker needs it and none of it is worth pickling per component. Set in
# filter_duplicate_file_vertices() from what main() assembled; None means the
# trimming methods are not in use. Keys:
#   method, radius_mode, anchor_ppm, anchor_rt_cutoff  - anchor_* methods
#   window_ppm, window_rt, window_max_candidates       - window method
#   bfs_max_depth, bfs_ppm_budget, bfs_rt_budget       - bfs method
#   range_x, range_y, mz_rt_weight_ratio               - the shared ranking metric
anchor_config = None

# Scaling constants as written by the pairing step, stashed by
# load_edges_from_paired_features(). A module global rather than a return value so
# the loader's signature stays as every existing caller expects it.
pairing_scaling_constants = None
pairing_cutoff_params = None

# Components handled so far BY THIS PROCESS. Workers are forked, so each gets its
# own copy of the counter and samples its own log lines independently.
_components_seen = 0


# ============================================================================
# SECTION 3 - INPUT: LOADING PAIRED-FEATURE EDGES
# ============================================================================

def load_edges_from_paired_features(paired_features_file: str) -> Tuple[List[Dict], Dict]:
    """
    Load edges and the file index map in a SINGLE file read operation.
    This avoids loading the pickle file twice and reduces peak memory usage.

    The pickle also carries a 'cutoff_params' block, but it is not read here: the
    pairing step already applied every cutoff before emitting an edge, so an edge
    existing at all means it passed. Re-filtering here could only ever be a no-op.

    Args:
        paired_features_file: Path to pickle or JSON file

    Returns:
        Tuple of (edges_list, file_index_map_dict)
        - edges_list: List of edge dictionaries (with index-based fields for memory efficiency)
        - file_index_map_dict: Dict mapping file indices to filenames (or empty if not found)
    """
    # Read the file exactly once; everything below works off `data`.
    if paired_features_file.endswith('.pkl'):
        with open(paired_features_file, 'rb') as f:
            data = pickle.load(f)
    elif paired_features_file.endswith('.json'):
        with open(paired_features_file, 'r') as f:
            data = json.load(f)
    else:
        raise ValueError("File must be .pkl or .json")

    # file_index_map: int file index -> filename. Edges reference files by index
    # to keep the pickle small; this is what turns them back into names. JSON
    # round-trips dict keys as strings, hence the conversion.
    file_index_map = {}
    if isinstance(data, dict):
        file_index_map = data.get('file_index_map', {})
        if file_index_map and isinstance(list(file_index_map.keys())[0], str):
            file_index_map = {int(k): v for k, v in file_index_map.items()}

    # The scaled space pairing matched in, for the anchor filter. Published as a
    # module global (SECTION 2) rather than returned, so this function's signature
    # is unchanged. Absent from pkls written before pairing recorded it, in which
    # case the anchor filter recomputes and warns.
    global pairing_scaling_constants, pairing_cutoff_params
    if isinstance(data, dict):
        pairing_scaling_constants = data.get('scaling_constants')
        # Not used to filter edges - pairing already applied these before emitting
        # any - but the anchor filter defaults its m/z tolerance to the pairing ppm.
        pairing_cutoff_params = data.get('cutoff_params')

    # Edges arrive in one of three historical shapes; flatten all of them to a
    # single flat list. The current pairing script emits the first one.
    edges = []
    if isinstance(data, dict):
        if 'network_edges' in data:
            # Wrapped: {'network_edges': {file_idx: [edges]} or [edges], ...}
            network_edges = data['network_edges']
            if isinstance(network_edges, dict):
                for file_idx, file_edges in network_edges.items():
                    if isinstance(file_edges, list):
                        edges.extend(file_edges)
            elif isinstance(network_edges, list):
                edges = network_edges
        elif all(isinstance(v, list) for v in data.values()):
            # Bare edges dictionary: {file_idx: [edges]} (int or str keys)
            for file_idx, file_edges in data.items():
                edges.extend(file_edges)
    elif isinstance(data, list):
        # Already a flat edge list
        edges = data

    return edges, file_index_map


# ============================================================================
# SECTION 4 - GRAPH CONSTRUCTION
# ============================================================================

def build_network_graph(edges: List[Dict], ident_lookup: Dict = None,
                        file_index_map: Dict = None) -> Tuple[ig.Graph, Dict, Dict]:
    """
    Build an igraph network graph from the paired-feature edges.

    Every edge is used. The pairing step already enforced the m/z and RT cutoffs
    before emitting an edge, so there is nothing left to filter here.

    Runs over the edge list TWICE, and has to: pass 1 discovers the vertex set so
    ids can be assigned in sorted (deterministic) order, pass 2 needs those ids to
    emit the edge list. See the note above pass 2 on the parsing asymmetry.

    Args:
        edges: List of edge dictionaries with index-based format (memory efficient):
               {'ref_file_idx': int, 'match_file_idx': int, 'current_feature_idx': int,
                'matched_feature_idx': int, 'distance': float, ...}
               Also supports legacy nested format:
               {'current_file': {filename, feature_idx}, 'matched_file': {...}, 'distance': float, ...}
        ident_lookup: Combined ident lookup for pep_ident/prot_ident metadata
        file_index_map: Maps file index -> filename

    Returns:
        Tuple of (igraph.Graph, vertex_id_map, vertex_attrs)
    """
    print(f"Loaded {len(edges)} total edges")

    if not edges:
        print("No edges to build a graph from!")
        return None, {}, {}

    # A vertex is one feature in one file, keyed by (filename, feature_idx).
    vertices = set()
    vertex_id_map = {}  # (filename, feature_idx) -> vertex_id
    vertex_attrs = {}   # vertex_id -> {'filename': str, 'feature_idx': int, 'label': str, ...}

    if file_index_map is None:
        file_index_map = {}

    # ---- Pass 1: collect the vertex set -----------------------------------
    for edge in edges:
        # Index-based format (memory efficient - what pairing emits today)
        if 'ref_file_idx' in edge and 'match_file_idx' in edge:
            ref_file_idx = edge['ref_file_idx']
            match_file_idx = edge['match_file_idx']
            current_feat_idx = edge.get('current_feature_idx', -1)
            matched_feat_idx = edge.get('matched_feature_idx', -1)

            if ref_file_idx in file_index_map and match_file_idx in file_index_map:
                file1_name = file_index_map[ref_file_idx]
                file2_name = file_index_map[match_file_idx]
            else:
                # Fallback: synthesise names so the graph still builds
                file1_name = f"file_{ref_file_idx}"
                file2_name = f"file_{match_file_idx}"

            vertex1_key = (file1_name, current_feat_idx)
            vertex2_key = (file2_name, matched_feat_idx)

        # Legacy nested format
        elif 'current_file' in edge and 'matched_file' in edge:
            file1_info = edge['current_file']
            file2_info = edge['matched_file']
            vertex1_key = (file1_info['filename'], file1_info['feature_idx'])
            vertex2_key = (file2_info['filename'], file2_info['feature_idx'])

        # Older 'file1'/'file2' format
        elif 'file1' in edge and 'file2' in edge:
            file1_info = edge['file1']
            file2_info = edge['file2']
            vertex1_key = (file1_info['filename'], file1_info['feature_idx'])
            vertex2_key = (file2_info['filename'], file2_info['feature_idx'])

        else:
            print(f"Warning: Edge with unknown format (missing file identifiers): {list(edge.keys())}")
            continue

        vertices.add(vertex1_key)
        vertices.add(vertex2_key)

    # ---- Assign vertex ids and attach per-vertex metadata ------------------
    # sorted() is what makes vertex ids deterministic: a set's iteration order is
    # not, and every downstream id (component membership, deleted lists) derives
    # from this numbering.
    for idx, vertex_key in enumerate(sorted(vertices)):
        vertex_id_map[vertex_key] = idx
        filename, feature_idx = vertex_key
        # Defaults first, then overwrite from the lookup. Keeping every vertex on
        # the same key set means the exported JSON has uniform records even for
        # features that were never identified.
        vertex_attrs[idx] = {
            'filename': filename,
            'feature_idx': feature_idx,
            'label': f"{Path(filename).stem}_{feature_idx}",
            'pep_ident': [],
            'prot_ident': [],
            'x_center': 0.0,   # m/z
            'y_center': 0.0,   # RT
            'intensity': 0.0,
            'openms_fid': '',
            'ms2_scans': '',
            'charge': 0
        }

        if ident_lookup:
            lookup_key = f"{filename}_{feature_idx}"
            if lookup_key in ident_lookup:
                ident_data = ident_lookup[lookup_key]
                vertex_attrs[idx]['pep_ident'] = ident_data.get('pep_ident', [])
                vertex_attrs[idx]['prot_ident'] = ident_data.get('prot_ident', [])
                vertex_attrs[idx]['x_center'] = ident_data.get('x_center', 0.0)
                vertex_attrs[idx]['y_center'] = ident_data.get('y_center', 0.0)
                vertex_attrs[idx]['intensity'] = ident_data.get('intensity', 0.0)
                vertex_attrs[idx]['openms_fid'] = ident_data.get('openms_fid', '')
                vertex_attrs[idx]['ms2_scans'] = ident_data.get('ms2_scans', '')
                vertex_attrs[idx]['charge'] = ident_data.get('charge', 0)

    print(f"Created {len(vertices)} unique vertices from edges")
    sys.stdout.flush()

    # How many vertices actually picked up a peptide identification, plus a couple
    # of samples - the usual first thing to check when identifications go missing.
    if ident_lookup:
        pep_ident_loaded_count = 0
        for idx, attrs in vertex_attrs.items():
            pep = attrs.get('pep_ident', [])
            if pep and len(pep) > 0:
                pep_ident_loaded_count += 1

        print(f"  ✓ Loaded pep_ident for: {pep_ident_loaded_count}/{len(vertex_attrs)} vertices")
        sys.stdout.flush()

        shown = 0
        for idx, attrs in vertex_attrs.items():
            pep = attrs.get('pep_ident', [])
            if pep and len(pep) > 0:
                print(f"    Vertex {idx}: {attrs['label']} -> pep_ident={pep}")
                shown += 1
                if shown >= 3:
                    break
        sys.stdout.flush()

    # ---- Materialise the graph --------------------------------------------
    graph = ig.Graph(len(vertices))

    for vertex_id, attrs in vertex_attrs.items():
        for attr_name, attr_value in attrs.items():
            graph.vs[vertex_id][attr_name] = attr_value

    # ---- Pass 2: translate edges into (vertex_id, vertex_id) pairs ---------
    # ASYMMETRY WITH PASS 1 (pre-existing, deliberately left alone): this pass is
    # laxer - it takes filenames via file_index_map.get() per index and treats
    # "no ref_file_idx and no file1" as the legacy nested format rather than
    # skipping. On well-formed input (index format + complete map) the two passes
    # agree exactly; they can only diverge on malformed edges, where pass 1 skips
    # and this one raises KeyError.
    edge_list = []
    edge_weights = []

    for edge in edges:
        if 'ref_file_idx' in edge and 'match_file_idx' in edge:
            ref_file_idx = edge['ref_file_idx']
            match_file_idx = edge['match_file_idx']
            file1_name = file_index_map.get(ref_file_idx, f'file_{ref_file_idx}')
            file2_name = file_index_map.get(match_file_idx, f'file_{match_file_idx}')
            vertex1_key = (file1_name, edge['current_feature_idx'])
            vertex2_key = (file2_name, edge['matched_feature_idx'])
        elif 'file1' in edge:
            file1_info = edge['file1']
            file2_info = edge['file2']
            vertex1_key = (file1_info['filename'], file1_info['feature_idx'])
            vertex2_key = (file2_info['filename'], file2_info['feature_idx'])
        else:
            file1_info = edge['current_file']
            file2_info = edge['matched_file']
            vertex1_key = (file1_info['filename'], file1_info['feature_idx'])
            vertex2_key = (file2_info['filename'], file2_info['feature_idx'])

        edge_list.append((vertex_id_map[vertex1_key], vertex_id_map[vertex2_key]))
        edge_weights.append(edge['distance'])

    # Bulk-add: one add_edges() call is far cheaper than N individual ones.
    graph.add_edges(edge_list)

    for edge_id, weight in enumerate(edge_weights):
        graph.es[edge_id]['distance'] = weight

    print(f"Created graph with {graph.vcount()} vertices and {graph.ecount()} edges")
    print("  Edge direction: current_file (source) → matched_file (target)")

    return graph, vertex_id_map, vertex_attrs


# ============================================================================
# SECTION 5 - DIAGNOSTIC ANALYSIS
#
# Everything in this section is DIAGNOSTIC OUTPUT ONLY. It is computed, printed
# and optionally written to JSON, but no downstream step of the clusterless
# workflow reads any of it - it exists for manual inspection of a run.
# --skip-analysis turns the whole section off.
# ============================================================================

def analyze_graph(graph: ig.Graph, vertex_attrs: Dict) -> Dict:
    """
    Summary statistics for the whole graph: degree distribution, components,
    clustering coefficient, diameter.

    DIAGNOSTIC OUTPUT ONLY - see the SECTION 5 header.

    COST NOTE: computes graph.diameter() on the largest connected component, which
    is the expensive part on big cycles. This is the main reason --skip-analysis
    exists.
    """

    if graph is None or graph.ecount() == 0:
        return {}

    # Global degree statistics
    degrees = graph.degree()
    analysis = {
        'num_vertices': graph.vcount(),
        'num_edges': graph.ecount(),
        'density': graph.density(),
        'num_components': len(graph.components()),
        'avg_degree': float(np.mean(degrees)),
        'std_degree': float(np.std(degrees)),
        'max_degree': int(max(degrees)),
        'min_degree': int(min(degrees)),
    }

    # Component size distribution
    components = graph.components()
    analysis['components'] = {
        'num_components': len(components),
        'component_sizes': [len(c) for c in components],
        'largest_component_size': max(len(c) for c in components) if components else 0
    }

    # Same degree statistics, broken down per input file - a file whose features
    # connect far less than the others usually means an alignment problem.
    file_degrees = {}
    for vertex_id in range(graph.vcount()):
        filename = vertex_attrs[vertex_id]['filename']
        degree = graph.degree(vertex_id)
        if filename not in file_degrees:
            file_degrees[filename] = []
        file_degrees[filename].append(degree)

    analysis['degree_by_file'] = {
        filename: {
            'num_vertices': len(degrees),
            'avg_degree': float(np.mean(degrees)),
            'std_degree': float(np.std(degrees)),
            'max_degree': int(max(degrees)) if degrees else 0,
            'min_degree': int(min(degrees)) if degrees else 0
        }
        for filename, degrees in file_degrees.items()
    }

    # Both of these can fail or be undefined on degenerate graphs; a missing
    # diagnostic must never take the run down, hence the bare excepts.
    try:
        cc = graph.transitivity_undirected()
        analysis['clustering_coefficient'] = cc
    except:
        analysis['clustering_coefficient'] = None

    try:
        if graph.is_connected():
            analysis['diameter'] = graph.diameter()
        else:
            # Disconnected: report the diameter of the largest component instead
            largest_comp = graph.induced_subgraph(components[0])
            analysis['diameter'] = largest_comp.diameter() if len(components[0]) > 1 else 0
    except:
        analysis['diameter'] = None

    return analysis


def analyze_graph_composition(graph: ig.Graph, vertex_attrs: Dict) -> Dict:
    """
    Break down what the components are MADE OF: how many features, from how many
    distinct files, and with what multiplicity.

    DIAGNOSTIC OUTPUT ONLY - see the SECTION 5 header.

    This is the view that tells you whether pairing behaved. A healthy run is
    dominated by components of size N drawn from N distinct files (key "N_N");
    a component keyed "6_4" holds 6 features from only 4 files, so two files
    contributed duplicates that the filtering stage will have to resolve.

    Returns:
    {
        "component_composition": {
            "3_3": 5,  # 5 components with 3 features from 3 files (aggregated)
            ...
        },
        "composition_summary": {... overall stats ...},
        "file_match_uniqueness": {
            "2": {  # keyed by component_size (2 features per component)
                "all_unique": 4,  # components where each file has exactly 1 feature
                "with_duplicates": 1,  # components where any file has 2+ features
                "total": 5,
                "file_composition": {
                    "file1.tsv, file2.tsv": 4,  # which files co-occur -> count
                    "file1.tsv, file3.tsv": 1,
                },
                "multiplicity_distribution": {
                    "1x2": 4,  # 4 components with pattern: each file has 1 feature
                    "2x1": 1,  # 1 component: one file has 2
                },
                "file_appearance": {
                    "file1.tsv": 5,  # how many components of size 2 include this file
                    "file2.tsv": 4,
                    "file3.tsv": 1,
                }
            },
            "3": {...},
            ...
        }
    }
    """
    if graph is None or graph.vcount() == 0:
        return {
            'component_composition': {},
            'composition_summary': {},
            'file_match_uniqueness': {}
        }

    components = graph.components()
    composition_map = {}                # "{num_features}_{num_files}" -> count
    file_match_uniqueness_by_size = {}   # str(component_size) -> detailed breakdown
    features_per_component = []
    files_per_component = []
    components_by_size = {}              # str(num_features) -> count

    for component_vertices in components:
        # How many features did each file contribute to this component?
        file_counts = {}
        for vertex_id in component_vertices:
            filename = vertex_attrs[vertex_id]['filename']
            file_counts[filename] = file_counts.get(filename, 0) + 1

        num_features = len(component_vertices)
        num_files = len(file_counts)

        # Aggregate shape counter, e.g. "6_4" = 6 features drawn from 4 files
        comp_key = f"{num_features}_{num_files}"
        composition_map[comp_key] = composition_map.get(comp_key, 0) + 1

        features_per_component.append(num_features)
        files_per_component.append(num_files)

        size_key = str(num_features)
        components_by_size[size_key] = components_by_size.get(size_key, 0) + 1

        if size_key not in file_match_uniqueness_by_size:
            file_match_uniqueness_by_size[size_key] = {
                'all_unique': 0,
                'with_duplicates': 0,
                'total': 0,
                'file_composition': {},
                'multiplicity_distribution': {},
                'file_appearance': {}
            }

        file_match_uniqueness_by_size[size_key]['total'] += 1

        # "all unique" = already satisfies the one-feature-per-file rule, so the
        # filtering stage has nothing to do for this component.
        multiplicities = sorted(file_counts.values(), reverse=True)
        is_all_unique = all(count == 1 for count in multiplicities)

        if is_all_unique:
            file_match_uniqueness_by_size[size_key]['all_unique'] += 1
        else:
            file_match_uniqueness_by_size[size_key]['with_duplicates'] += 1

        # Which set of files this component spans
        files_in_group = sorted(set(file_counts.keys()))
        file_comp_str = ', '.join(files_in_group)
        if file_comp_str not in file_match_uniqueness_by_size[size_key]['file_composition']:
            file_match_uniqueness_by_size[size_key]['file_composition'][file_comp_str] = 0
        file_match_uniqueness_by_size[size_key]['file_composition'][file_comp_str] += 1

        # Multiplicity shape, e.g. "2x1,1x2" = one file gave 2 features and two
        # files gave 1 each. Read as "<features per file>x<number of such files>".
        multiplicity_pattern = ','.join([f"{mult}x{multiplicities.count(mult)}"
                                        for mult in sorted(set(multiplicities), reverse=True)])
        if multiplicity_pattern not in file_match_uniqueness_by_size[size_key]['multiplicity_distribution']:
            file_match_uniqueness_by_size[size_key]['multiplicity_distribution'][multiplicity_pattern] = 0
        file_match_uniqueness_by_size[size_key]['multiplicity_distribution'][multiplicity_pattern] += 1

        # How often each file turns up in components of this size
        for fname in file_counts.keys():
            if fname not in file_match_uniqueness_by_size[size_key]['file_appearance']:
                file_match_uniqueness_by_size[size_key]['file_appearance'][fname] = 0
            file_match_uniqueness_by_size[size_key]['file_appearance'][fname] += 1

    summary = {
        'total_components': len(components),
        'total_features': sum(features_per_component),
        'avg_features_per_component': float(np.mean(features_per_component)) if features_per_component else 0,
        'avg_files_per_component': float(np.mean(files_per_component)) if files_per_component else 0,
        'max_features_in_component': max(features_per_component) if features_per_component else 0,
        'max_files_in_component': max(files_per_component) if files_per_component else 0,
        'min_features_in_component': min(features_per_component) if features_per_component else 0,
        'min_files_in_component': min(files_per_component) if files_per_component else 0,
        'components_by_size': {k: components_by_size[k] for k in sorted(components_by_size.keys(), key=int)},
    }

    return {
        'component_composition': composition_map,
        'composition_summary': summary,
        'file_match_uniqueness': file_match_uniqueness_by_size
    }


def print_analysis(analysis: Dict):
    """
    Pretty print analyze_graph() results to stdout.

    DIAGNOSTIC OUTPUT ONLY - see the SECTION 5 header.
    """
    if not analysis:
        print("No graph data to analyze")
        return

    print(f"\n{'='*70}")
    print("Network Graph Analysis")
    print(f"{'='*70}")

    print("\nBasic Statistics:")
    print(f"  Vertices: {analysis.get('num_vertices', 'N/A')}")
    print(f"  Edges: {analysis.get('num_edges', 'N/A')}")
    print(f"  Density: {analysis.get('density', 'N/A'):.4f}")
    print(f"  Average Degree: {analysis.get('avg_degree', 'N/A'):.2f} ± {analysis.get('std_degree', 'N/A'):.2f}")
    print(f"  Max Degree: {analysis.get('max_degree', 'N/A')}")
    print(f"  Min Degree: {analysis.get('min_degree', 'N/A')}")

    if 'components' in analysis:
        print("\nComponents:")
        print(f"  Number of Components: {analysis['components']['num_components']}")
        print(f"  Largest Component: {analysis['components']['largest_component_size']} vertices")
        if len(analysis['components']['component_sizes']) <= 10:
            print(f"  Component Sizes: {analysis['components']['component_sizes']}")

    if 'clustering_coefficient' in analysis:
        print(f"\nClustering Coefficient: {analysis['clustering_coefficient']}")

    if 'diameter' in analysis:
        print(f"Diameter: {analysis['diameter']}")

    if 'degree_by_file' in analysis:
        print("\nDegree Statistics by File:")
        for filename, stats in analysis['degree_by_file'].items():
            print(f"  {Path(filename).stem}:")
            print(f"    Vertices: {stats['num_vertices']}")
            print(f"    Avg Degree: {stats['avg_degree']:.2f}")
            print(f"    Max Degree: {stats['max_degree']}")
            print(f"    Min Degree: {stats['min_degree']}")


# ============================================================================
# SECTION 6 - DUPLICATE-VERTEX FILTERING
#
# THE RULE: a component may hold at most one feature per input file. Where a file
# contributed several, exactly one is kept and the rest are deleted (and re-paired
# in the next cycle).
#
# THE CHOICE: which one to keep. Two methods, picked per component by candidate
# count vs --mega_complexity_threshold:
#   METHOD 1  EXHAUSTIVE - score every combination, take the joint best
#             (in filter_one_component)
#   METHOD 2  HEURISTIC  - decide each duplicate set on its own, greedily
#             (in generate_filter_candidates)
# Each is documented in full at its implementation.
#
# THE MACHINERY, bottom-up in reading order:
#   6.1  resolve_worker_count / _init_worker    - how much parallelism, pinned where
#   6.2  count_filter_candidates / generate_filter_candidates  - what to choose between
#   6.3  filter_one_component                   - decide one component
#   6.4  process_component_chunk                - a worker's whole share, memory-paced
#   6.5  filter_duplicate_file_vertices         - orchestrator: scan, schedule, collect
#
# ============================================================================

# ---------------------------------------------------------------------------
# 6.1 Parallelism setup
# ---------------------------------------------------------------------------

def resolve_worker_count(requested: int = None) -> Tuple[int, List[int]]:
    """
    SINGLE SOURCE OF TRUTH for how much parallelism the filtering stage uses.

    One process per worker, each pinned to one CPU - so the worker count is the
    only number this stage needs, and every caller gets it from here.

    Worker count does NOT affect results. Output is invariant to how the work is
    distributed (components are sorted at the export boundary), so this is purely
    a performance knob and is safe to tune.

    A count of 1 runs the filtering in-process, taking the same path as
    --disable-multiprocessing. That is deliberate: a 1-worker Pool still pays the
    fork and IPC cost with nothing to overlap, and measured ~23% SLOWER than
    running in-process (3.24s vs 2.64s on cycle 3).

    Args:
        requested: Explicit worker count (--num_workers). None/<=0 -> DEFAULT_NUM_WORKERS
            (1, i.e. in-process).

    Returns:
        (num_workers, cpu_list) - cpu_list comes from sched_getaffinity so it
        respects cgroup/container limits, and is what workers are pinned to.
    """
    cpu_list = sorted(os.sched_getaffinity(0)) or [0]
    available = len(cpu_list)

    if requested is None or requested <= 0:
        num_workers, source = DEFAULT_NUM_WORKERS, f"default={DEFAULT_NUM_WORKERS}"
    else:
        num_workers, source = requested, "--num_workers"

    if num_workers > available:
        print(f"[CPU] Requested {num_workers} workers but only {available} CPU(s) "
              f"available to this process - clamping to {available}", flush=True)
        num_workers = available

    print(f"[CPU] Filtering with {num_workers} worker process(es) ({source}); "
          f"CPUs available: {cpu_list}", flush=True)
    return num_workers, cpu_list


def _init_worker(worker_index, cpu_list):
    """
    Pool initializer, run once inside each freshly forked worker.

    Its only job is CPU pinning: one worker per CPU, so the OS cannot migrate two
    workers onto the same core. Nothing else needs passing in - the large shared
    structures (SECTION 2) are already inherited through fork.

    Args:
        worker_index: Unique index for this worker (0-based), handed out by an
            atomic counter in the parent so two workers cannot claim the same CPU.
        cpu_list: CPUs this process may use, from resolve_worker_count(). Taken
            from sched_getaffinity, so pinning stays inside cgroup/container limits.
    """
    pid = os.getpid()
    cpu_to_pin = cpu_list[worker_index % len(cpu_list)]

    try:
        os.sched_setaffinity(0, {cpu_to_pin})
        print(f"[WORKER {pid}] pinned to CPU {cpu_to_pin}", file=sys.stderr, flush=True)
    except Exception as e:
        print(f"[WORKER {pid}] CPU pinning failed: {e}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# 6.2 Candidate generation - what a component gets to choose between
# ---------------------------------------------------------------------------

def component_edge_table(comp_vertices: List[int]) -> Tuple[List[int], Dict, Dict]:
    """
    The component's DISTINCT connections, for scoring a kept set.

    Returns (sorted_vertices, position_map, {(pos_a, pos_b): distance}).

    Parallel edges are collapsed. A mutual nearest-neighbour pair is emitted twice
    by pairing - once per direction, with identical distance - and
    induced_subgraph hands back both. Counting records instead of connections
    would weight a mutual match double, which measurably trades real connections
    away for mutuality, so the pair is counted once here.
    """
    cs = sorted(comp_vertices)
    pos = {v: i for i, v in enumerate(cs)}
    table = {}
    if graph is not None and len(cs) > 1:
        sub = graph.induced_subgraph(cs)
        for e in sub.es:
            a, b = e.source, e.target
            table[(a, b) if a < b else (b, a)] = e['distance']
    return cs, pos, table


def score_kept_set(kept, pos: Dict, edge_table: Dict) -> Tuple[int, float]:
    """
    Score a candidate: (retained connections, -total retained distance).

    A connection is retained only when BOTH its endpoints are kept, which is what
    makes the score non-separable and the joint search meaningful: dropping one
    vertex destroys its partner's edge, so the duplicate sets cannot be decided
    independently. Compared lexicographically - connections first, total distance
    only as a tie-break, which is needed because the count alone ties in roughly
    70% of components and would otherwise be settled by candidate emission order.
    """
    kp = {pos[v] for v in kept}
    inside = [d for (a, b), d in edge_table.items() if a in kp and b in kp]
    return len(inside), -sum(inside)


# ---------------------------------------------------------------------------
# 6.2b Anchor filtering - trim a component to what fits around a centre
# ---------------------------------------------------------------------------

def _anchor_scaled(vertex_id: int) -> Tuple[float, float]:
    """Vertex in pairing's scaled space, used only to RANK candidates by closeness.

    Always normalises, even when pairing itself ran unnormalised: a raw-unit
    euclidean over (m/z in Da, RT in s) is dominated by RT and would make the
    weight ratio meaningless.

    No origin offset. The only consumer is _anchor_distance, which subtracts two of
    these, so any common offset cancels - the axis ranges and the weight ratio are
    the whole of what shapes the ranking.
    """
    cfg = anchor_config
    a = vertex_attrs[vertex_id]
    return (a['x_center'] / cfg['range_x'],
            a['y_center'] / cfg['range_y'] * cfg['mz_rt_weight_ratio'])


def _anchor_distance(u: int, v: int) -> float:
    ux, uy = _anchor_scaled(u)
    vx, vy = _anchor_scaled(v)
    return math.hypot(ux - vx, uy - vy)


def _in_anchor_range(anchor: int, v: int, ppm: float, rt: float) -> bool:
    """Is v inside the anchor's tolerance box?

    The m/z tolerance is derived from the REFERENCE feature's own m/z, exactly as
    pairing does it - and here the anchor is the reference, so the test is
    asymmetric by construction.
    """
    a, b = vertex_attrs[anchor], vertex_attrs[v]
    return (abs(a['x_center'] - b['x_center']) <= a['x_center'] * (ppm / 1e6)
            and abs(a['y_center'] - b['y_center']) <= rt)


def _pair_conforms(u: int, v: int, ppm: float, rt: float) -> bool:
    """Symmetric conformance, for the 'prune' radius mode."""
    a, b = vertex_attrs[u], vertex_attrs[v]
    tol = min(a['x_center'], b['x_center']) * (ppm / 1e6)
    return (abs(a['x_center'] - b['x_center']) <= tol
            and abs(a['y_center'] - b['y_center']) <= rt)


def build_anchor_pool(anchor: int, comp_vertices: List[int],
                      ppm: float, rt: float, radius_mode: str) -> List[int]:
    """
    Vertices this anchor would keep in play. Everything outside is evicted.

    Radius modes trade coverage against pairwise conformance. Being within the
    cutoff OF THE ANCHOR does not put two pool members within the cutoff of EACH
    OTHER - they can sit on opposite sides and end up 2x the cutoff apart:
      full  - radius = the cutoff. Most permissive, no pairwise guarantee.
      half  - radius = cutoff/2. Triangle inequality then guarantees every pair
              conforms, at the cost of evicting substantially more.
      prune - full radius, then drop members that clash with an already-accepted
              one, nearest to the anchor first so the tightest survive.
    """
    if radius_mode == 'half':
        ppm, rt = ppm / 2.0, rt / 2.0
    pool = [v for v in comp_vertices if v == anchor or _in_anchor_range(anchor, v, ppm, rt)]

    if radius_mode == 'prune' and len(pool) > 2:
        accepted = []
        for v in sorted(pool, key=lambda x: (_anchor_distance(anchor, x), x)):
            if all(_pair_conforms(v, u, ppm, rt) for u in accepted):
                accepted.append(v)
        pool = sorted(accepted)
    return pool


def resolve_pool_nearest(pool: List[int], anchor: int) -> Tuple[List[int], float]:
    """One vertex per file: the one closest to the anchor. Score = -total distance."""
    by_file = {}
    for v in pool:
        by_file.setdefault(vertex_attrs[v]['filename'], []).append(v)
    kept = sorted(min(vids, key=lambda v: (_anchor_distance(anchor, v), v))
                  for vids in by_file.values())
    return kept, -sum(_anchor_distance(anchor, v) for v in kept)


def _local_search(seed: List[int], by_file: Dict, pos: Dict, edge_table: Dict,
                  max_rounds: int = 12) -> Tuple[List[int], Tuple]:
    """
    Improve a kept set by swapping one file's choice at a time.

    Used when a pool has too many combinations to enumerate. The alternative -
    the greedy heuristic - ranks vertices by their FULL-GRAPH degree, which counts
    edges to vertices the same candidate is deleting; that is separable, so it
    optimises something other than what survives. Hill-climbing on the real score
    fixes that and can only improve on its seed.

    Scored INCREMENTALLY. Swapping one vertex only disturbs the edges incident to
    the two vertices involved, so each trial costs their degree rather than a scan
    of the whole component - the difference between a component-sized rescan per
    trial and a handful of lookups, which is what makes this affordable on the big
    components that reach it.

    Deterministic: files and their vertices are visited in sorted order and a swap
    is taken only on a strict improvement.
    """
    adj = {}
    for (a, b), d in edge_table.items():
        adj.setdefault(a, {})[b] = d
        adj.setdefault(b, {})[a] = d

    current = set(pos[v] for v in seed)
    n_edges = 0
    sum_d = 0.0
    for p in current:
        for q, d in adj.get(p, {}).items():
            if q in current and q > p:
                n_edges += 1
                sum_d += d

    def delta(p, kept_without):
        """Edges and distance p contributes against the rest of the kept set."""
        n = 0
        s = 0.0
        for q, d in adj.get(p, {}).items():
            if q in kept_without:
                n += 1
                s += d
        return n, s

    held_of = {vertex_attrs[v]['filename']: pos[v] for v in seed}
    for _ in range(max_rounds):
        improved = False
        for filename in sorted(by_file):
            vids = by_file[filename]
            if len(vids) < 2:
                continue
            held = held_of.get(filename)
            if held is None:
                continue
            rest = current - {held}
            hn, hs = delta(held, rest)
            for alt in sorted(vids):
                ap = pos[alt]
                if ap == held:
                    continue
                an, a_s = delta(ap, rest)
                if (an, -a_s) > (hn, -hs):
                    current = rest | {ap}
                    n_edges += an - hn
                    sum_d += a_s - hs
                    held_of[filename] = ap
                    held, rest, hn, hs = ap, current - {ap}, an, a_s
                    improved = True
        if not improved:
            break

    inv = {p: v for v, p in pos.items()}
    return sorted(inv[p] for p in current), (n_edges, -sum_d)


def resolve_pool_connectivity(pool: List[int], comp_idx: int, pos: Dict,
                              edge_table: Dict, mega_complexity_threshold: float,
                              anchor: int = None) -> Tuple[List[int], Tuple]:
    """One vertex per file, chosen jointly by score_kept_set over the pool."""
    by_file = {}
    for v in pool:
        by_file.setdefault(vertex_attrs[v]['filename'], []).append(v)

    candidates, _, use_heuristic = generate_filter_candidates(
        comp_idx, by_file, graph, mega_complexity_threshold)

    if use_heuristic:
        # Too many combinations to enumerate. Seed from nearest-to-anchor - a good
        # starting point, and it makes this method's result at least as connected
        # as anchor_nearest's on the same pool - then hill-climb.
        seed = (resolve_pool_nearest(pool, anchor)[0] if anchor is not None
                else sorted(candidates[0]))
        return _local_search(seed, by_file, pos, edge_table)

    best_key = best_kept = None
    for candidate in candidates:
        kept = sorted(candidate)
        key = score_kept_set(kept, pos, edge_table)
        if best_key is None or key > best_key:
            best_key, best_kept = key, kept
    if best_kept is None:
        best_kept = sorted(pool)
        best_key = score_kept_set(best_kept, pos, edge_table)
    return best_kept, best_key


def trims_components(cfg: Dict = None) -> bool:
    """
    Does the configured method reduce a component beyond de-duplicating it?

    True for every method except duplicate_only. Such methods may evict a vertex
    for sitting too far from the centre even when no file contributed a duplicate,
    so they must be shown EVERY component - not just the ones the duplicate filter
    would have work to do on.
    """
    cfg = cfg if cfg is not None else anchor_config
    return bool(cfg and cfg.get('method', 'duplicate_only') != 'duplicate_only')


# ---------------------------------------------------------------------------
# 6.2c BFS filtering - grow outward from a start vertex, one vertex per file
# ---------------------------------------------------------------------------

def grow_bfs(start: int, adj: Dict, max_depth, ppm_budget, rt_budget) -> Tuple[List[int], float]:
    """
    Breadth-first growth from `start`, claiming at most one vertex per file.

    EVERY vertex occupies its file, intermediate path nodes included, so the
    result is a connected subgraph with no file repeated - a path that "passes
    through" a file has spent it.

    Within a layer candidates are taken in ascending edge distance and claim their
    file as they go, so where several vertices of one file sit in the same layer
    the nearest wins.

    The budgets are on the ACCUMULATED offset along the path from the start, which
    is what makes them useful: by the triangle inequality the direct offset never
    exceeds the accumulated one, so bounding it at B puts every member within B of
    the start and therefore every PAIR within 2B. Setting B to half the pairing
    cutoff makes the finished component pairwise-conforming - measured at exactly
    zero violations, where an unbudgeted graph gives ~14% of pairs breaking.

    Returns (kept_vertices, total_edge_distance).
    """
    tree = {start}
    used_files = {vertex_attrs[start]['filename']}
    path_rt = {start: 0.0}
    path_ppm = {start: 0.0}
    total_d = 0.0
    frontier = [start]
    depth = 0

    while frontier and (max_depth is None or depth < max_depth):
        depth += 1
        candidates = []
        for u in frontier:
            for v, d in adj.get(u, {}).items():
                if v not in tree and vertex_attrs[v]['filename'] not in used_files:
                    candidates.append((d, v, u))
        candidates.sort()

        next_frontier = []
        for d, v, u in candidates:
            # Re-check: an earlier candidate in this same layer may have taken the file.
            if v in tree or vertex_attrs[v]['filename'] in used_files:
                continue
            a, b = vertex_attrs[u], vertex_attrs[v]
            nrt = path_rt[u] + abs(b['y_center'] - a['y_center'])
            if rt_budget is not None and nrt > rt_budget:
                continue
            nppm = path_ppm[u] + (abs(b['x_center'] - a['x_center'])
                                  / min(a['x_center'], b['x_center']) * 1e6)
            if ppm_budget is not None and nppm > ppm_budget:
                continue
            tree.add(v)
            used_files.add(b['filename'])
            path_rt[v] = nrt
            path_ppm[v] = nppm
            total_d += d
            next_frontier.append(v)
        if not next_frontier:
            break
        frontier = next_frontier

    return sorted(tree), total_d


# ---------------------------------------------------------------------------
# 6.2d Window filtering - the exact "every pair inside the cutoffs" formulation
# ---------------------------------------------------------------------------

def component_windows(comp_vertices: List[int], ppm: float, rt: float) -> List[List[int]]:
    """
    Every maximal (m/z, RT) window of a component.

    "Every pair within cutoff" reads like a max-clique constraint, but the
    compatibility relation here is axis-aligned and one-dimensional per axis:

        all pairs conform  <=>  (mz_max - mz_min) <= tol AND (rt_max - rt_min) <= cut

    so it is really a constraint on the BOUNDING BOX of the kept set. Every feasible
    set is therefore contained in one of these windows, which makes enumerating them
    EXACT - no feasible selection is missed - where a radius or budget rule only
    approximates the constraint from the conservative side.

    Generated with two pointers over the m/z order and then the RT order inside each
    m/z span, so windows are contiguous runs rather than sets rebuilt from scratch.
    Sort keys carry the vertex id so the order is total and the output deterministic.
    """
    by_mz = sorted(comp_vertices, key=lambda v: (vertex_attrs[v]['x_center'], v))
    n = len(by_mz)
    seen = set()
    out = []
    j = 0
    for i in range(n):
        lo = vertex_attrs[by_mz[i]]['x_center']
        if j < i:
            j = i
        while (j + 1 < n
               and (vertex_attrs[by_mz[j + 1]]['x_center'] - lo) / lo * 1e6 <= ppm):
            j += 1
        span = sorted(by_mz[i:j + 1], key=lambda v: (vertex_attrs[v]['y_center'], v))
        m = len(span)
        b = 0
        for a in range(m):
            base = vertex_attrs[span[a]]['y_center']
            if b < a:
                b = a
            while b + 1 < m and vertex_attrs[span[b + 1]]['y_center'] - base <= rt:
                b += 1
            w = frozenset(span[a:b + 1])
            if w not in seen:
                seen.add(w)
                out.append(sorted(w))
    return out


def _mean_pair_distance(sel: List[int]) -> float:
    """Mean distance over all pairs of a selection, in the scaled space."""
    n = len(sel)
    if n < 2:
        return 0.0
    total = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            total += _anchor_distance(sel[i], sel[j])
    return total / (n * (n - 1) / 2.0)


def select_tightest_in_window(win: List[int],
                              mega_complexity_threshold: float) -> Tuple[List[int], float]:
    """
    One vertex per file inside a window, minimising MEAN pairwise distance.

    The pairwise cutoff is already satisfied by every member of the window, and the
    file count is fixed by the window, so the only thing left to choose is which
    vertex represents each file. Mean rather than total because total scales with the
    pair count and is only comparable at equal coverage - which is all that is ever
    compared here, so the two rank identically; mean is kept because it is the
    number worth reporting.

    Exhaustive while the combination count is small, hill-climbing beyond that.
    """
    by_file = {}
    for v in win:
        by_file.setdefault(vertex_attrs[v]['filename'], []).append(v)
    files = sorted(by_file)

    combos = 1
    for f in files:
        combos *= len(by_file[f])

    if combos <= mega_complexity_threshold:
        best_sel = best_score = None
        for combo in product(*[sorted(by_file[f]) for f in files]):
            sel = sorted(combo)
            score = _mean_pair_distance(sel)
            if best_score is None or score < best_score:
                best_score, best_sel = score, sel
        return best_sel, best_score

    # Seed each file with the vertex closest to the rest of the window, then swap
    # one file at a time while it strictly improves.
    current = sorted(min(by_file[f],
                         key=lambda v: (sum(_anchor_distance(v, u) for u in win if u != v), v))
                     for f in files)
    best_score = _mean_pair_distance(current)
    for _ in range(8):
        improved = False
        for f in files:
            if len(by_file[f]) < 2:
                continue
            held = next(v for v in current if vertex_attrs[v]['filename'] == f)
            for alt in sorted(by_file[f]):
                if alt == held:
                    continue
                trial = sorted([v for v in current if v != held] + [alt])
                score = _mean_pair_distance(trial)
                if score < best_score:
                    best_score, current, held, improved = score, trial, alt, True
        if not improved:
            break
    return current, best_score


def window_filter_component(comp_idx: int, comp_vertices: List[int], pos: Dict,
                            edge_table: Dict,
                            mega_complexity_threshold: float) -> Tuple[List[int], int, int]:
    """
    Keep the tightest window that covers the most files.

    Ranked by (files covered, -mean pairwise distance, -lowest vertex id). Because
    coverage is the primary key and a window's coverage is fixed by its contents, only
    the windows achieving the MAXIMUM coverage need the expensive inner selection -
    the rest are discarded on a set-length comparison. That matters: window counts are
    heavily skewed (median 11 per component, but up to ~42k on the largest).

    Returns (kept_vertices, representative_vertex, n_tied_windows) - the third value
    keeps the signature aligned with the anchor and BFS filters.
    """
    cfg = anchor_config
    ppm, rt = cfg['window_ppm'], cfg['window_rt']

    wins = component_windows(comp_vertices, ppm, rt)
    if not wins:
        return sorted(comp_vertices), comp_vertices[0], 0

    counts = [len({vertex_attrs[v]['filename'] for v in w}) for w in wins]
    top = max(counts)
    best_windows = [w for w, c in zip(wins, counts) if c == top]

    cap = cfg.get('window_max_candidates') or 0
    if cap > 0 and len(best_windows) > cap:
        # Deterministic truncation: windows are already emitted in sorted order.
        best_windows = best_windows[:cap]

    best_key = best_sel = None
    for w in best_windows:
        sel, score = select_tightest_in_window(w, mega_complexity_threshold)
        key = (-score, -min(sel))
        if best_key is None or key > best_key:
            best_key, best_sel = key, sel
    return best_sel, min(best_sel), len(best_windows)


def bfs_filter_component(comp_idx: int, comp_vertices: List[int], pos: Dict,
                         edge_table: Dict) -> Tuple[List[int], int, int]:
    """
    Try every vertex as the BFS start; keep the best result.

    Ranked by (files covered, -total edge distance, -vertex id). Files first
    because coverage is the primary goal; total distance because among equally
    complete solutions the tightest one is wanted; vertex id last so the choice is
    total and cannot depend on iteration order.

    Note depth is a weak knob in practice: with a budget of half the cutoff, a
    2-hop path has to fit both hops inside it, which almost never happens - going
    from depth 1 to depth 2 moved average linkage by 0.002 on the cycle-3 sample.
    The budgets are what actually control the outcome.

    Returns (kept_vertices, chosen_start, n_tied_starts).
    """
    cfg = anchor_config
    max_depth = cfg.get('bfs_max_depth')
    ppm_budget = cfg.get('bfs_ppm_budget')
    rt_budget = cfg.get('bfs_rt_budget')

    # Adjacency over the component, in original vertex ids, weighted by the
    # stored pairing distance. Built from the de-duplicated edge table so a mutual
    # match is one edge, not two.
    adj = {}
    for (a, b), d in edge_table.items():
        u, w = comp_vertices[a], comp_vertices[b]
        adj.setdefault(u, {})[w] = d
        adj.setdefault(w, {})[u] = d

    best_key = best_kept = best_start = None
    n_tied = 0
    for start in comp_vertices:
        kept, dist = grow_bfs(start, adj, max_depth, ppm_budget, rt_budget)
        key = (len(kept), -dist, -start)
        if best_key is None or key > best_key:
            best_key, best_kept, best_start = key, kept, start
            n_tied = 1
        elif key[:2] == best_key[:2]:
            n_tied += 1
    return best_kept, best_start, n_tied


def anchor_filter_component(comp_idx: int, comp_vertices: List[int], pos: Dict,
                            edge_table: Dict,
                            mega_complexity_threshold: float) -> Tuple[List[int], int, int]:
    """
    Pick an anchor, trim to its pool, resolve duplicates inside it.

    Anchors that tie on file coverage are RESOLVED BY EVALUATION, not by an
    arbitrary key: each is scored by the method in use and the best result wins.
    Ties there fall back to the lowest vertex id, because the module's determinism
    contract (see the header) admits no dependence on iteration order.

    Returns (kept_vertices, chosen_anchor, tied_anchor_count).
    """
    cfg = anchor_config
    ppm, rt, mode = cfg['anchor_ppm'], cfg['anchor_rt_cutoff'], cfg['radius_mode']

    pools = {a: build_anchor_pool(a, comp_vertices, ppm, rt, mode) for a in comp_vertices}
    file_counts = {a: len({vertex_attrs[v]['filename'] for v in p}) for a, p in pools.items()}
    top = max(file_counts.values())
    tied = sorted(a for a in comp_vertices if file_counts[a] == top)

    if cfg['method'] == 'anchor_nearest':
        best = None
        for a in tied:
            kept, score = resolve_pool_nearest(pools[a], a)
            key = (score, -a)
            if best is None or key > best[0]:
                best = (key, kept, a)
        return best[1], best[2], len(tied)

    # Connectivity: the enumeration is the expensive part and depends only on the
    # POOL, so anchors sharing a pool are evaluated once. Measured on cycle 3 this
    # collapses 47809 anchor evaluations to 2769 - one component with 2209 tied
    # anchors resolved to a single distinct pool.
    # Resolving a pool is the expensive step and it depends almost entirely on the
    # POOL, so anchors sharing one are resolved once. Measured on cycle 3 this
    # collapses 47809 anchor evaluations to 2769 - one component with 2209 tied
    # anchors resolves to a single distinct pool. `tied` is sorted, so the anchor
    # that seeds a shared pool is always its lowest id: deterministic, and the same
    # anchor the final tie-break would prefer anyway.
    seen = {}
    best = None
    for a in tied:
        key_pool = frozenset(pools[a])
        cached = seen.get(key_pool)
        if cached is None:
            cached = resolve_pool_connectivity(
                pools[a], comp_idx, pos, edge_table, mega_complexity_threshold, a)
            seen[key_pool] = cached
        kept, score = cached
        key = (score, -a)
        if best is None or key > best[0]:
            best = (key, kept, a)
    return best[1], best[2], len(tied)


def count_filter_candidates(file_vertex_map: Dict) -> int:
    """
    Exact number of keep-one-vertex-per-file combinations for a component.

    Filtering must keep exactly one vertex from each duplicate set, so the candidate
    count is the product of the duplicate-set sizes. This is the single source of
    truth for the heuristic-vs-exhaustive decision: it is used both to schedule the
    work (filter_duplicate_file_vertices) and to make the actual choice
    (generate_filter_candidates), so the two can never disagree.

    Args:
        file_vertex_map: {filename: [vertex_ids]} for one component

    Returns:
        Exact combination count (>= 1). May be a very large int - clamp before
        converting to float.
    """
    total = 1
    for vertex_ids in file_vertex_map.values():
        if len(vertex_ids) > 1:
            total *= len(vertex_ids)
    return total


def generate_filter_candidates(comp_idx: int, file_vertex_map: Dict,
                               graph: ig.Graph = None,
                               mega_complexity_threshold: float = 10000) -> Tuple[Iterable, int, bool]:
    """
    Produce the candidate vertex sets to choose between when filtering a component.

    Filtering must keep exactly ONE vertex per source file. Files that contributed a
    single vertex are therefore kept unconditionally; every file that contributed
    several forms a "duplicate set" we must pick from.

    A component IS the unit here - the workflow is clusterless, so the choice is
    made once per component rather than once per sub-cluster.

    Args:
        comp_idx: Component index (for debugging)
        file_vertex_map: {filename: [vertex_ids]} for this component
        graph: igraph Graph, used to score vertices by degree in the heuristic
        mega_complexity_threshold: at or above this candidate count, use the heuristic

    Returns:
        (candidates, combo_estimate, use_heuristic)
        - candidates: a single-element list (heuristic) or a GENERATOR producing one
          set at a time (exhaustive - millions of candidates must not be materialised)
        - combo_estimate: exact number of combinations
        - use_heuristic: which method produced the candidates
    """
    # Split the component's files into "no choice to make" and "must pick one".
    # Sorted by filename so candidate order - and therefore tie-breaking in the
    # scorer, which keeps the first best-scoring candidate - is deterministic.
    always_keep = set()
    duplicate_sets = []
    for filename in sorted(file_vertex_map.keys()):
        vertex_ids = file_vertex_map[filename]
        if len(vertex_ids) > 1:
            duplicate_sets.append(vertex_ids)
        else:
            always_keep.add(vertex_ids[0])

    combo_estimate = count_filter_candidates(file_vertex_map)
    use_heuristic = combo_estimate >= mega_complexity_threshold

    if use_heuristic:
        # ------------------------------------------------------------------
        # SELECTION METHOD 2 of 2: HEURISTIC (GREEDY PER DUPLICATE SET)
        #
        # Used when enumerating every combination would be too expensive. Each
        # duplicate set is decided on its own: keep the highest-degree vertex.
        # Greedy, so it can miss the combination that is best overall - that is
        # the trade against METHOD 1 (see filter_one_component), and why
        # raising mega_complexity_threshold changes which vertices survive.
        #
        # Plain degree is the right measure: every neighbour of a component
        # member is itself a component member, so there is no "outside" to
        # discount against.
        # ------------------------------------------------------------------
        def _degree(vertex_id):
            # Without the graph every degree is 0, so max() falls through to
            # the first vertex of the set.
            return graph.degree(vertex_id) if graph is not None else 0

        candidate = set(always_keep)
        for vertex_ids in duplicate_sets:
            candidate.add(max(vertex_ids, key=_degree))  # max() keeps the FIRST maximum
        return [candidate], combo_estimate, use_heuristic

    # ----------------------------------------------------------------------
    # Exhaustive: yield one combination at a time. Scoring happens in
    # filter_one_component (SELECTION METHOD 1); candidates are produced on
    # demand so millions of them are never held in memory at once.
    # ----------------------------------------------------------------------
    def _combinations(set_idx, chosen):
        if set_idx == len(duplicate_sets):
            yield always_keep.union(chosen)
            return
        for vertex_id in duplicate_sets[set_idx]:
            yield from _combinations(set_idx + 1, chosen + [vertex_id])

    return _combinations(0, []), combo_estimate, use_heuristic


# ---------------------------------------------------------------------------
# 6.3 Per-component filtering
# ---------------------------------------------------------------------------

def filter_one_component(component_data: Dict) -> Dict:
    """
    Decide which vertices ONE component keeps, and which it deletes.

    This is the unit of work. It runs unchanged whether the caller is a Pool
    worker or the in-process path, and never raises: any failure comes back as a
    result with status 'error' so one bad component cannot take the run down.

    Args:
        component_data: Dictionary containing:
            - comp_idx: Component index
            - comp_vertices: List of vertex IDs in component
            - vertex_attrs: Vertex attributes (subset needed for this component)
            - mega_complexity_threshold: Threshold for using heuristic
            NOTE: the graph is read from the module global (SECTION 2), not passed
            in - it is far too large to ship per component.

    Returns:
        Result dict: comp_idx, status ('skipped' | 'processed' | 'error'),
        kept_component ({'vertices': [...]}), deleted_vertices, combo_estimate,
        plus __-prefixed timing/attribution fields consumed by the orchestrator.
    """
    # Set before the try block: the except handler reports elapsed time, so it must
    # be bound even if the failure happens in the very first statements.
    start_time = time.time()
    global _components_seen
    try:
        worker_pid = os.getpid()

        # Log this component only on the sampling cadence (see LOG_EVERY_N_COMPONENTS).
        _components_seen += 1
        verbose = (_components_seen == 1 or _components_seen % LOG_EVERY_N_COMPONENTS == 0)

        # Affinity at start, to verify pinning held for the duration of the task.
        try:
            cpu_list_start = sorted(os.sched_getaffinity(0))
        except Exception:
            cpu_list_start = None

        comp_idx = component_data['comp_idx']
        comp_vertices = component_data['comp_vertices']
        vertex_attrs = component_data['vertex_attrs']
        mega_complexity_threshold = component_data.get('mega_complexity_threshold', 10000)

        # --- Early exits: nothing to choose between -------------------------
        # A lone vertex cannot violate the one-per-file rule.
        if len(comp_vertices) <= 1:
            elapsed = time.time() - start_time
            return {
                'comp_idx': comp_idx,
                'status': 'skipped',
                'kept_component': {'vertices': comp_vertices},
                'deleted_vertices': [],
                'combo_estimate': 0,
                '__worker_pid': worker_pid,
                '__elapsed_ms': int(elapsed * 1000)
            }

        # Group this component's vertices by source file. Clusterless: a component
        # IS the unit, so there is no cluster layer to iterate over.
        file_vertex_map = {}
        for vertex_id in comp_vertices:
            file_vertex_map.setdefault(vertex_attrs[vertex_id]['filename'], []).append(vertex_id)

        # Already one-per-file: the component passes through untouched. This is the
        # common case - the orchestrator screens most of them out before dispatch,
        # but the check is repeated here so the function is correct standalone.
        #
        # The anchor methods are the exception: they may evict a vertex for sitting
        # outside the anchor's tolerance even when no file contributed a duplicate,
        # so they must see every component.
        if (not any(len(vids) > 1 for vids in file_vertex_map.values())
                and not trims_components()):
            return {
                'comp_idx': comp_idx,
                'status': 'skipped',
                'kept_component': {'vertices': comp_vertices},
                'deleted_vertices': [],
                'combo_estimate': 0
            }

        # --- Choose which vertices to keep ----------------------------------
        # Inherited through fork (SECTION 2), not passed in.
        graph = globals().get('graph', None)

        # --- ANCHOR METHODS -------------------------------------------------
        # A different contract from the duplicate filter: instead of only choosing
        # among duplicates, pick a centre and keep what fits around it. Vertices
        # outside the anchor's tolerance are evicted, which is what lets a
        # component shrink to something that actually respects the pairing cutoff
        # pairwise - a duplicate-only filter cannot do that, because the vertices
        # driving most violations come from files that contributed just one.
        if trims_components():
            cs, pos, edge_table = component_edge_table(comp_vertices)
            if anchor_config['method'] == 'bfs':
                kept_list, chosen_anchor, n_tied = bfs_filter_component(
                    comp_idx, cs, pos, edge_table)
            elif anchor_config['method'] == 'window':
                kept_list, chosen_anchor, n_tied = window_filter_component(
                    comp_idx, cs, pos, edge_table, mega_complexity_threshold)
            else:
                kept_list, chosen_anchor, n_tied = anchor_filter_component(
                    comp_idx, cs, pos, edge_table, mega_complexity_threshold)
            best_candidate = set(kept_list)
            retained_edges, neg_distance = score_kept_set(kept_list, pos, edge_table)

            deleted_vertices_list = []
            for vertex_id in comp_vertices:
                if vertex_id in best_candidate:
                    continue
                same_file = file_vertex_map[vertex_attrs[vertex_id]['filename']]
                _m = anchor_config['method']
                reason_text = (f'{_m} - duplicate file vertex' if len(same_file) > 1
                               else f'{_m} - outside tolerance / unreachable')
                deleted_vertices_list.append({
                    'component_id': comp_idx,
                    'vertex_id': vertex_id,
                    'filename': vertex_attrs[vertex_id]['filename'],
                    'feature_idx': vertex_attrs[vertex_id]['feature_idx'],
                    'label': vertex_attrs[vertex_id]['label'],
                    'pep_ident': vertex_attrs[vertex_id].get('pep_ident', []),
                    'prot_ident': vertex_attrs[vertex_id].get('prot_ident', []),
                    'x_center': vertex_attrs[vertex_id].get('x_center', 0.0),
                    'y_center': vertex_attrs[vertex_id].get('y_center', 0.0),
                    'intensity': vertex_attrs[vertex_id].get('intensity', 0.0),
                    'openms_fid': vertex_attrs[vertex_id].get('openms_fid', ''),
                    'ms2_scans': vertex_attrs[vertex_id].get('ms2_scans', ''),
                    'charge': vertex_attrs[vertex_id].get('charge', 0),
                    'reason': reason_text,
                    'anchor_vertex': chosen_anchor,
                    'retained_edges': retained_edges,
                    'total_distance': -neg_distance,
                })

            if verbose:
                print(f"[COMP {comp_idx}] start v{chosen_anchor} ({n_tied} tied): "
                      f"{len(comp_vertices)} -> {len(kept_list)} vertices, "
                      f"{retained_edges} connections", flush=True)

            elapsed = time.time() - start_time
            return {
                'comp_idx': comp_idx,
                'status': 'processed',
                'kept_component': {'vertices': sorted(best_candidate)},
                'deleted_vertices': deleted_vertices_list,
                'combo_estimate': 0,
                '__worker_pid': worker_pid,
                '__elapsed_ms': int(elapsed * 1000),
                '__cpu_start': cpu_list_start,
            }

        if graph is None:
            # Always logged: this is an anomaly, not progress. Without the graph
            # every degree is 0 and both methods degenerate to "keep the first".
            print(f"[COMP {comp_idx}] ⚠️  WARNING: graph is None! Heuristic and induced_subgraph will fall back", flush=True)

        candidates, combo_estimate, use_heuristic = generate_filter_candidates(
            comp_idx, file_vertex_map, graph, mega_complexity_threshold)

        if use_heuristic:
            # METHOD 2 already made the choice inside generate_filter_candidates -
            # there is exactly one candidate and nothing left to score.
            if verbose:
                print(f"[COMP {comp_idx}] heuristic (complexity {combo_estimate:.2e} >= {mega_complexity_threshold:.2e})", flush=True)
            best_candidate = candidates[0] if candidates else set(comp_vertices)
            best_score = None
            candidates_evaluated = 0
            deleted_vertices_list = []
        else:
            best_candidate = None
            best_score = None
            candidates_evaluated = 0
            deleted_vertices_list = []

            if verbose:
                print(f"[COMP {comp_idx}] exhaustive (complexity {combo_estimate:.2e})", flush=True)

            # ----------------------------------------------------------------
            # SELECTION METHOD 1 of 2: EXHAUSTIVE (JOINT OPTIMUM)
            #
            # Used when the candidate count is < mega_complexity_threshold. Every
            # "keep exactly one vertex per file" combination is generated and
            # scored, and the best-scoring one wins. Candidates arrive from a
            # generator, so only one is ever held in memory at a time.
            #
            # The score is what the candidate LEAVES CONNECTED:
            #
            #     score = ( connections with both endpoints kept , -total distance )
            #
            # compared lexicographically. Both parts matter:
            #
            #   Primary. A connection survives only if BOTH its features are kept,
            #   so dropping one vertex destroys its partner's edge. That coupling
            #   is what makes the duplicate sets a JOINT decision and gives the
            #   enumeration something to find - a locally worse vertex can be kept
            #   because it keeps the overall set better connected. Scoring each
            #   vertex on its own full-graph degree instead would be separable, and
            #   a separable score reduces exactly to the greedy heuristic.
            #
            #   Secondary. The count alone ties in roughly 70% of components, which
            #   without a tie-break would leave those decided by candidate emission
            #   order, i.e. by filename sort order. Total retained distance settles
            #   all but a handful.
            #
            # Ties keep the FIRST candidate seen (strict `>` below). Candidate
            # order is deterministic - filenames are iterated sorted - so the
            # outcome is reproducible.
            # ----------------------------------------------------------------
            _, _pos, _edge_table = component_edge_table(comp_vertices)

            for candidate in candidates:
                candidates_evaluated += 1

                # Long exhaustive runs would otherwise look hung
                if candidates_evaluated % 100000 == 0:
                    print(f"[COMP {comp_idx}]   evaluated {candidates_evaluated} candidates, best={best_score}", flush=True)

                kept_vertices = [v for v in comp_vertices if v in candidate]
                if len(kept_vertices) <= 1:
                    continue

                if graph is not None:
                    score = score_kept_set(kept_vertices, _pos, _edge_table)
                else:
                    # Without the graph there are no connections to count; fall
                    # back to candidate size so the comparison stays well-defined.
                    score = (len(kept_vertices), 0.0)

                if best_score is None or score > best_score:
                    best_score = score
                    best_candidate = candidate

                # The losing candidate is garbage-collected on the next iteration,
                # which is what keeps peak memory flat over millions of them.

            if verbose:
                print(f"[COMP {comp_idx}]   done: {candidates_evaluated} candidates evaluated, best={best_score}", flush=True)

        # --- Record the outcome ---------------------------------------------
        # Deleted vertices (same for both methods).
        #
        # Only two keys are actually consumed: 'vertex_id' (to look the vertex up in
        # vertex_attrs on export) and 'component_id' (the sort key that makes the
        # deleted list deterministic). The rest is DIAGNOSTIC bookkeeping, kept
        # because it makes a deleted vertex readable in a debugger / crash dump
        # without a second lookup.
        #
        # CAVEAT: these diagnostic fields never reach disk. deleted_vertices.json is
        # written by export_deleted_vertices_for_next_cycle(), which takes only
        # 'vertex_id' and re-reads every attribute from vertex_attrs. So 'reason' and
        # 'best_score' - the two fields that exist ONLY here - are currently
        # discarded. Surface them in the export if they are wanted for analysis.
        if best_candidate:
            for vertex_id in comp_vertices:
                if vertex_id not in best_candidate:
                    reason_text = ('Duplicate file vertex - filtered via heuristic' if use_heuristic
                                   else 'Duplicate file vertex - filtered for best degree score')
                    deleted_vertices_list.append({
                        'component_id': comp_idx,
                        'vertex_id': vertex_id,
                        'filename': vertex_attrs[vertex_id]['filename'],
                        'feature_idx': vertex_attrs[vertex_id]['feature_idx'],
                        'label': vertex_attrs[vertex_id]['label'],
                        'pep_ident': vertex_attrs[vertex_id].get('pep_ident', []),
                        'prot_ident': vertex_attrs[vertex_id].get('prot_ident', []),
                        'x_center': vertex_attrs[vertex_id].get('x_center', 0.0),
                        'y_center': vertex_attrs[vertex_id].get('y_center', 0.0),
                        'intensity': vertex_attrs[vertex_id].get('intensity', 0.0),
                        'openms_fid': vertex_attrs[vertex_id].get('openms_fid', ''),
                        'ms2_scans': vertex_attrs[vertex_id].get('ms2_scans', ''),
                        'charge': vertex_attrs[vertex_id].get('charge', 0),
                        'reason': reason_text,
                        # Split rather than a single number: the exhaustive score is
                        # a (connections, -distance) pair and float() would raise.
                        'retained_edges': best_score[0] if best_score is not None else None,
                        'total_distance': -best_score[1] if best_score is not None else None,
                    })

            # sorted(), not list(): a set's iteration order depends on how it
            # was built, so listing it directly would make the emitted vertex
            # order track the candidate construction sequence. Sorting pins it to
            # ascending vertex id - same reasoning as the component ordering in
            # export_accepted_components.
            kept_component = {'vertices': sorted(best_candidate)}
        else:
            # No candidate scored (e.g. every combination collapsed to <=1 kept
            # vertex): leave the component exactly as it came in.
            kept_component = {'vertices': comp_vertices}

        elapsed = time.time() - start_time

        # Verify CPU pinning held for the whole task. A CHANGE is always reported;
        # the "held" confirmation only on the sampling cadence.
        try:
            cpu_list_end = sorted(os.sched_getaffinity(0))
            if cpu_list_start and cpu_list_end != cpu_list_start:
                print(f"[WORKER {worker_pid}] ⚠️  AFFINITY CHANGED during task: "
                      f"{cpu_list_start} -> {cpu_list_end}", file=sys.stderr, flush=True)
            elif verbose:
                print(f"[WORKER {worker_pid}] affinity held on {cpu_list_end} "
                      f"(component #{_components_seen}, {elapsed:.2f}s)", file=sys.stderr, flush=True)
        except Exception:
            pass

        return {
            'comp_idx': comp_idx,
            'status': 'processed',
            'kept_component': kept_component,
            'deleted_vertices': deleted_vertices_list,
            'combo_estimate': combo_estimate,
            '__worker_pid': worker_pid,
            '__elapsed_ms': int(elapsed * 1000)
        }

    except Exception as e:
        # Never propagate: the orchestrator keeps the unfiltered component and
        # carries on with the rest of the cycle.
        elapsed = time.time() - start_time
        return {
            'comp_idx': component_data.get('comp_idx', -1),
            'status': 'error',
            'error': str(e),
            'kept_component': {'vertices': component_data.get('comp_vertices', [])},
            'deleted_vertices': [],
            'combo_estimate': -1,
            '__worker_pid': os.getpid(),
            '__elapsed_ms': int(elapsed * 1000)
        }


# ---------------------------------------------------------------------------
# 6.4 Worker chunk driver
# ---------------------------------------------------------------------------

def process_component_chunk(metadata_stream):
    """
    Process one worker's entire share of components, in memory-paced batches.

    This is the function handed to Pool.imap_unordered (and called directly on the
    in-process path). What crosses the process boundary is only LIGHTWEIGHT
    METADATA - the heavy per-component work items are built here, just in time,
    a batch at a time, and dropped again once that batch is done. Without this the
    orchestrator would have to hold every work item at once.

    Batch boundary: accumulate components until their summed load weight would
    exceed COMBINATIONS_BATCH_LIMIT, then build, process, free and gc.collect().

    Args:
        metadata_stream: iterator of metadata dicts for this worker's components.

    Returns:
        List of result dicts from filter_one_component, tagged with the worker pid
        and batch number.
    """
    # Read from the global inherited via fork (SECTION 2); freed at the end.
    global vertex_attrs

    worker_start = time.time()
    worker_pid = os.getpid()

    def build_work_item(metadata):
        """Expand lightweight metadata into a full work item for filter_one_component."""
        comp_idx = metadata['comp_idx']
        comp_vertices = metadata['comp_vertices']

        # Component-scoped slice of the shared attrs, so the worker function never
        # needs the full graph-wide dict.
        comp_vertex_attrs = {v: vertex_attrs[v] for v in comp_vertices}

        # Rough size of this work item, purely to flag pathological components in
        # the log. Per-entry costs are hand-estimated, not measured.
        comp_vertex_count = len(comp_vertices)
        estimated_item_size_bytes = (
            comp_vertex_count * 8 +      # comp_vertices list
            comp_vertex_count * 100      # comp_vertex_attrs dict
        )

        if estimated_item_size_bytes > 50_000_000:  # >50MB
            print(f"[PID {worker_pid}] ⚠ LARGE COMPONENT DETECTED: comp_idx={comp_idx}, "
                  f"vertices={comp_vertex_count}, "
                  f"est_size={estimated_item_size_bytes / 1_000_000:.1f}MB", flush=True)

        return {
            'comp_idx': comp_idx,
            'comp_vertices': comp_vertices,
            'vertex_attrs': comp_vertex_attrs,
            'mega_complexity_threshold': metadata.get('mega_complexity_threshold', 10000)
        }

    def run_batch(batch_metadata, batch_num):
        """Build, run and free one batch; returns its tagged results."""
        batch_start = time.time()
        batch_work_items = [build_work_item(m) for m in batch_metadata]
        batch_results = [filter_one_component(item) for item in batch_work_items]
        batch_elapsed = time.time() - batch_start

        for result in batch_results:
            result['__worker_pid'] = worker_pid
            result['__batch_num'] = batch_num
            result['__batch_elapsed_ms'] = int(batch_elapsed * 1000)

        del batch_work_items
        gc.collect()
        return batch_results, batch_elapsed

    metadata_list = list(metadata_stream)

    results = []
    batch_times = []
    batch_count = 0
    total_metadata_processed = 0
    current_batch = []
    current_weight = 0
    processing_start = time.time()

    for metadata in metadata_list:
        load_weight = metadata.get('load_weight', 0)

        # Would adding this component push the batch over the memory-pacing limit?
        if current_weight + load_weight > COMBINATIONS_BATCH_LIMIT and current_batch:
            batch_count += 1
            batch_results, batch_elapsed = run_batch(current_batch, batch_count)
            results.extend(batch_results)
            batch_times.append(batch_elapsed)
            current_batch = []
            current_weight = 0

        current_batch.append(metadata)
        current_weight += load_weight
        total_metadata_processed += 1

    if current_batch:
        batch_count += 1
        batch_results, batch_elapsed = run_batch(current_batch, batch_count)
        results.extend(batch_results)
        batch_times.append(batch_elapsed)

    processing_end = time.time()

    # Safe to drop the shared attrs now: build_work_item will not run again.
    try:
        del vertex_attrs
    except:
        pass
    gc.collect()

    # One line per worker: component count, batch count, and the split between
    # time spent filtering and total time alive.
    proc = processing_end - processing_start
    total = processing_end - worker_start
    slowest = f", slowest batch {max(batch_times):.2f}s" if batch_times else ""
    print(f"[WORKER {worker_pid}] {total_metadata_processed} components in {batch_count} "
          f"batch(es): {proc:.2f}s processing, {total:.2f}s total{slowest}", flush=True)

    return results


# ---------------------------------------------------------------------------
# 6.5 Orchestrator
# ---------------------------------------------------------------------------

def filter_duplicate_file_vertices(graph: ig.Graph, components, vertex_attrs: Dict,
                                   num_workers: int = None,
                                   use_multiprocessing: bool = True,
                                   mega_complexity_threshold: float = 10000,
                                   anchor_cfg: Dict = None) -> Tuple[Dict, List]:
    """
    Enforce the one-feature-per-file rule across every component.

    Four phases:
      PHASE 1    scan every component, pass through the ones that already satisfy
                 the rule, and collect the ones that do not
      PHASE 1.5  build lightweight metadata for those, including a load weight
      PHASE 2    dispatch - either in-process or across a pinned worker Pool
      PHASE 3    merge results back and report

    Args:
        graph: igraph Graph object
        components: the connected components, as returned by graph.components().
            Indexed by component id; components[i] is already a plain list of
            vertex ids, so nothing is copied out of it here.
        vertex_attrs: {vertex_id: {'filename': str, 'feature_idx': int, ...}}
        num_workers: Worker processes for parallel filtering (see resolve_worker_count).
            None -> DEFAULT_NUM_WORKERS (1 = in-process). Does not affect results, only runtime.
        use_multiprocessing: False (--disable-multiprocessing) forces the in-process path.
        mega_complexity_threshold: Candidate-combination count at or above which a component
            is filtered by the greedy heuristic instead of exhaustive evaluation.
            Set in production by params.mmf_mega_complexity_threshold in
            map_mzml_features_v2_clusterless.nf, threaded through the orchestrator.
            Raising it sends more components down the exhaustive path, which is
            slower and CAN CHANGE WHICH VERTICES ARE KEPT: the exhaustive path
            optimises the retained set jointly, the heuristic decides greedily.

    Returns:
        Tuple of (accepted_components, deleted_vertices_list) where
        accepted_components is {comp_idx: {'vertices': [...]}} covering EVERY
        component - passed through untouched or with duplicates removed.
    """
    # Publish the shared structures as module globals so forked workers inherit
    # them copy-on-write instead of receiving pickled copies. Must happen before
    # the Pool is created. See SECTION 2.
    current_module = sys.modules[__name__]
    current_module.graph = graph
    current_module.vertex_attrs = vertex_attrs
    current_module.anchor_config = anchor_cfg

    using_anchor = trims_components(anchor_cfg)

    print("\n" + "="*70)
    if using_anchor and anchor_cfg['method'] == 'window':
        print(f"FILTERING - window ({anchor_cfg['window_ppm']} ppm x "
              f"{anchor_cfg['window_rt']}s, set-wise; every pair conforms by construction)")
    elif using_anchor and anchor_cfg['method'] == 'bfs':
        print(f"FILTERING - bfs (depth {anchor_cfg['bfs_max_depth'] or 'unlimited'}, "
              f"budgets {anchor_cfg['bfs_ppm_budget']} ppm / "
              f"{anchor_cfg['bfs_rt_budget']}s accumulated)")
    elif using_anchor:
        print(f"FILTERING - {anchor_cfg['method']} "
              f"(radius {anchor_cfg['radius_mode']}, "
              f"{anchor_cfg['anchor_ppm']} ppm / {anchor_cfg['anchor_rt_cutoff']}s)")
    else:
        print("FILTERING DUPLICATE FILE VERTICES (Parallel Processing)")
    print("="*70)

    num_workers, cpu_list = resolve_worker_count(num_workers)

    deleted_vertices_list = []
    accepted_components = {}    # comp_idx -> {'vertices': [...]}, filled below

    total_components = len(components)
    components_with_dups = 0
    components_metadata = []

    # ---- PHASE 1: find the components that actually need filtering ---------
    # Everything else is written straight into accepted_components and never
    # dispatched. Only indices and small maps are retained here - building full
    # work items for every component would exhaust RAM on large cycles.
    print(f"PHASE 1: Identifying components with duplicate files ({total_components} total)...", flush=True)
    sys.stdout.flush()

    components_needing_filter = []  # (comp_idx, file_vertex_map, comp_vertices)

    for comp_idx, comp_vertices in enumerate(components):
        if comp_idx and comp_idx % 5000 == 0:
            print(f"  scanned {comp_idx}/{total_components} components...", flush=True)

        # A lone vertex cannot violate the rule
        if len(comp_vertices) <= 1:
            accepted_components[comp_idx] = {'vertices': comp_vertices}
            continue

        # Group by source file; a component needs filtering only if some file
        # contributed more than one vertex.
        file_vertex_map = {}
        for vertex_id in comp_vertices:
            file_vertex_map.setdefault(vertex_attrs[vertex_id]['filename'], []).append(vertex_id)

        # accepts components without conflicts if not in anchor mode and adds them to accepted components
        # onyl fires if anchor mode is off, since otherwsie we have to always look into each component, because the decision isn't only about duplicates any more
        if not any(len(vids) > 1 for vids in file_vertex_map.values()) and not using_anchor:
            accepted_components[comp_idx] = {'vertices': comp_vertices}
            continue

        components_with_dups += 1
        components_needing_filter.append((comp_idx, file_vertex_map, comp_vertices))

    # ---- PHASE 1.5: metadata + load weights --------------------------------
    # Workers get only this; they expand it into work items themselves.
    if components_needing_filter:
        print(f"PHASE 1.5: Preparing metadata for {len(components_needing_filter)} components with duplicates...", flush=True)

        for comp_idx, file_vertex_map, comp_vertices in components_needing_filter:
            # Exact candidate count for this component - the same value
            # generate_filter_candidates will compute in the worker, so
            # will_use_heuristic below is a guarantee, not a guess.
            combo_count = count_filter_candidates(file_vertex_map)

            # load_weight drives BOTH the greedy bin-packing across workers and
            # the batch boundary inside a worker. It is a cost estimate, not a
            # candidate count - the two diverge for heuristic components.
            #
            # NOTE: compare the UNCLAMPED count - clamping first would misclassify
            # a huge component whenever the threshold is set above the clamp.
            if using_anchor:
                # Anchor cost is dominated by building one pool per vertex, which
                # is quadratic in component size, plus at most one enumeration per
                # DISTINCT pool. Size is the usable predictor; the combination
                # count is not, because most anchor components have no duplicates
                # at all and would all weigh 1.
                load_weight = min(float(len(comp_vertices)) ** 2, 10 ** 6)
                will_use_heuristic = combo_count >= mega_complexity_threshold
            elif combo_count >= mega_complexity_threshold:
                # The heuristic is a greedy pass over the component's vertices,
                # NOT an exponential search, so weight it by size rather than by
                # combination count.
                load_weight = min(float(len(comp_vertices)) * 10, 1000)
                will_use_heuristic = True
            else:
                # Exhaustive: cost really is the number of combinations. Clamp as
                # an int first - a huge product would overflow float().
                load_weight = float(min(combo_count, 10 ** 15))
                will_use_heuristic = False

            components_metadata.append({
                'comp_idx': comp_idx,
                'load_weight': load_weight,
                'will_use_heuristic': will_use_heuristic,  # for the summary line below
                'comp_vertices': comp_vertices,
                'mega_complexity_threshold': mega_complexity_threshold
            })

        num_heuristic = sum(1 for m in components_metadata if m['will_use_heuristic'])
        weights = [m['load_weight'] for m in components_metadata]
        if weights:
            print(f"  {num_heuristic} heuristic / {len(components_metadata) - num_heuristic} exhaustive; "
                  f"candidate counts min {min(weights):.1e}, max {max(weights):.1e}, "
                  f"avg {sum(weights)/len(weights):.1e}", flush=True)

    # Per-process attribution, filled from the __-prefixed result fields
    components_per_pid = {}   # pid -> component count
    component_times_ms = []

    filter_start_time = time.time()

    # ---- PHASE 2: dispatch --------------------------------------------------
    if components_metadata:

        results = []
        filter_timing = {}
        sorted_metadata = []
        metadata_chunks = []

        # A single worker has nothing to overlap, so skip the Pool entirely and run
        # in-process - fork + IPC would be pure overhead (measurably slower).
        if not use_multiprocessing or num_workers <= 1:
            reason = "--disable-multiprocessing" if not use_multiprocessing else "num_workers=1"
            print(f"Processing {len(components_metadata)} components in-process ({reason})...\n", flush=True)
            results = process_component_chunk(iter(components_metadata))
        else:
            # Split the work across workers by GREEDY BIN-PACKING on load_weight:
            # sort heaviest-first, then repeatedly give the next component to the
            # least-loaded chunk. Component count alone would be a poor split -
            # one exhaustive component can outweigh thousands of heuristic ones.
            sorted_metadata = sorted(components_metadata, key=lambda m: m.get('load_weight', 0), reverse=True)

            chunk_weights = [0] * num_workers
            metadata_chunks = [[] for _ in range(num_workers)]

            for metadata in sorted_metadata:
                min_idx = chunk_weights.index(min(chunk_weights))
                metadata_chunks[min_idx].append(metadata)
                chunk_weights[min_idx] += metadata.get('load_weight', 0)

            print(f"  split across {num_workers} workers: "
                  + ", ".join(str(len(c)) for c in metadata_chunks)
                  + f" components (complexity {min(chunk_weights):.1e}-{max(chunk_weights):.1e})", flush=True)
            sys.stdout.flush()

            filter_timing['pool_start'] = time.time()

            try:
                parent_affinity = os.sched_getaffinity(0)
                parent_pid = os.getpid()
                print(f"  parent PID {parent_pid} on CPUs {sorted(parent_affinity)}", flush=True)

                pool = None
                filter_timing['pool_create_start'] = time.time()
                try:
                    # Hand each worker a distinct index so _init_worker can pin them
                    # to distinct CPUs. A shared counter under a lock is the simplest
                    # way to do that - Pool gives the initializer no worker identity.
                    worker_counter = Value('i', 0)
                    counter_lock = Lock()

                    def init_worker_wrapper(counter, lock, cpus):
                        with lock:
                            worker_index = counter.value
                            counter.value += 1
                        _init_worker(worker_index, cpus)

                    # No data in initargs beyond the counter and CPU list: the large
                    # structures come through fork (SECTION 2), saving ~2GB of copies.
                    pool = Pool(processes=num_workers,
                                initializer=init_worker_wrapper,
                                initargs=(worker_counter, counter_lock, cpu_list))
                    filter_timing['pool_create_end'] = time.time()

                    # One task per worker (chunksize=1 over exactly num_workers
                    # chunks), consumed as it completes so results are not all held
                    # in flight at once.
                    filter_timing['result_collection_start'] = time.time()

                    results = []
                    results_received = 0

                    for chunk_results in pool.imap_unordered(
                        process_component_chunk,
                        [iter(chunk) for chunk in metadata_chunks],
                        chunksize=1
                    ):
                        for result in chunk_results:
                            results.append(result)
                            results_received += 1

                            if results_received % 1000 == 0:
                                print(f"[MAIN] Streamed {results_received} results from workers...", flush=True)

                    filter_timing['result_collection_end'] = time.time()

                finally:
                    if pool:
                        pool.close()
                        pool.join()
                        filter_timing['pool_end'] = time.time()

            except Exception as e:
                # Pool creation or a worker died: redo the whole stage in-process
                # rather than lose the cycle.
                import traceback
                print(f"⚠ Multiprocessing failed ({type(e).__name__}: {e})", flush=True)
                print(f"Traceback: {traceback.format_exc()}", flush=True)
                print("Falling back to serial processing...", flush=True)
                results = process_component_chunk(iter(components_metadata))
                filter_timing['pool_end'] = time.time()

        # ---- PHASE 3: merge results back ---------------------------------
        filter_timing['result_processing_start'] = time.time()

        del components_metadata
        del sorted_metadata
        del metadata_chunks
        gc.collect()

        for result in results:
            worker_pid = result.get('__worker_pid')
            elapsed_ms = result.get('__elapsed_ms', 0)

            if worker_pid:
                components_per_pid[worker_pid] = components_per_pid.get(worker_pid, 0) + 1
            if elapsed_ms:
                component_times_ms.append(elapsed_ms)

            if result.get('status') == 'error':
                print(f"⚠ Component {result['comp_idx']} processing error: {result.get('error', 'Unknown')}", flush=True)
                # Keep the component unfiltered rather than dropping it
                accepted_components[result['comp_idx']] = result['kept_component']
                continue

            accepted_components[result['comp_idx']] = result['kept_component']
            if result.get('status') != 'skipped':
                deleted_vertices_list.extend(result['deleted_vertices'])

        filter_timing['result_processing_end'] = time.time()

        del results
        gc.collect()

    # ---- Report -------------------------------------------------------------
    # Wall clock, where that time went (fork / workers / collection), how the
    # components split across processes, and per-component averages.
    if components_per_pid:
        elapsed_time = time.time() - filter_start_time
        print(f"\n{'-'*70}")
        print(f"Filtering ran on {len(components_per_pid)} process(es) in {elapsed_time:.2f}s wall clock")

        if filter_timing and 'pool_start' in filter_timing:
            def _span(a, b):
                return filter_timing.get(b, 0) - filter_timing.get(a, 0) if a in filter_timing else 0.0
            pool_create = _span('pool_create_start', 'pool_create_end')
            collection = _span('result_collection_start', 'result_collection_end')
            post = _span('result_processing_start', 'result_processing_end')
            print(f"  pool fork {pool_create:.2f}s | workers {collection:.2f}s | "
                  f"collect {post:.2f}s")

        # Per-process component counts, so an unbalanced split is visible at a glance.
        counts = sorted(components_per_pid.items())
        print("  per process: " + " | ".join(f"{pid}:{n}" for pid, n in counts))

        if component_times_ms:
            import statistics
            avg_time = statistics.mean(component_times_ms)
            max_time = max(component_times_ms)
            print(f"  per component: avg {avg_time:.0f}ms, max {max_time:.0f}ms")
            if max_time > avg_time * 5:
                print(f"  ⚠ load imbalance: slowest component is {max_time/avg_time:.0f}x the average")
        print(f"{'-'*70}")

    print(f"\n{'='*70}")
    print("FILTERING SUMMARY")
    print(f"{'='*70}")
    print(f"Components processed: {components_with_dups}")
    print(f"Total vertices deleted: {len(deleted_vertices_list)}")
    print(f"{'='*70}\n")

    # Deterministic order for the JSON export - results arrive in worker
    # completion order, which is a race.
    deleted_vertices_list.sort(key=lambda x: (x['component_id'], x['vertex_id']))

    return accepted_components, deleted_vertices_list


# ============================================================================
# SECTION 7 - OUTPUT
#
# Two files matter and are both consumed by orchestrate_clusterless_cycles.py:
#   accepted_components.json  -> feature groups finalised this cycle
#   deleted_vertices.json     -> vertices handed back for re-pairing next cycle
# The graph exports below it are optional debugging artefacts.
# ============================================================================

def export_accepted_components(
    accepted_components: Dict,
    vertex_attrs: Dict,
    output_path: str,
    cycle_num: int
):
    """
    Export components that passed filtering (accepted for final results).

    These components have ORIGINAL coordinates (not aligned or weighted).

    Components are emitted in ascending component_id order. accepted_components is
    filled in worker-completion order when filtering runs in parallel, so
    iterating it directly would make the output order race-dependent - and the
    orchestrator assigns final component IDs by position in this array.

    Args:
        accepted_components: {comp_idx: {'vertices': [vertex_id, ...]}}
        vertex_attrs: Dictionary with vertex attributes (original coordinates)
        output_path: Path to save accepted components JSON
        cycle_num: Current cycle number
    """
    output = {
        'cycle': cycle_num,
        'components': []
    }

    for comp_id, comp_data in sorted(accepted_components.items(), key=lambda item: item[0]):
        # NOTE: 'pep_ident'/'prot_ident' are component-level fields that nothing
        # in this script ever sets, so they always serialise as null. They are
        # emitted because the consumer expects the key to exist; the real
        # identifications live on each vertex below.
        component = {
            'component_id': str(comp_id),
            'num_vertices': len(comp_data.get('vertices', [])),
            'pep_ident': comp_data.get('pep_ident'),
            'prot_ident': comp_data.get('prot_ident'),
            'vertices': []
        }

        # Full attribute record per vertex, with its id folded in
        for vertex_id in comp_data.get('vertices', []):
            if vertex_id in vertex_attrs:
                v = vertex_attrs[vertex_id].copy()
                v['vertex_id'] = vertex_id
                component['vertices'].append(v)

        output['components'].append(component)

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n✓ Exported {len(output['components'])} accepted components to: {os.path.basename(output_path)}")
    print(f"  Total vertices in accepted components: {sum(c['num_vertices'] for c in output['components'])}")


def export_deleted_vertices_for_next_cycle(
    deleted_list: List[Dict],
    vertex_attrs: Dict,
    output_path: str
):
    """
    Export deleted vertices with ORIGINAL coordinates for re-pairing in next cycle.

    These vertices were filtered out but may form valid components in subsequent cycles.

    Only 'vertex_id' is taken from each record - every attribute is re-read from
    vertex_attrs, so the diagnostic fields the worker attached ('reason',
    'best_score') do not reach disk. Order comes from the caller, which sorted
    the list by (component_id, vertex_id).

    Args:
        deleted_list: List of deleted vertex dictionaries
        vertex_attrs: Dictionary with vertex attributes (original coordinates)
        output_path: Path to save deleted vertices JSON
    """
    deleted_vertices = []

    for deleted in deleted_list:
        vertex_id = deleted.get('vertex_id')
        if vertex_id is not None and vertex_id in vertex_attrs:
            deleted_vertices.append(vertex_attrs[vertex_id].copy())

    output = {
        'deleted_vertices': deleted_vertices,
        'total_deleted': len(deleted_vertices)
    }

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n✓ Exported {len(deleted_vertices)} deleted vertices to: {os.path.basename(output_path)}")


def save_graph(graph: ig.Graph, output_path: str, graph_format: str = 'graphml'):
    """
    Dump the raw graph for external inspection (Cytoscape, Gephi, ...).

    DIAGNOSTIC ONLY - written only when the matching --output_* flag is given, and
    read by no downstream step.
    """
    if graph_format == 'graphml':
        graph.write_graphml(output_path)
        print(f"✓ Saved graph (GraphML): {output_path}")
    elif graph_format == 'gml':
        graph.write_gml(output_path)
        print(f"✓ Saved graph (GML): {output_path}")
    elif graph_format == 'edgelist':
        graph.write_edgelist(output_path)
        print(f"✓ Saved graph (EdgeList): {output_path}")
    else:
        raise ValueError(f"Unsupported format: {graph_format}")


# ============================================================================
# SECTION 8 - CLI ENTRY POINT
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Build network graph and filter components (clusterless version). "
                    "Input is edges dictionary from pairing script."
    )
    # --- Required I/O ---
    parser.add_argument("--input_pkl", required=True,
                       help="Input pickle file: edges dictionary {file_idx: [edges]} or paired features file")
    parser.add_argument("--cycle_num", type=int, default=1,
                       help="Current cycle number for tracking")
    parser.add_argument("--output_accepted_components", required=True,
                       help="Output JSON file for accepted components")
    parser.add_argument("--output_deleted_vertices", required=True,
                       help="Output JSON file for deleted vertices (input for next cycle)")

    # --- Optional diagnostic outputs (SECTION 5 / save_graph) ---
    parser.add_argument("--output_graphml", help="Output GraphML file path (optional)")
    parser.add_argument("--output_gml", help="Output GML file path (optional)")
    parser.add_argument("--output_edgelist", help="Output edge list file path (optional)")
    parser.add_argument("--output_analysis",
                       help="DIAGNOSTIC ONLY: graph statistics JSON, for manual inspection. "
                            "Not read by any downstream step (optional)")
    parser.add_argument("--output_composition",
                       help="DIAGNOSTIC ONLY: component composition JSON, for manual inspection. "
                            "Not read by any downstream step (optional)")
    parser.add_argument("--skip-analysis", action='store_true',
                       help="Skip graph analysis computation and statistics (speeds up execution)")

    # --- Filtering (core functionality, SECTION 6) ---
    parser.add_argument("--filter-duplicate-file-vertices", action='store_true',
                       help="Filter duplicate vertices from same file within components")
    parser.add_argument("--mega_complexity_threshold", type=float, default=5000,
                       help="Candidate-combination count at or above which a component is filtered "
                            "by the greedy heuristic instead of exhaustive evaluation. Set in "
                            "production by params.mmf_mega_complexity_threshold in the .nf")
    parser.add_argument("--filter_method", type=str, default='duplicate_only',
                        choices=['duplicate_only', 'anchor_connectivity', 'anchor_nearest',
                                 'bfs', 'window'],
                        help="How a component is reduced. 'duplicate_only' keeps one vertex per "
                             "file and evicts nothing else. The anchor methods pick the vertex "
                             "with the most distinct files inside its tolerance, keep only what "
                             "fits around it, and resolve the remaining duplicates either by "
                             "connectivity (most retained connections, then shortest total "
                             "distance) or by nearest-to-anchor. 'bfs' instead grows outward "
                             "from a start vertex along EDGES, claiming one vertex per file, "
                             "bounded by the budgets below. 'window' instead solves the pairwise "
                             "constraint exactly: it keeps the (m/z, RT) window covering the most "
                             "files, which by construction has every pair inside the cutoffs, then "
                             "picks the tightest one vertex per file inside it.")
    parser.add_argument("--anchor_radius_mode", type=str, default='full',
                        choices=['full', 'half', 'prune'],
                        help="Anchor tolerance. 'full' = the cutoff (no pairwise guarantee - two "
                             "members can sit on opposite sides and be 2x the cutoff apart). "
                             "'half' = cutoff/2, which guarantees every pair conforms. 'prune' = "
                             "full radius, then drop members that clash pairwise.")
    parser.add_argument("--bfs_max_depth", type=int, default=1,
                        help="BFS only. Maximum hops from the start vertex. 1 keeps direct "
                             "neighbours only. Depth is a weak knob: with a half-cutoff budget a "
                             "2-hop path must fit both hops inside it, so going 1 -> 2 moved "
                             "average linkage by 0.002 on a full cycle. 0 or negative = unlimited.")
    parser.add_argument("--bfs_ppm_budget", type=float, default=None,
                        help="BFS only. Maximum ACCUMULATED m/z offset (ppm) along the path from "
                             "the start. By the triangle inequality this bounds each member's "
                             "direct offset, hence every pair at twice it - so half the pairing "
                             "ppm yields a pairwise-conforming component (measured: zero "
                             "violations). Default: half the pairing ppm from the input pkl.")
    parser.add_argument("--bfs_rt_budget", type=float, default=None,
                        help="BFS only. Maximum accumulated RT offset (seconds) along the path "
                             "from the start. Same halving logic as --bfs_ppm_budget. "
                             "Default: half the pairing RT cutoff from the input pkl.")
    parser.add_argument("--window_ppm", type=float, default=None,
                        help="Window only. m/z width of the window, in ppm. Applied to the SET "
                             "(mz_max-mz_min), which is what makes every pair conform - so unlike "
                             "the BFS budget this takes the FULL pairing cutoff rather than half "
                             "of it. Default: the pairing ppm from the input pkl.")
    parser.add_argument("--window_rt", type=float, default=None,
                        help="Window only. RT height of the window, in seconds. Same set-wise "
                             "semantics as --window_ppm. Default: the pairing RT cutoff.")
    parser.add_argument("--window_max_candidates", type=int, default=0,
                        help="Window only. Cap on how many maximum-coverage windows get the inner "
                             "selection, for pathological components (window counts are skewed: "
                             "median ~11 per component, up to ~42k on the largest). 0 = no cap. "
                             "Truncation is deterministic - windows are emitted in sorted order.")
    parser.add_argument("--anchor_ppm", type=float, default=None,
                        help="m/z tolerance for the anchor, in ppm, applied from the anchor's own "
                             "m/z exactly as pairing does it. Defaults to the pairing ppm recorded "
                             "in the input pkl.")
    parser.add_argument("--mz_rt_weight_ratio", type=float, default=1.0,
                        help="RT weight relative to m/z in the scaled space the anchor and "
                             "window methods RANK candidates in (1.0 = equal, < 1.0 "
                             "emphasises m/z). Must match the value pairing ran with - both "
                             "come from params.mmf_mz_rt_weight_ratio in the .nf. Unused by "
                             "--filter_method bfs and duplicate_only, which never rank in "
                             "that space.")
    parser.add_argument("--anchor_rt_cutoff", type=float, default=150.0,
                        help="RT tolerance for the anchor, in seconds. Deliberately wider than the "
                             "pairing RT cutoff: the exported y_center is PRE-alignment, so a "
                             "coordinate-based RT test carries the alignment offset (measured at "
                             "up to ~30s per edge).")
    parser.add_argument("--disable-multiprocessing", action='store_true',
                       help="Disable multiprocessing for filtering (equivalent to running "
                            "everything in one process; overrides --num_workers)")
    parser.add_argument("--num_workers", type=int, default=None,
                       help="Worker processes for duplicate-vertex filtering "
                            f"(default: {DEFAULT_NUM_WORKERS} = run in-process, same as "
                            "--disable-multiprocessing). Clamped to the CPUs actually available "
                            "to the process. Performance only - results are identical for any "
                            "worker count.")

    # --- Metadata inputs ---
    parser.add_argument("--ident-lookup", type=str, default=None,
                       help="JSON file with ident lookup from pairing script")

    args = parser.parse_args()

    print(f"\n{'='*70}")
    print("Building Network Graph (Clusterless)")
    print(f"{'='*70}")

    # Phase timestamps, reported as a breakdown at the end of the run
    timing_marks = {}
    timing_marks['start'] = time.time()

    # ---- Load edges -------------------------------------------------------
    print(f"\nLoading edges and parameters from: {args.input_pkl}")
    timing_marks['load_edges_start'] = time.time()
    edges, file_index_map = load_edges_from_paired_features(args.input_pkl)
    timing_marks['load_edges_end'] = time.time()
    if not edges:
        print("✗ No edges found in input file!")
        return

    print(f"  Loaded {len(edges)} total edges")

    # ---- Build the graph: turn edges into vertices + connections -----------
    # Produces the graph STRUCTURE only. It knows nothing about components yet -
    # they are a property of that structure, extracted two steps below.
    print("\nBuilding graph...")
    timing_marks['build_graph_start'] = time.time()

    # Per-feature metadata (identifications, coordinates, intensity) produced by
    # the pairing step. Optional: without it vertices keep their default attrs.
    ident_lookup = {}
    if args.ident_lookup and os.path.exists(args.ident_lookup):
        try:
            with open(args.ident_lookup, 'r') as f:
                ident_lookup = json.load(f)
            print(f"Loaded ident_lookup from {os.path.basename(args.ident_lookup)}")
            print(f"  Items in ident_lookup: {len(ident_lookup)}")
            for i, (key, value) in enumerate(list(ident_lookup.items())[:3]):
                pep_ident = value.get('pep_ident', [])
                print(f"    {key}: pep_ident={pep_ident}")
        except Exception as e:
            print(f"Warning: Could not load ident_lookup: {e}")
    else:
        print("No ident_lookup file provided or file not found")

    graph, vertex_id_map, vertex_attrs = build_network_graph(edges,
                                                             ident_lookup=ident_lookup,
                                                             file_index_map=file_index_map)
    timing_marks['build_graph_end'] = time.time()

    # Both are multi-GB on large cycles and fully absorbed into the graph and
    # vertex_attrs by now. Releasing them here keeps peak RSS below the filtering
    # stage's fork, which is where memory pressure actually bites.
    print("[MEMORY] Freeing temporary data structures (ident_lookup, edges)...", flush=True)
    del ident_lookup
    del edges
    gc.collect()
    print("[MEMORY] ✓ Freed temporary structures", flush=True)

    if graph is None:
        return

    # ---- Diagnostic analysis (SECTION 5) ----------------------------------
    analysis = None
    graph_composition = None
    # Neither graph_analysis.json nor graph_composition.json is read by any
    # downstream step - they exist for manual inspection. --skip-analysis turns
    # the whole block off (analyze_graph's diameter() is the expensive part).
    if not args.skip_analysis:
        print("\nAnalyzing graph...")
        timing_marks['analyze_start'] = time.time()
        analysis = analyze_graph(graph, vertex_attrs)
        print_analysis(analysis)

        print("\nAnalyzing graph composition...")
        graph_composition = analyze_graph_composition(graph, vertex_attrs)
        timing_marks['analyze_end'] = time.time()

        comp_summary = graph_composition.get('composition_summary', {})
        if comp_summary:
            print(f"  Total connected components: {comp_summary.get('total_components')}")
            print(f"  Total features in components: {comp_summary.get('total_features')}")
            print(f"  Average features per component: {comp_summary.get('avg_features_per_component', 0):.2f}")
            print(f"  Average files per component: {comp_summary.get('avg_files_per_component', 0):.2f}")
            print(f"  Components by size: {comp_summary.get('components_by_size', {})}")
    else:
        print("Skipping graph analysis (--skip-analysis enabled)")

    # ---- Find the connected components -------------------------------------
    # This is the clusterless step: each connected component is reduced into one feature group via filtering of vertices from the same file origin. 
    # Compleate Filtered Groups are saved as such and removed from the active pool.
    # Vertices removed by the filter are collected and used as input for the next cycle (orchestrator script), where they are paired (pairing script) again.
    # No community detection runs on top to split them further.
    #
    # components[i] is a plain list of vertex ids, so it is handed straight to
    # filtering - nothing is copied into a parallel structure.
    #
    # The graph stays in use after this: graph.degree() is the score behind both
    # selection methods, which is why it is shared into the worker processes
    # (SECTION 2). What is settled here is MEMBERSHIP. components() is not called
    # again and the graph is never mutated, so deleting a vertex cannot split a
    # component - only which vertices survive is still open.
    print("\n" + "="*70)
    print(f"Component Building (Clusterless) - Cycle {args.cycle_num}")
    print("="*70)
    timing_marks['components_start'] = time.time()

    components = graph.components()
    print(f"  Found {len(components)} connected components")

    timing_marks['components_end'] = time.time()

    # ---- Duplicate-vertex filtering (SECTION 6) ---------------------------
    if args.filter_duplicate_file_vertices:
        print("\n" + "="*70)
        print("Filtering Duplicate File Vertices")
        print("="*70)
        timing_marks['filter_start'] = time.time()

        # ---- Assemble the anchor configuration -----------------------------
        # The ranking metric needs one axis scale per dimension plus the weight
        # ratio. The scales come from the pkl, where pairing derived them from EVERY
        # input feature; the weight ratio comes straight from the workflow parameter
        # via --mz_rt_weight_ratio. A pkl without the scales falls back to
        # recomputing from this graph's own vertices - close, but not identical,
        # since the graph only ever sees the PAIRED features.
        anchor_cfg = None
        if args.filter_method != 'duplicate_only':
            sc = pairing_scaling_constants or {}
            range_x, range_y = sc.get('range_x'), sc.get('range_y')
            if not range_x or not range_y:
                print("  ⚠ Input pkl carries no scaling constants; recomputing from "
                      "paired vertices only.")
                xs = [a['x_center'] for a in vertex_attrs.values()]
                ys = [a['y_center'] for a in vertex_attrs.values()]
                range_x = (max(xs) - min(xs)) or 1.0
                range_y = (max(ys) - min(ys)) or 1.0

            # Coordinates reach this script only through the ident lookup. Without
            # it every x_center/y_center defaults to 0.0, which would leave the
            # anchor comparing every vertex against every other as identical -
            # silently meaningless rather than obviously broken. Refuse instead.
            # Checked against the VERTICES, not the scales: a pkl carrying good
            # scales says nothing about whether this run got its coordinates.
            if len({round(a['x_center'], 6) for a in vertex_attrs.values()}) <= 1:
                print("  ✗ ERROR: --filter_method needs feature coordinates, and every "
                      "vertex has the same x_center. Pass --ident-lookup so coordinates "
                      "are loaded, or use --filter_method duplicate_only.")
                sys.exit(1)

            anchor_ppm = args.anchor_ppm
            if anchor_ppm is None:
                anchor_ppm = (pairing_cutoff_params or {}).get('ppm')
            if anchor_ppm is None:
                print("  ✗ ERROR: --filter_method needs an m/z tolerance, and the input pkl "
                      "records no pairing ppm. Pass --anchor_ppm explicitly.")
                sys.exit(1)

            # BFS budgets default to HALF the pairing cutoffs, which is the setting
            # that makes finished components pairwise-conforming rather than merely
            # bounded at twice the cutoff.
            pairing_ppm = (pairing_cutoff_params or {}).get('ppm')
            pairing_rt = (pairing_cutoff_params or {}).get('rt_cutoff')
            bfs_ppm = args.bfs_ppm_budget
            if bfs_ppm is None and pairing_ppm:
                bfs_ppm = float(pairing_ppm) / 2.0
            bfs_rt = args.bfs_rt_budget
            if bfs_rt is None and pairing_rt:
                bfs_rt = float(pairing_rt) / 2.0
            if args.filter_method == 'bfs' and (bfs_ppm is None or bfs_rt is None):
                print("  ✗ ERROR: --filter_method bfs needs m/z and RT budgets, and the input "
                      "pkl records no pairing cutoffs to halve. Pass --bfs_ppm_budget and "
                      "--bfs_rt_budget explicitly.")
                sys.exit(1)

            win_ppm = args.window_ppm if args.window_ppm is not None else pairing_ppm
            win_rt = args.window_rt if args.window_rt is not None else pairing_rt
            if args.filter_method == 'window' and (win_ppm is None or win_rt is None):
                print("  ✗ ERROR: --filter_method window needs m/z and RT widths, and the input "
                      "pkl records no pairing cutoffs to fall back on. Pass --window_ppm and "
                      "--window_rt explicitly.")
                sys.exit(1)

            anchor_cfg = {
                'method': args.filter_method,
                'window_ppm': float(win_ppm) if win_ppm is not None else None,
                'window_rt': float(win_rt) if win_rt is not None else None,
                'window_max_candidates': args.window_max_candidates,
                'bfs_max_depth': (args.bfs_max_depth if args.bfs_max_depth
                                  and args.bfs_max_depth > 0 else None),
                'bfs_ppm_budget': bfs_ppm,
                'bfs_rt_budget': bfs_rt,
                'radius_mode': args.anchor_radius_mode,
                'anchor_ppm': float(anchor_ppm),
                'anchor_rt_cutoff': float(args.anchor_rt_cutoff),
                'range_x': range_x, 'range_y': range_y,
                'mz_rt_weight_ratio': args.mz_rt_weight_ratio,
            }

        accepted_components, deleted_vertices_list = filter_duplicate_file_vertices(
            graph,
            components,
            vertex_attrs,
            args.num_workers,
            not args.disable_multiprocessing,
            args.mega_complexity_threshold,
            anchor_cfg
        )

        print(f"  ✓ Filtered {len(deleted_vertices_list)} duplicate vertices")
        # Every vertex belongs to exactly one component and can be deleted at
        # most once, so the survivors are simply the remainder.
        print(f"  ✓ Retained {graph.vcount() - len(deleted_vertices_list)} vertices")

        timing_marks['filter_end'] = time.time()
    else:
        print("\nFiltering disabled (use --filter-duplicate-file-vertices to enable)")
        timing_marks['filter_start'] = None
        timing_marks['filter_end'] = None

        # Nothing filtered: every component is accepted exactly as found
        accepted_components = {i: {'vertices': v} for i, v in enumerate(components)}
        deleted_vertices_list = []

    # ---- Export (SECTION 7) -----------------------------------------------
    print("\n" + "="*70)
    print("Exporting Results")
    print("="*70)
    timing_marks['export_start'] = time.time()

    export_accepted_components(
        accepted_components,
        vertex_attrs,
        args.output_accepted_components,
        args.cycle_num
    )

    export_deleted_vertices_for_next_cycle(
        deleted_vertices_list,
        vertex_attrs,
        args.output_deleted_vertices
    )

    # Optional debugging artefacts
    if args.output_graphml:
        save_graph(graph, args.output_graphml, graph_format='graphml')
        print(f"✓ Saved GraphML: {args.output_graphml}")

    if args.output_gml:
        save_graph(graph, args.output_gml, graph_format='gml')
        print(f"✓ Saved GML: {args.output_gml}")

    if args.output_edgelist:
        save_graph(graph, args.output_edgelist, graph_format='edgelist')
        print(f"✓ Saved edge list: {args.output_edgelist}")

    if args.output_analysis and analysis is not None:
        with open(args.output_analysis, 'w') as f:
            json.dump(analysis, f, indent=2)
        print(f"✓ Saved analysis: {args.output_analysis}")

    if args.output_composition and graph_composition is not None:
        with open(args.output_composition, 'w') as f:
            json.dump(graph_composition, f, indent=2)
        print(f"✓ Saved graph composition: {args.output_composition}")

    timing_marks['export_end'] = time.time()

    # ---- Timing breakdown --------------------------------------------------
    print(f"\n{'='*70}")
    print("TIMING SUMMARY")
    print(f"{'='*70}")

    timing_phases = [
        ('Loading edges', 'load_edges_start', 'load_edges_end'),
        ('Building graph', 'build_graph_start', 'build_graph_end'),
        ('Graph analysis', 'analyze_start', 'analyze_end'),
        ('Component building', 'components_start', 'components_end'),
        ('Filtering duplicates', 'filter_start', 'filter_end'),
        ('Export/Save', 'export_start', 'export_end'),
    ]

    total_elapsed = time.time() - timing_marks['start']
    phase_times = {}

    for phase_name, start_key, end_key in timing_phases:
        if start_key in timing_marks and end_key in timing_marks and timing_marks[start_key] and timing_marks[end_key]:
            elapsed = timing_marks[end_key] - timing_marks[start_key]
            phase_times[phase_name] = elapsed
            pct = (elapsed / total_elapsed * 100) if total_elapsed > 0 else 0
            print(f"  {phase_name:.<40} {elapsed:>8.2f}s ({pct:>5.1f}%)")
        else:
            # The two phases that can legitimately be turned off
            if phase_name in ['Graph analysis', 'Filtering duplicates']:
                print(f"  {phase_name:.<40} {'SKIPPED':>8s}")

    print(f"{'-'*70}")
    print(f"  {'TOTAL':.<40} {total_elapsed:>8.2f}s")
    print(f"{'='*70}\n")

    if phase_times:
        longest_phase = max(phase_times, key=phase_times.get)
        longest_time = phase_times[longest_phase]
        pct = (longest_time / total_elapsed * 100) if total_elapsed > 0 else 0
        print(f"⏱️  BOTTLENECK: '{longest_phase}' takes {longest_time:.2f}s ({pct:.1f}% of total)\n")

    print(f"✅ Cycle {args.cycle_num} complete!")
    print(f"   Accepted components: {len(accepted_components)}")
    print(f"   Deleted vertices: {len(deleted_vertices_list)}")

    return graph, vertex_attrs, analysis


if __name__ == "__main__":
    main()

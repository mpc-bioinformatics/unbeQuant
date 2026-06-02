#!/usr/bin/env python3
"""
Build a network graph from paired feature edges using igraph.
Analyze and cluster feature connections across files.

Refactored version with visualization code removed for reduced memory footprint.
"""

import argparse
import pickle
import json
import numpy as np
import os
import shutil
import sys
import math
import gc
from pathlib import Path
from typing import Dict, List, Tuple
from pdb import set_trace as bp
import time
from multiprocessing import Pool, get_context, Value, Lock
import random

try:
    import igraph as ig
except ImportError:
    print("igraph not found. Install with: pip install python-igraph")
    exit(1)

try:
    from scipy.optimize import minimize_scalar
except ImportError:
    minimize_scalar = None

# ============================================================================
# MODULE-LEVEL GLOBALS (shared via fork/copy-on-write with worker processes)
# These are populated in filter_duplicate_file_vertices() before creating Pool
# Workers inherit them automatically - NO need to pass via initargs (saves 2GB!)
# ============================================================================
edge_index = None
vertex_attrs = None
vertex_to_cluster_map = None
filter_method = None
component_cluster_info = None


def load_edges_and_cutoff_params_from_paired_features(paired_features_file: str) -> Tuple[List[Dict], Dict]:
    """
    Load BOTH edges and cutoff parameters in a SINGLE file read operation.
    This avoids loading the pickle file twice and reduces peak memory usage.
    
    Args:
        paired_features_file: Path to pickle or JSON file
    
    Returns:
        Tuple of (edges_list, cutoff_params_dict)
        - edges_list: List of edge dictionaries
        - cutoff_params_dict: Dict with keys mz_cutoff, rt_cutoff, edges_cutoff (or empty if not found)
    """
    # Load file once
    if paired_features_file.endswith('.pkl'):
        with open(paired_features_file, 'rb') as f:
            data = pickle.load(f)
    elif paired_features_file.endswith('.json'):
        with open(paired_features_file, 'r') as f:
            data = json.load(f)
    else:
        raise ValueError("File must be .pkl or .json")
    
    # Extract edges
    edges = []
    if isinstance(data, dict):
        # Check if it's the new wrapped edges format {network_edges: {file_idx: [edges]}, cutoff_params: {...}}
        if 'network_edges' in data and 'cutoff_params' in data:
            # Wrapped format with cutoff parameters from --skip-matchfinder output
            network_edges = data['network_edges']
            if all(isinstance(v, list) for v in network_edges.values()):
                # Dictionary format: flatten all edges (keys can be int or str)
                for file_idx, file_edges in network_edges.items():
                    edges.extend(file_edges)
            else:
                edges = network_edges
        # Check if it's the edges dictionary format {file_idx: [edges]} (int or str keys)
        elif all(isinstance(v, list) for v in data.values()):
            # Dictionary format: flatten all edges (keys can be int or str)
            for file_idx, file_edges in data.items():
                edges.extend(file_edges)
        else:
            # Paired features format with metadata: data['network_edges']
            edges = data.get('network_edges', [])
    elif isinstance(data, list):
        # If it's a list itself
        edges = data
    
    # Extract cutoff parameters
    cutoff_params = {}
    if isinstance(data, dict) and 'cutoff_params' in data:
        cutoff_params = data['cutoff_params']
    
    return edges, cutoff_params


def build_network_graph(edges: List[Dict], edge_cutoff: float = float('inf'), 
                       mz_cutoff: float = None, rt_cutoff: float = None,
                       ident_lookup: Dict = None) -> Tuple[ig.Graph, Dict]:
    """
    Build an igraph network graph from edges with flexible cutoff filtering.
    
    Args:
        edges: List of edge dictionaries with format:
               {'file1': {filename, feature_idx}, 'file2': {...}, 'distance': float, 'file1_idx': int, 'file2_idx': int}
               or old format:
               {'current_file': {filename, feature_idx}, 'matched_file': {...}, 'distance': float, ...}
        edge_cutoff: Maximum euclidean distance to include edges. Default: no cutoff.
        mz_cutoff: Maximum m/z (x-axis) coordinate difference. If specified with rt_cutoff, overrides edge_cutoff.
        rt_cutoff: Maximum RT (y-axis) coordinate difference. If specified with mz_cutoff, overrides edge_cutoff.
        ident_lookup: Combined ident lookup for pep_ident/prot_ident metadata
    
    Returns:
        Tuple of (igraph.Graph, vertex_id_map, vertex_attrs)
    """
    
    # Filter edges based on cutoff mode
    if mz_cutoff is not None and rt_cutoff is not None:
        # Coordinate-based filtering: filter by m/z and RT distances instead of euclidean distance
        print(f"Using coordinate-based filtering: mz_cutoff={mz_cutoff}, rt_cutoff={rt_cutoff}")
        filtered_edges = [e for e in edges 
                         if e.get('mz_distance', 0) <= mz_cutoff and e.get('rt_distance', 0) <= rt_cutoff]
        print(f"Loaded {len(edges)} total edges, {len(filtered_edges)} pass coordinate cutoff (mz≤{mz_cutoff}, rt≤{rt_cutoff})")
    else:
        # Distance-based filtering (default)
        filtered_edges = [e for e in edges if e.get('distance', 0) <= edge_cutoff]
        print(f"Loaded {len(edges)} total edges, {len(filtered_edges)} pass distance cutoff {edge_cutoff}")
    
    if not filtered_edges:
        print("No edges pass the cutoff threshold!")
        return None, {}, {}
    
    # Build vertex set - each vertex represents a unique (filename, feature_idx) pair
    vertices = set()
    vertex_id_map = {}  # (filename, feature_idx) -> vertex_id
    vertex_attrs = {}   # vertex_id -> {'filename': str, 'feature_idx': int, 'label': str}
    
    for edge in filtered_edges:
        # Handle both old and new edge formats
        if 'file1' in edge:
            file1_info = edge['file1']
            file2_info = edge['file2']
        else:
            file1_info = edge['current_file']
            file2_info = edge['matched_file']
        
        # Create unique identifiers for vertices
        vertex1_key = (file1_info['filename'], file1_info['feature_idx'])
        vertex2_key = (file2_info['filename'], file2_info['feature_idx'])
        
        vertices.add(vertex1_key)
        vertices.add(vertex2_key)
    
    # Assign vertex IDs
    for idx, vertex_key in enumerate(sorted(vertices)):
        vertex_id_map[vertex_key] = idx
        filename, feature_idx = vertex_key
        vertex_attrs[idx] = {
            'filename': filename,
            'feature_idx': feature_idx,
            'label': f"{Path(filename).stem}_{feature_idx}",
            'pep_ident': [],  # Default empty, will be populated from lookup if available
            'prot_ident': [],  # Default empty, will be populated from lookup if available
            'x_center': 0.0,  # RT value, will be populated from lookup if available
            'y_center': 0.0,  # mz value, will be populated from lookup if available
            'intensity': 0.0,
            'openms_fid': '',
            'ms2_scans': '',
            'charge': 0
        }
        
        # Load identifications and coordinates from combined lookup if available
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
    
    # Debug: Show pep_ident loading status
    pep_ident_loaded_count = 0
    if ident_lookup:
        for idx, attrs in vertex_attrs.items():
            pep = attrs.get('pep_ident', [])
            if pep and len(pep) > 0:
                pep_ident_loaded_count += 1
        
        print(f"  ✓ Loaded pep_ident for: {pep_ident_loaded_count}/{len(vertex_attrs)} vertices")
        sys.stdout.flush()
        
        # Show first 3 vertices with pep_ident
        shown = 0
        for idx, attrs in vertex_attrs.items():
            pep = attrs.get('pep_ident', [])
            if pep and len(pep) > 0:
                print(f"    Vertex {idx}: {attrs['label']} -> pep_ident={pep}")
                shown += 1
                if shown >= 3:
                    break
        sys.stdout.flush()
    
    # Create graph with vertices
    graph = ig.Graph(len(vertices))
    
    # Add vertex attributes
    for vertex_id, attrs in vertex_attrs.items():
        for attr_name, attr_value in attrs.items():
            graph.vs[vertex_id][attr_name] = attr_value
    
    # Add edges with distance weights
    edge_list = []
    edge_weights = []
    
    for edge in filtered_edges:
        # Handle both old and new edge formats
        if 'file1' in edge:
            file1_info = edge['file1']
            file2_info = edge['file2']
        else:
            # Old format: explicitly directed current_file → matched_file
            file1_info = edge['current_file']
            file2_info = edge['matched_file']
        
        vertex1_key = (file1_info['filename'], file1_info['feature_idx'])
        vertex2_key = (file2_info['filename'], file2_info['feature_idx'])
        
        vertex1_id = vertex_id_map[vertex1_key]
        vertex2_id = vertex_id_map[vertex2_key]
        
        # Edge points from file1/current_file (source) to file2/matched_file (target)
        edge_list.append((vertex1_id, vertex2_id))
        edge_weights.append(edge['distance'])
    
    # Add all edges to graph at once
    graph.add_edges(edge_list)
    
    # Add edge weights
    for edge_id, weight in enumerate(edge_weights):
        graph.es[edge_id]['distance'] = weight
    
    print(f"Created graph with {graph.vcount()} vertices and {graph.ecount()} edges")
    print(f"  Edge direction: current_file (source) → matched_file (target)")
    
    return graph, vertex_id_map, vertex_attrs


def analyze_graph(graph: ig.Graph, vertex_attrs: Dict) -> Dict:
    """Analyze network graph properties."""
    
    if graph is None or graph.ecount() == 0:
        return {}
    
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
    
    # Get connected components
    components = graph.components()
    analysis['components'] = {
        'num_components': len(components),
        'component_sizes': [len(c) for c in components],
        'largest_component_size': max(len(c) for c in components) if components else 0
    }
    
    # Get degree distribution per file
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
    
    # Clustering coefficient
    try:
        cc = graph.transitivity_undirected()
        analysis['clustering_coefficient'] = cc
    except:
        analysis['clustering_coefficient'] = None
    
    # Shortest path statistics
    try:
        if graph.is_connected():
            analysis['diameter'] = graph.diameter()
        else:
            # For disconnected graphs, compute diameter for largest component
            largest_comp = graph.induced_subgraph(components[0])
            analysis['diameter'] = largest_comp.diameter() if len(components[0]) > 1 else 0
    except:
        analysis['diameter'] = None
    
    return analysis


def analyze_graph_composition(graph: ig.Graph, vertex_attrs: Dict) -> Dict:
    """
    Analyze graph composition with detailed file and feature multiplicity information.
    
    Returns enriched statistics including file uniqueness, multiplicity patterns, and composition breakdown
    needed for interactive visualization:
    {
        "component_composition": {
            "3_3": 5,  # 5 components with 3 features from 3 files (aggregated)
            ...
        },
        "composition_summary": {... overall stats ...},
        "file_match_uniqueness": {
            "2": {  # component_size (2 features per component)
                "all_unique": 4,  # components where each file has exactly 1 feature
                "with_duplicates": 1,  # components where any file has 2+ features
                "total": 5,
                "file_composition": {
                    "file1.tsv, file2.tsv": 4,  # composition pattern -> count
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
    composition_map = {}  # Key: f"{num_features}_{num_files}", Value: count
    file_match_uniqueness_by_size = {}  # component_size -> detailed analysis
    features_per_component = []
    files_per_component = []
    components_by_size = {}  # Key: num_features, Value: count
    
    for component_vertices in components:
        # Count files and analyze file multiplicity in this component
        file_counts = {}  # filename -> count of features from that file
        feature_ids = []  # Track all feature IDs in component
        
        for vertex_id in component_vertices:
            filename = vertex_attrs[vertex_id]['filename']
            feature_id = vertex_attrs[vertex_id].get('feature_id', f'v{vertex_id}')
            file_counts[filename] = file_counts.get(filename, 0) + 1
            feature_ids.append(feature_id)
        
        num_features = len(component_vertices)
        num_files = len(file_counts)
        
        # Create composition key for aggregate counting
        comp_key = f"{num_features}_{num_files}"
        composition_map[comp_key] = composition_map.get(comp_key, 0) + 1
        
        features_per_component.append(num_features)
        files_per_component.append(num_files)
        
        # Track components by size
        size_key = str(num_features)
        components_by_size[size_key] = components_by_size.get(size_key, 0) + 1
        
        # Initialize uniqueness tracking for this component size if needed
        if size_key not in file_match_uniqueness_by_size:
            file_match_uniqueness_by_size[size_key] = {
                'all_unique': 0,
                'with_duplicates': 0,
                'total': 0,
                'file_composition': {},
                'multiplicity_distribution': {},
                'file_appearance': {}
            }
        
        # Track this component
        file_match_uniqueness_by_size[size_key]['total'] += 1
        
        # Check if all files are unique (each appears exactly once)
        multiplicities = sorted(file_counts.values(), reverse=True)
        is_all_unique = all(count == 1 for count in multiplicities)
        
        if is_all_unique:
            file_match_uniqueness_by_size[size_key]['all_unique'] += 1
        else:
            file_match_uniqueness_by_size[size_key]['with_duplicates'] += 1
        
        # Track file composition (which files are in this component)
        files_in_group = sorted(set(file_counts.keys()))
        file_comp_str = ', '.join(files_in_group)
        if file_comp_str not in file_match_uniqueness_by_size[size_key]['file_composition']:
            file_match_uniqueness_by_size[size_key]['file_composition'][file_comp_str] = 0
        file_match_uniqueness_by_size[size_key]['file_composition'][file_comp_str] += 1
        
        # Track multiplicity pattern (e.g., "1x3" for 3 files each with 1 feature, or "2x1,1x1")
        # Format: count of unique multiplicities, repeated as needed
        multiplicity_pattern = ','.join([f"{mult}x{multiplicities.count(mult)}" 
                                        for mult in sorted(set(multiplicities), reverse=True)])
        if multiplicity_pattern not in file_match_uniqueness_by_size[size_key]['multiplicity_distribution']:
            file_match_uniqueness_by_size[size_key]['multiplicity_distribution'][multiplicity_pattern] = 0
        file_match_uniqueness_by_size[size_key]['multiplicity_distribution'][multiplicity_pattern] += 1
        
        # Track file appearance (which files appear in this component size)
        for fname in file_counts.keys():
            if fname not in file_match_uniqueness_by_size[size_key]['file_appearance']:
                file_match_uniqueness_by_size[size_key]['file_appearance'][fname] = 0
            file_match_uniqueness_by_size[size_key]['file_appearance'][fname] += 1
    
    # Build summary statistics
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



def _compute_cluster_edge_weights(subgraph: ig.Graph, weight_mode: str) -> List[float]:
    """Compute edge weights for clustering from the 'distance' attribute."""
    weights = []
    for edge in subgraph.es:
        distance = edge['distance'] if 'distance' in edge.attributes() else 1.0
        if weight_mode == 'distance':
            weight = float(distance)
        else:
            # Default: inverse distance (stronger weight for closer nodes)
            weight = 1.0 / (1.0 + float(distance))
        weights.append(weight)
    return weights


def detect_conflicting_pep_idents_in_cluster(vertex_list: List[int], vertex_attrs: Dict) -> bool:
    """
    Detect if a cluster has conflicting pep_idents.
    
    A conflict exists if vertices with pep_idents don't share at least one common pep_ident.
    Uses same logic as generate_pairing_report.py.
    """
    pep_ident_sets = []
    
    for vertex_id in vertex_list:
        if vertex_id in vertex_attrs:
            pep_ident = vertex_attrs[vertex_id].get('pep_ident')
            if pep_ident:
                if isinstance(pep_ident, list):
                    pep_set = set(pep_ident) if pep_ident else set()
                else:
                    pep_set = {pep_ident}
                
                if pep_set:
                    pep_ident_sets.append(pep_set)
    
    # If 0 or 1 vertices have pep_idents, no conflict
    if len(pep_ident_sets) <= 1:
        return False
    
    # Find intersection of all pep_ident sets
    common_pep_idents = pep_ident_sets[0]
    for pep_set in pep_ident_sets[1:]:
        common_pep_idents = common_pep_idents.intersection(pep_set)
        if not common_pep_idents:
            break
    
    return not common_pep_idents


def detect_cluster_pep_ident_duplicates_in_component(comp_info: Dict, vertex_attrs: Dict) -> bool:
    """
    Detect if clusters within a component share identical pep_idents.
    Returns True if identical pep_idents are found in different clusters.
    """
    cluster_pep_idents = {}  # Maps cluster_id -> set of pep_idents
    
    for cluster_id, vertex_list in comp_info.get('clusters', {}).items():
        pep_idents = set()
        for vertex_id in vertex_list:
            if vertex_id in vertex_attrs:
                pep_ident = vertex_attrs[vertex_id].get('pep_ident')
                if pep_ident:
                    if isinstance(pep_ident, list):
                        pep_idents.update(pep_ident)
                    else:
                        pep_idents.add(pep_ident)
        
        if pep_idents:
            cluster_pep_idents[cluster_id] = pep_idents
    
    # If fewer than 2 clusters have pep_idents, no duplicates possible
    if len(cluster_pep_idents) <= 1:
        return False
    
    # Check if any pep_idents are shared between different clusters
    cluster_ids = list(cluster_pep_idents.keys())
    for i in range(len(cluster_ids)):
        for j in range(i + 1, len(cluster_ids)):
            shared = cluster_pep_idents[cluster_ids[i]].intersection(cluster_pep_idents[cluster_ids[j]])
            if shared:
                return True  # Found shared pep_idents
    
    return False


def create_cluster_size_distribution(cluster_sizes: List[int]) -> Dict[int, int]:
    """
    Create a distribution dictionary from a list of cluster sizes.
    
    Args:
        cluster_sizes: List of cluster sizes
    
    Returns:
        Dictionary mapping cluster size -> count of clusters of that size
    """
    distribution = {}
    for size in cluster_sizes:
        distribution[size] = distribution.get(size, 0) + 1
    return distribution


def recover_deleted_vertices_edge_based(graph: ig.Graph, deleted_vertices_list: List[Dict],
                                        component_cluster_info: Dict, vertex_attrs: Dict,
                                        vertex_to_cluster_map: Dict) -> Tuple[Dict, int]:
    """
    Recover deleted vertices using edge-based logic (proper recovery).
    
    For each deleted vertex:
    1. Find all edges to vertices in other clusters
    2. Sum edge weights by cluster (cluster with strongest connection wins)
    3. Add to that cluster IF it doesn't already have a vertex from the same file
    
    Args:
        graph: igraph Graph object
        deleted_vertices_list: List of deleted vertex dictionaries
        component_cluster_info: Component clustering information
        vertex_attrs: Vertex attributes dictionary
        vertex_to_cluster_map: Vertex to (comp_idx, cluster_id) mapping
    
    Returns:
        Tuple of (updated_component_cluster_info, total_recovered_count)
    """
    if not deleted_vertices_list:
        return component_cluster_info, 0
    
    recovered_cluster_info = {}
    
    # Deep copy the cluster structure
    for comp_idx, info in component_cluster_info.items():
        recovered_cluster_info[comp_idx] = {
            'num_clusters': info['num_clusters'],
            'clusters': {cid: list(vlist) for cid, vlist in info['clusters'].items()},
            'method': info.get('method'),
            'modularity': info.get('modularity'),
            'vertices': list(info.get('vertices', [])),
            'details': info.get('details', {})
        }
    
    total_recovered = 0
    
    for deleted_vertex_info in deleted_vertices_list:
        deleted_vertex_id = deleted_vertex_info['vertex_id']
        deleted_filename = deleted_vertex_info['filename']
        orig_comp_idx = deleted_vertex_info['component_id']
        orig_cluster_id = deleted_vertex_info['cluster_id']
        
        # Find edges from this deleted vertex to other clusters
        neighbors = graph.neighbors(deleted_vertex_id)
        
        # Build map: (comp_idx, cluster_id) -> weighted_edge_sum
        cluster_weighted_edges = {}
        
        for neighbor_id in neighbors:
            # Check if neighbor is still in graph
            if neighbor_id not in vertex_to_cluster_map:
                continue
            
            neighbor_comp_idx, neighbor_cluster_id = vertex_to_cluster_map[neighbor_id]
            
            # Only for different clusters within same component
            if neighbor_comp_idx != orig_comp_idx or neighbor_cluster_id == orig_cluster_id:
                continue
            
            # Get edge weight
            edge_id = graph.get_eid(deleted_vertex_id, neighbor_id, directed=False, error=False)
            if edge_id == -1:
                continue
            
            edge_weight = graph.es[edge_id]['weight'] if 'weight' in graph.es[edge_id].attributes() else 1.0
            
            # Sum weighted edges to this cluster
            cluster_key = (neighbor_comp_idx, neighbor_cluster_id)
            cluster_weighted_edges[cluster_key] = cluster_weighted_edges.get(cluster_key, 0.0) + edge_weight
        
        if not cluster_weighted_edges:
            # No edges to other clusters - skip recovery
            continue
        
        # Find cluster with strongest weighted connection
        best_cluster_key = max(cluster_weighted_edges, key=cluster_weighted_edges.get)
        target_comp_idx, target_cluster_id = best_cluster_key
        
        # Check if target cluster already has vertex from same file
        target_vertices = recovered_cluster_info[target_comp_idx]['clusters'][target_cluster_id]
        file_already_present = any(
            vertex_attrs[v]['filename'] == deleted_filename
            for v in target_vertices
        )
        
        if file_already_present:
            # Skip recovery - file conflict
            continue
        
        # Recover the vertex!
        recovered_cluster_info[target_comp_idx]['clusters'][target_cluster_id].append(deleted_vertex_id)
        recovered_cluster_info[target_comp_idx]['vertices'].append(deleted_vertex_id)
        
        # Update vertex map
        vertex_to_cluster_map[deleted_vertex_id] = (target_comp_idx, target_cluster_id)
        vertex_attrs[deleted_vertex_id]['recovered'] = True
        vertex_attrs[deleted_vertex_id]['original_cluster'] = orig_cluster_id
        
        total_recovered += 1
    
    return recovered_cluster_info, total_recovered


def evaluate_resolution(graph: ig.Graph, vertex_attrs: Dict,
                       resolution: float,
                       method: str = 'leiden',
                       use_weights: bool = False,
                       weight_mode: str = 'inverse',
                       objective_function: str = 'CPM',
                       apply_filtering: bool = False,
                       filter_method: str = 'simplified-modularity',
                       total_features_after_filter: int = None,
                       unpaired_vertices: List[Dict] = None,
                       num_cores: int = None,
                       use_multiprocessing: bool = True,
                       initial_seed: int = 42,
                       mega_complexity_threshold: float = 10000,
                       iteration_output_dir: str = None,
                       iteration_num: int = None,
                       pre_computed_analysis: Dict = None,
                       pre_computed_composition: Dict = None) -> Dict:
    """
    Evaluate clustering quality at a given resolution parameter.
    
    Args:
        graph: igraph Graph object
        vertex_attrs: Dictionary of vertex attributes
        resolution: Resolution parameter to test
        method: Clustering method (default: leiden)
        use_weights: Whether to use edge weights
        weight_mode: Weight mode (inverse or distance)
        objective_function: Objective function for leiden (CPM or modularity)
        apply_filtering: Whether to apply duplicate vertex filtering
        filter_method: Method for filtering (simplified-modularity or modularity)
        total_features_after_filter: Total features after filtering for percentage calculation
        unpaired_vertices: List of unpaired vertex dictionaries (default: None/empty)
    
    Returns:
        Dictionary with metrics: {
            'resolution': float,
            'vertices_deleted': int,
            'avg_cluster_size': float,
            'score': float,
            'iteration_dir': str
        }
    """
    # Reset random seed for this resolution to ensure deterministic clustering
    # Use the same seed unmodified in each iteration (as provided via --random-seed)
    import random
    random.seed(initial_seed)
    
    # Perform clustering at this resolution
    vertex_to_cluster_map, component_cluster_info = cluster_graph_independent_components(
        graph,
        vertex_attrs,
        method=method,
        use_weights=use_weights,
        weight_mode=weight_mode,
        resolution_parameter=resolution,
        objective_function=objective_function,
        random_seed=initial_seed
    )
    
    # Calculate initial (before filtering) cluster sizes for distribution
    initial_cluster_sizes = []
    for comp_id, comp_info in component_cluster_info.items():
        for cluster_id, vertex_list in comp_info.get('clusters', {}).items():
            initial_cluster_sizes.append(len(vertex_list))
    initial_cluster_distribution = create_cluster_size_distribution(initial_cluster_sizes)
    
    # Apply filtering if enabled
    deleted_vertices = 0
    deleted_list = []
    filtered_cluster_info = None
    filtered_vertex_map = None
    filtered_cluster_sizes = initial_cluster_sizes.copy()  # Default: same as initial
    
    if apply_filtering and vertex_to_cluster_map:
        filtered_vertex_map, filtered_cluster_info, deleted_list = filter_duplicate_file_vertices(
            graph, vertex_to_cluster_map, component_cluster_info, vertex_attrs, filter_method, num_cores, use_multiprocessing,
            mega_complexity_threshold=mega_complexity_threshold
        )
        deleted_vertices = len(deleted_list)
        component_cluster_info = filtered_cluster_info  # Use filtered results
        
        # Calculate filtered (after filtering, before recovery) cluster sizes
        filtered_cluster_sizes = []
        for comp_id, comp_info in component_cluster_info.items():
            for cluster_id, vertex_list in comp_info.get('clusters', {}).items():
                filtered_cluster_sizes.append(len(vertex_list))
    
    filtered_cluster_distribution = create_cluster_size_distribution(filtered_cluster_sizes)
    
    # RECOVERY PHASE: Use real edge-based recovery
    # Deleted vertices are added back to clusters they're most strongly connected to
    total_recovered = 0
    recovery_map = None
    if deleted_list:
        # Use the filtered map if available, otherwise build from current component_cluster_info
        recovery_map = filtered_vertex_map
        if not recovery_map:
            recovery_map = {}
            for comp_id, comp_info in component_cluster_info.items():
                for cluster_id, vertex_ids in comp_info['clusters'].items():
                    for v_id in vertex_ids:
                        recovery_map[v_id] = (comp_id, cluster_id)
        
        # Perform edge-based recovery
        component_cluster_info, total_recovered = recover_deleted_vertices_edge_based(
            graph, deleted_list, component_cluster_info, vertex_attrs, recovery_map
        )
    
    # Calculate overall average cluster size (all clusters, AFTER recovery)
    # This includes both filtered clusters and recovered vertices added to existing clusters
    all_cluster_sizes = []
    for comp_id, comp_info in component_cluster_info.items():
        for cluster_id, vertex_list in comp_info['clusters'].items():
            all_cluster_sizes.append(len(vertex_list))
    
    # Calculate unrecovered deleted vertices (these become singleton clusters)
    unrecovered_deleted = deleted_vertices - total_recovered
    
    # Count pep_ident conflicts and duplicates
    clusters_with_conflicting_pep_idents = 0
    components_with_cluster_pep_ident_duplicates = 0
    
    for comp_id, comp_info in component_cluster_info.items():
        # Check for inter-cluster pep_ident duplicates
        if comp_info['num_clusters'] > 1:  # Only for multi-cluster components
            if detect_cluster_pep_ident_duplicates_in_component(comp_info, vertex_attrs):
                components_with_cluster_pep_ident_duplicates += 1
        
        # Check for intra-cluster conflicts
        for cluster_id, vertex_list in comp_info.get('clusters', {}).items():
            if detect_conflicting_pep_idents_in_cluster(vertex_list, vertex_attrs):
                clusters_with_conflicting_pep_idents += 1
    
    # Create cluster recovered distribution (after recovery)
    cluster_recovered_final = create_cluster_size_distribution(all_cluster_sizes)
    
    # Calculate unpaired/unrecovered singletons count for metrics
    unpaired_count = len(unpaired_vertices) if unpaired_vertices else 0
    total_singletons = unpaired_count + unrecovered_deleted
    
    # Include unpaired singletons AND unrecovered deleted vertices in average cluster size calculation
    # Both become singleton (size 1) clusters
    total_size_sum = sum(all_cluster_sizes) + total_singletons  
    total_cluster_count = len(all_cluster_sizes) + total_singletons 
    
    # avg_cluster_size is calculated AFTER recovery - includes recovered single-vertex clusters
    # Also includes unpaired singletons and unrecovered deleted vertices (all treated as clusters of size 1)
    avg_cluster_size = total_size_sum / total_cluster_count if total_cluster_count > 0 else 0.0
    
    # Calculate score (maximize avg_cluster_size, minimize vertices_deleted)
    # score = avg_cluster_size - λ * (vertices_deleted / total_features)
    vertices_deleted_ratio = 0.0
    if total_features_after_filter and total_features_after_filter > 0:
        vertices_deleted_ratio = deleted_vertices / total_features_after_filter
    
    # Calculate average modularity across components
    modularity_values = []
    for comp_id, comp_info in component_cluster_info.items():
        modularity = comp_info.get('modularity')
        if modularity is not None:
            modularity_values.append(modularity)
    
    avg_modularity = sum(modularity_values) / len(modularity_values) if modularity_values else 0.0
    
    # DEBUG: Print singleton and cluster count metrics
    print(f"[DEBUG] total_recovered={total_recovered}, deleted_vertices={deleted_vertices}, unrecovered_deleted={unrecovered_deleted}, total_singletons={total_singletons}, total_size_sum={total_size_sum}, total_cluster_count={total_cluster_count}")
    
    # Add singleton clusters (unpaired + unrecovered) to the distribution for metrics
    if total_singletons > 0:
        cluster_recovered_final[1] = cluster_recovered_final.get(1, 0) + total_singletons
    
    # Note: lambda will be passed in main function for scoring
    # For now, just return the metrics
    
    # SAVE recovered state BEFORE adding singletons (for clusters_recovered.json export)
    import copy
    component_cluster_info_recovered = copy.deepcopy(component_cluster_info)
    
    # Also save filtered state (in case there's no actual filtering, use this for clusters_filtered.json)
    component_cluster_info_filtered = copy.deepcopy(component_cluster_info) if not filtered_cluster_info else None
    
    # ===== PHASE 2a: Add real unpaired singletons to cluster structures (for consistency with metrics) =====
    # The metrics already count unpaired singletons in avg_cluster_size calculation
    # So the JSON files should also include them to be consistent
    # Use the existing function to add real unpaired vertices as singleton clusters
    if unpaired_vertices:
        create_singleton_clusters_for_unpaired(
            unpaired_vertices,
            component_cluster_info,
            vertex_to_cluster_map,
            vertex_attrs,
            next_vertex_id=max(vertex_to_cluster_map.keys()) + 1 if vertex_to_cluster_map else 0,
            next_component_id=max(component_cluster_info.keys()) + 1 if component_cluster_info else 0
        )
    
    # ===== PHASE 2: Save clustering results to iteration directory (if requested) =====
    iteration_dir = None
    if iteration_output_dir is not None and iteration_num is not None:
        iteration_dir = f"{iteration_output_dir}/iteration_{iteration_num}"
        try:
            os.makedirs(iteration_dir, exist_ok=True)
            
            # Export clusters (before filtering/recovery applied)
            export_clusters_to_json(vertex_to_cluster_map, component_cluster_info, vertex_attrs,
                                   f"{iteration_dir}/clusters_fraction1.0.json")
            
            # Export filtered clusters (after filtering, before recovery)
            if filtered_cluster_info:
                export_filtered_clusters_to_json(filtered_vertex_map if filtered_vertex_map else vertex_to_cluster_map,
                                                filtered_cluster_info, vertex_attrs,
                                                f"{iteration_dir}/clusters.json")
            else:
                # If no filtering, export from pre-singleton copy
                export_filtered_clusters_to_json(vertex_to_cluster_map, component_cluster_info_filtered, vertex_attrs,
                                                f"{iteration_dir}/clusters.json")
            
            # Export recovered clusters (after recovery, before unpaired singletons)
            if 'component_cluster_info_recovered' in locals():
                export_clusters_to_json(recovery_map if recovery_map else vertex_to_cluster_map, 
                                       component_cluster_info_recovered, vertex_attrs,
                                       f"{iteration_dir}/clusters_recovered.json")
            
            # Save deleted vertices report
            if deleted_list:
                save_deleted_vertices_report(deleted_list, 
                                            output_dir=f"{iteration_dir}/deleted_vertices.json",
                                            filter_method=filter_method,
                                            num_clusters_with_recoveries=0,  # Will be calculated if needed
                                            total_recovered=total_recovered)
            else:
                # Create empty deleted_vertices.json if no deletions
                with open(f"{iteration_dir}/deleted_vertices.json", 'w') as f:
                    json.dump({'summary': {'total_deleted': 0, 'total_recovered': 0}}, f, indent=2)
            
            # Save graph analysis (pre-computed, just copy to keep consistency)
            if pre_computed_analysis:
                analysis_path = f"{iteration_dir}/graph_analysis_fraction1.0.json"
                with open(analysis_path, 'w') as f:
                    json.dump(pre_computed_analysis, f, indent=2)
            
            # Save graph composition (pre-computed, just copy to keep consistency)
            if pre_computed_composition:
                composition_path = f"{iteration_dir}/graph_composition_fraction1.0.json"
                with open(composition_path, 'w') as f:
                    json.dump(pre_computed_composition, f, indent=2)
            
        except Exception as e:
            print(f"\n✗ ERROR: Failed to save iteration {iteration_num} results to {iteration_dir}")
            print(f"  Error: {str(e)}")
            raise  # Hard stop as per user specification
    
    return {
        'resolution': resolution,
        'vertices_deleted': deleted_vertices,
        'recovered_vertices': total_recovered,
        'avg_cluster_size': avg_cluster_size,
        'avg_modularity': avg_modularity,
        'vertices_deleted_ratio': vertices_deleted_ratio,
        'num_vertices': graph.vcount(),
        'num_clusters': total_cluster_count,
        'clusters_with_conflicting_pep_idents': clusters_with_conflicting_pep_idents,
        'components_with_cluster_pep_ident_duplicates': components_with_cluster_pep_ident_duplicates,
        'initial_cluster_distribution': initial_cluster_distribution,
        'filtered_cluster_distribution': filtered_cluster_distribution,
        'cluster_recovered_final': cluster_recovered_final,
        'iteration_dir': iteration_dir if iteration_dir else None
    }


def select_optimal_resolution(graph: ig.Graph, vertex_attrs: Dict,
                             method: str = 'leiden',
                             use_weights: bool = False,
                             weight_mode: str = 'inverse',
                             objective_function: str = 'CPM',
                             apply_filtering: bool = False,
                             filter_method: str = 'simplified-modularity',
                             total_features_after_filter: int = None,
                             initial_resolution: float = 0.1,
                             lambda_param: float = 0.1,
                             max_iterations: int = 15,
                             resolution_tolerance: float = 1e-5,
                             optimization_metric: str = 'combined',
                             output_iterations_path: str = None,
                             unpaired_vertices: List[Dict] = None,
                             num_cores: int = None,
                             initial_seed: int = 42,
                             mega_complexity_threshold: float = 10000,
                             iteration_output_dir: str = None,
                             pre_computed_analysis: Dict = None,
                             pre_computed_composition: Dict = None) -> Tuple[float, Dict, List, Dict, int]:
    """
    Automatically select optimal resolution parameter using Brent's method (scipy).
    
    Optimization metrics:
    - 'combined': Maximize avg_cluster_size while minimizing vertices_deleted
      Score = avg_cluster_size - λ * (vertices_deleted / total_features)
    - 'modularity': Maximize average modularity directly
    - 'avg_cluster_size': Maximize overall average cluster size (all clusters, after filtering & recovery)
    
    Args:
        graph: igraph Graph object
        vertex_attrs: Dictionary of vertex attributes
        method: Clustering method (default: leiden)
        use_weights: Whether to use edge weights
        weight_mode: Weight mode (inverse or distance)
        objective_function: Objective function for leiden
        apply_filtering: Whether to apply duplicate filtering
        filter_method: Filtering method
        total_features_after_filter: Total features for ratio calculation
        initial_resolution: Initial resolution to start optimization from
        lambda_param: Weight of deletion penalty in score function (only for 'combined' metric)
        max_iterations: Maximum number of iterations (passed to scipy, may not be fully respected)
        resolution_tolerance: Convergence tolerance for resolution parameter (default: 1e-5, lower=more precise)
        output_iterations_path: Path to save iteration history JSON
    
    Returns:
        Tuple of (optimal_resolution, best_metrics, iteration_history)
    """
    print(f"\n{'='*70}")
    print("RESOLUTION OPTIMIZATION (Brent's Method via scipy.optimize)")
    print(f"{'='*70}")
    print(f"  Initial resolution: {initial_resolution}")
    print(f"  Resolution bounds: [0.01, 2.0]")
    print(f"  Optimization metric: {optimization_metric.upper()}")
    if optimization_metric == 'combined':
        print(f"  Lambda (deletion penalty): {lambda_param}")
        print(f"  Objective: Maximize avg_cluster_size - λ * (vertices_deleted / total_features)")
    elif optimization_metric == 'modularity':
        print(f"  Objective: Maximize average modularity")
    elif optimization_metric == 'avg_cluster_size':
        print(f"  Objective: Maximize overall average cluster size (all clusters, after filtering & recovery)")
    print(f"  Convergence tolerance (xatol): {resolution_tolerance}")
    print(f"  Max iterations: {max_iterations}\n")
    
    if minimize_scalar is None:
        print("ERROR: scipy not installed. Cannot use Brent's method optimization.")
        print("Install with: pip install scipy")
        return initial_resolution, {}, []
    
    def compute_score(metrics_dict):
        """Helper to compute score from metrics based on optimization metric."""
        if optimization_metric == 'modularity':
            return metrics_dict.get('avg_modularity', -float('inf')) if metrics_dict.get('avg_modularity') is not None else -float('inf')
        elif optimization_metric == 'avg_cluster_size':
            return metrics_dict.get('avg_cluster_size', 0.0)
        else:  # 'combined'
            return metrics_dict['avg_cluster_size'] - lambda_param * metrics_dict['vertices_deleted_ratio']
    
    # Storage for iteration history and best iteration tracking
    iteration_history = []
    call_count = [0]  # Use list to allow modification in nested function
    best_iteration_num = None  # Track which iteration (1-indexed) is best
    
    def objective_function_wrapper(resolution):
        """
        Wrapper that evaluates the score at a given resolution.
        Negates the score because scipy minimizes by default.
        Tracks iteration history.
        """
        call_count[0] += 1
        iteration = call_count[0]
        
        print(f"  [{iteration:3d}] res={resolution:.4f}...", end='', flush=True)
        
        # Evaluate current resolution
        metrics = evaluate_resolution(
            graph, vertex_attrs,
            resolution=resolution,
            method=method,
            use_weights=use_weights,
            weight_mode=weight_mode,
            objective_function=objective_function,
            apply_filtering=apply_filtering,
            filter_method=filter_method,
            total_features_after_filter=total_features_after_filter,
            unpaired_vertices=unpaired_vertices,
            num_cores=num_cores,
            initial_seed=initial_seed,
            mega_complexity_threshold=mega_complexity_threshold,
            iteration_output_dir=iteration_output_dir,
            iteration_num=iteration,
            pre_computed_analysis=pre_computed_analysis,
            pre_computed_composition=pre_computed_composition
        )
        score = compute_score(metrics)
        metrics['score'] = score
        
        # Format modularity for display
        modularity_val = metrics.get('avg_modularity')
        modularity_str = f"{modularity_val:.3f}" if isinstance(modularity_val, float) else 'N/A'
        
        print(f" size={metrics['avg_cluster_size']:.3f} score={score:.4f} modularity={modularity_str}")
        
        # Track iteration
        iteration_history.append({
            'iteration': iteration,
            'resolution': resolution,
            'vertices_deleted': metrics['vertices_deleted'],
            'vertices_recovered': metrics.get('recovered_vertices', 0),
            'avg_cluster_size': metrics['avg_cluster_size'],
            'avg_modularity': metrics.get('avg_modularity'),
            'score': score,
            'num_clusters': metrics['num_clusters'],
            'clusters_with_conflicting_pep_idents': metrics.get('clusters_with_conflicting_pep_idents', 0),
            'components_with_cluster_pep_ident_duplicates': metrics.get('components_with_cluster_pep_ident_duplicates', 0),
            'initial_cluster_distribution': metrics.get('initial_cluster_distribution', {}),
            'filtered_cluster_distribution': metrics.get('filtered_cluster_distribution', {}),
            'cluster_recovered_final': metrics.get('cluster_recovered_final', {})
        })
        
        # Return negative score (scipy minimizes, we want to maximize)
        return -score
    
    # Run Brent's method optimization
    result = minimize_scalar(
        objective_function_wrapper,
        bounds=(0.01, 2.0),
        method='bounded',  # Uses Brent's method internally with bounds
        options={'xatol': resolution_tolerance}
    )
    
    # Extract optimal resolution
    best_resolution = result.x
    
    # Find best metrics from iteration history (highest score)
    best_iteration = max(iteration_history, key=lambda x: x['score'])
    best_iteration_num = best_iteration['iteration']  # Extract iteration number (1-indexed)
    
    best_metrics = {
        'resolution': best_iteration['resolution'],
        'vertices_deleted': best_iteration['vertices_deleted'],
        'recovered_vertices': best_iteration['vertices_recovered'],
        'avg_cluster_size': best_iteration['avg_cluster_size'],
        'avg_modularity': best_iteration['avg_modularity'],
        'score': best_iteration['score'],
        'num_clusters': best_iteration['num_clusters'],
        'clusters_with_conflicting_pep_idents': best_iteration['clusters_with_conflicting_pep_idents'],
        'components_with_cluster_pep_ident_duplicates': best_iteration['components_with_cluster_pep_ident_duplicates'],
        'initial_cluster_distribution': best_iteration['initial_cluster_distribution'],
        'filtered_cluster_distribution': best_iteration['filtered_cluster_distribution'],
        'cluster_recovered_final': best_iteration['cluster_recovered_final']
    }
    
    # Get clustering data if available (for reuse to avoid re-clustering)
    # NOTE: In Phase 2 refactoring, clustering data is saved to disk per iteration.
    # Returning empty dict here; final processing will load from disk if needed.
    best_clustering_data = {}  # No longer stored in iteration_history
    
    # Print summary
    print(f"\n{'='*70}")
    print("OPTIMIZATION RESULT")
    print(f"{'='*70}")
    print(f"  Best resolution: {best_resolution:.4f}")
    print(f"  Metric: {optimization_metric.upper()}")
    print(f"  Best score: {best_iteration['score']:.4f}")
    print(f"  Avg cluster size: {best_metrics['avg_cluster_size']:.3f}")
    print(f"  Avg modularity: {best_metrics.get('avg_modularity', 'N/A')}")
    print(f"  Vertices deleted: {best_metrics['vertices_deleted']}")
    print(f"  Vertices recovered: {best_metrics.get('recovered_vertices', 0)}")
    print(f"  Number of clusters: {best_metrics['num_clusters']}")
    print(f"  Function evaluations: {call_count[0]}")
    print(f"  Scipy success: {result.success}")
    print(f"  Scipy message: {result.message}")
    print(f"\n")
    
    # Save iteration history if requested
    if output_iterations_path:
        iterations_output = {
            'optimization_config': {
                'method': 'brent_scipy',
                'clustering_method': method,
                'resolution_bounds': [0.01, 2.0],
                'initial_resolution': initial_resolution,
                'optimization_metric': optimization_metric,
                'lambda': lambda_param,
                'filter_enabled': apply_filtering,
                'filter_method': filter_method if apply_filtering else None
            },
            'iterations': iteration_history
        }
        
        with open(output_iterations_path, 'w') as f:
            json.dump(iterations_output, f, indent=2)
        print(f"  Saved iteration history: {output_iterations_path}")
    
    return best_resolution, best_metrics, iteration_history, best_clustering_data, best_iteration_num


def grid_search_resolution(graph: ig.Graph, vertex_attrs: Dict,
                          method: str = 'leiden',
                          use_weights: bool = False,
                          weight_mode: str = 'inverse',
                          objective_function: str = 'CPM',
                          apply_filtering: bool = False,
                          filter_method: str = 'simplified-modularity',
                          total_features_after_filter: int = None,
                          min_resolution: float = 0.01,
                          max_resolution: float = 2.0,
                          num_points: int = 20,
                          lambda_param: float = 0.1,
                          optimization_metric: str = 'combined',
                          output_iterations_path: str = None,
                          unpaired_vertices: List[Dict] = None,
                          num_cores: int = None,
                          initial_seed: int = 42,
                          mega_complexity_threshold: float = 10000,
                          iteration_output_dir: str = None,
                          pre_computed_analysis: Dict = None,
                          pre_computed_composition: Dict = None) -> Tuple[float, Dict, List, Dict, int]:
    """
    Automatically select optimal resolution parameter using grid search.
    
    Evaluates resolution at evenly spaced points across the specified range.
    
    Optimization metrics:
    - 'combined': Maximize avg_cluster_size while minimizing vertices_deleted
      Score = avg_cluster_size - λ * (vertices_deleted / total_features)
    - 'modularity': Maximize average modularity directly
    
    Args:
        graph: igraph Graph object
        vertex_attrs: Dictionary of vertex attributes
        method: Clustering method (default: leiden)
        use_weights: Whether to use edge weights
        weight_mode: Weight mode (inverse or distance)
        objective_function: Objective function for leiden
        apply_filtering: Whether to apply duplicate filtering
        filter_method: Filtering method
        total_features_after_filter: Total features for ratio calculation
        min_resolution: Minimum resolution value to test
        max_resolution: Maximum resolution value to test
        num_points: Number of evenly spaced points to evaluate
        lambda_param: Weight of deletion penalty in score function (only for 'combined' metric)
        optimization_metric: Optimization metric ('combined' or 'modularity')
        output_iterations_path: Path to save iteration history JSON
    
    Returns:
        Tuple of (optimal_resolution, best_metrics, iteration_history)
    """
    print(f"\n{'='*70}")
    print("RESOLUTION OPTIMIZATION (Grid Search)")
    print(f"{'='*70}")
    print(f"  Resolution range: {min_resolution:.3f} to {max_resolution:.3f}")
    print(f"  Number of points: {num_points}")
    print(f"  Optimization metric: {optimization_metric.upper()}")
    if optimization_metric == 'combined':
        print(f"  Lambda (deletion penalty): {lambda_param}")
        print(f"  Objective: Maximize avg_cluster_size - λ * (vertices_deleted / total_features)")
    elif optimization_metric == 'modularity':
        print(f"  Objective: Maximize average modularity")
    elif optimization_metric == 'avg_cluster_size':
        print(f"  Objective: Maximize overall average cluster size (all clusters, after filtering & recovery)")
    print(f"\n")
    
    # Generate evenly spaced resolution values
    resolution_values = np.linspace(min_resolution, max_resolution, num_points)
    
    iteration_history = []
    best_resolution = None
    best_score = None
    best_metrics = None
    best_iteration_num = None  # Track which iteration (1-indexed) is best
    # NOTE: In Phase 2 refactoring, clustering data is saved to disk per iteration.
    # No longer storing _clustering_data in metrics; will be loaded from disk as needed.
    best_clustering_data = {}  # Empty dict for consistency with select_optimal_resolution
    
    for iteration, current_resolution in enumerate(resolution_values):
        # iteration is 0-indexed from enumerate; convert to 1-indexed for consistency
        iteration_num = iteration + 1
        print(f"  [{iteration_num:2d}/{num_points}] Evaluating resolution={current_resolution:.4f}...", end='', flush=True)
        
        # Evaluate current resolution
        metrics = evaluate_resolution(
            graph, vertex_attrs,
            resolution=current_resolution,
            method=method,
            use_weights=use_weights,
            weight_mode=weight_mode,
            objective_function=objective_function,
            apply_filtering=apply_filtering,
            filter_method=filter_method,
            total_features_after_filter=total_features_after_filter,
            unpaired_vertices=unpaired_vertices,
            num_cores=num_cores,
            initial_seed=initial_seed,
            mega_complexity_threshold=mega_complexity_threshold,
            iteration_output_dir=iteration_output_dir,
            iteration_num=iteration_num,
            pre_computed_analysis=pre_computed_analysis,
            pre_computed_composition=pre_computed_composition
        )
        
        # Calculate score based on selected optimization metric
        if optimization_metric == 'modularity':
            score = metrics.get('avg_modularity', -float('inf')) if metrics.get('avg_modularity') is not None else -float('inf')
        elif optimization_metric == 'avg_cluster_size':
            score = metrics.get('avg_cluster_size', 0.0)  # Maximize overall average cluster size
        else:  # 'combined' (default)
            score = metrics['avg_cluster_size'] - lambda_param * metrics['vertices_deleted_ratio']
        
        metrics['score'] = score
        
        # Format modularity for display
        modularity_val = metrics.get('avg_modularity')
        modularity_str = f"{modularity_val:.3f}" if isinstance(modularity_val, float) else 'N/A'
        
        print(f" size={metrics['avg_cluster_size']:.3f} deleted={metrics['vertices_deleted']} modularity={modularity_str} score={score:.4f}")
        
        # Track iteration
        iteration_history.append({
            'iteration': iteration_num,
            'resolution': current_resolution,
            'vertices_deleted': metrics['vertices_deleted'],
            'vertices_recovered': metrics.get('recovered_vertices', 0),
            'avg_cluster_size': metrics['avg_cluster_size'],
            'avg_modularity': metrics.get('avg_modularity'),
            'score': score,
            'num_clusters': metrics['num_clusters'],
            'clusters_with_conflicting_pep_idents': metrics.get('clusters_with_conflicting_pep_idents', 0),
            'components_with_cluster_pep_ident_duplicates': metrics.get('components_with_cluster_pep_ident_duplicates', 0),
            'initial_cluster_distribution': metrics.get('initial_cluster_distribution', {}),
            'filtered_cluster_distribution': metrics.get('filtered_cluster_distribution', {}),
            'cluster_recovered_final': metrics.get('cluster_recovered_final', {})
        })
        
        # Check if this is the best so far (maximize score for all metrics)
        if best_score is None or score > best_score:
            best_score = score
            best_resolution = current_resolution
            best_metrics = metrics
            best_iteration_num = iteration_num  # Track best iteration number (1-indexed)
            # NOTE: In Phase 2 refactoring, clustering data is saved to disk per iteration.
            # No longer storing _clustering_data in metrics.
            print(f"    ✓ New best!")
        else:
            print(f"")
    
    # Print summary
    print(f"\n{'='*70}")
    print("OPTIMIZATION RESULT")
    print(f"{'='*70}")
    print(f"  Best resolution: {best_resolution:.4f}")
    print(f"  Metric: {optimization_metric.upper()}")
    print(f"  Best score: {best_score:.4f}")
    print(f"  Avg cluster size: {best_metrics['avg_cluster_size']:.3f}")
    print(f"  Avg modularity: {best_metrics.get('avg_modularity', 'N/A')}")
    print(f"  Vertices deleted: {best_metrics['vertices_deleted']}")
    print(f"  Vertices recovered: {best_metrics.get('recovered_vertices', 0)}")
    print(f"  Number of clusters: {best_metrics['num_clusters']}")
    print(f"  Resolution values tested: {num_points}\n")
    
    # Save iteration history if requested
    if output_iterations_path:
        iterations_output = {
            'optimization_config': {
                'method': 'grid_search',
                'clustering_method': method,
                'min_resolution': min_resolution,
                'max_resolution': max_resolution,
                'num_points': num_points,
                'optimization_metric': optimization_metric,
                'lambda': lambda_param,
                'filter_enabled': apply_filtering,
                'filter_method': filter_method if apply_filtering else None
            },
            'iterations': iteration_history
        }
        
        with open(output_iterations_path, 'w') as f:
            json.dump(iterations_output, f, indent=2)
        print(f"  Saved iteration history: {output_iterations_path}")
    
    return best_resolution, best_metrics, iteration_history, best_clustering_data, best_iteration_num


def cluster_graph_independent_components(graph: ig.Graph, vertex_attrs: Dict, 
                                          method: str = 'louvain',
                                          use_weights: bool = False,
                                          weight_mode: str = 'inverse',
                                          resolution_parameter: float = 0.5,
                                          objective_function: str = 'CPM',
                                          random_seed: int = None) -> Tuple[Dict, Dict]:
    """
    Perform clustering on each connected component independently.
    
    Args:
        graph: igraph Graph object
        vertex_attrs: Dictionary of vertex attributes
        method: Clustering method ('louvain', 'leiden', 'walktrap', 'label_propagation', 'edge_betweenness', 'spinglass', 'cluster_optimal', 'cluster_faster_greedy')
    
    Returns:
        Tuple of (vertex_to_cluster_map, component_cluster_info)
        - vertex_to_cluster_map: {vertex_id: ('component_id', 'cluster_id')}
        - component_cluster_info: {component_id: {'vertices': [...], 'clusters': {cluster_id: [vertices]}, 'method': ..., 'modularity': ...}}
    """
    weight_note = f", weighted ({weight_mode})" if use_weights else ""
    print(f"\nClustering graph using {method} method (independent components{weight_note})...")
    
    if graph is None or graph.vcount() == 0:
        print("✗ Cannot cluster empty graph")
        return {}, {}
    
    components = graph.components()
    print(f"  Found {len(components)} connected components")
    
    vertex_to_cluster_map = {}
    component_cluster_info = {}
    
    for comp_idx, component_vertices in enumerate(components):
        if len(component_vertices) == 1:
            # Single vertex component - assign to its own cluster
            vertex = component_vertices[0]
            vertex_to_cluster_map[vertex] = (comp_idx, 0)
            component_cluster_info[comp_idx] = {
                'vertices': component_vertices,
                'clusters': {0: component_vertices},
                'method': method,
                'modularity': None,
                'num_clusters': 1,
                'details': f"Single vertex component"
            }
            continue
        
        # Extract subgraph for this component
        subgraph = graph.induced_subgraph(component_vertices)
        
        # Perform clustering based on method
        try:
            # Reset random seed before clustering if provided (ensures deterministic results)
            if random_seed is not None:
                import random
                random.seed(random_seed)
                # Note: This resets Python's RNG; igraph will use it via ig.set_random_number_generator()
            
            weights = None
            if use_weights:
                weights = _compute_cluster_edge_weights(subgraph, weight_mode)
                subgraph.es['weight'] = weights

            if method.lower() == 'louvain':
                clustering = subgraph.community_multilevel(weights=weights) if use_weights else subgraph.community_multilevel()
            elif method.lower() == 'leiden':
                clustering = subgraph.community_leiden(weights=weights, objective_function=objective_function, n_iterations = 100, resolution = resolution_parameter, beta=0.1) if use_weights else subgraph.community_leiden(objective_function=objective_function, n_iterations = 100, resolution = resolution_parameter, beta=0.1)
            elif method.lower() == 'walktrap':
                clustering = subgraph.community_walktrap(weights=weights).as_clustering() if use_weights else subgraph.community_walktrap().as_clustering()
            elif method.lower() == 'label_propagation':
                clustering = subgraph.community_label_propagation(weights=weights) if use_weights else subgraph.community_label_propagation()
            elif method.lower() == 'edge_betweenness':
                betweenness_weights = weights
                if use_weights and weight_mode == 'inverse':
                    betweenness_weights = _compute_cluster_edge_weights(subgraph, 'distance')
                clustering = subgraph.community_edge_betweenness(weights=betweenness_weights).as_clustering() if use_weights else subgraph.community_edge_betweenness().as_clustering()
            elif method.lower() == 'spinglass':
                clustering = subgraph.community_spinglass(weights=weights) if use_weights else subgraph.community_spinglass()
            elif method.lower() == 'cluster_optimal':
                # Optimal modularity (works only on undirected graphs)
                if subgraph.is_directed():
                    subgraph = subgraph.as_undirected()
                clustering = subgraph.community_optimal_modularity()
            elif method.lower() == 'cluster_faster_greedy':
                # Fast greedy algorithm (works only on undirected graphs)
                if subgraph.is_directed():
                    subgraph = subgraph.as_undirected()
                clustering = subgraph.community_fastgreedy(weights=weights).as_clustering() if use_weights else subgraph.community_fastgreedy().as_clustering()
            else:
                print(f"  ✗ Unknown clustering method: {method}. Using Louvain.")
                clustering = subgraph.community_multilevel()
            
            # Extract modularity
            try:
                modularity = subgraph.modularity(clustering)
            except:
                modularity = None
            
            # Map cluster IDs back to original vertex IDs
            clusters_dict = {}
            for cluster_id, members in enumerate(clustering):
                original_vertex_ids = [component_vertices[local_idx] for local_idx in members]
                clusters_dict[cluster_id] = original_vertex_ids
                
                # Map each vertex to its cluster
                for vertex_id in original_vertex_ids:
                    vertex_to_cluster_map[vertex_id] = (comp_idx, cluster_id)
            
            component_cluster_info[comp_idx] = {
                'vertices': component_vertices,
                'clusters': clusters_dict,
                'method': method,
                'modularity': float(modularity) if modularity is not None else None,
                'num_clusters': len(clusters_dict),
                'details': f"{len(clusters_dict)} clusters detected, modularity={modularity:.3f}" if modularity is not None else f"{len(clusters_dict)} clusters detected"
            }
            
        except Exception as e:
            print(f"  ⚠ Error clustering component {comp_idx}: {e}")
            # Fall back to single cluster
            vertex_to_cluster_map[component_vertices[0]] = (comp_idx, 0)
            component_cluster_info[comp_idx] = {
                'vertices': component_vertices,
                'clusters': {0: component_vertices},
                'method': method,
                'modularity': None,
                'num_clusters': 1,
                'details': f"Clustering failed, using single cluster"
            }
    
    # Print clustering summary
    total_clusters = sum(info['num_clusters'] for info in component_cluster_info.values())
    print(f"  ✓ Clustering complete: {total_clusters} clusters across {len(components)} components")
    
    for comp_idx in sorted(component_cluster_info.keys())[:10]:  # Show first 10
        info = component_cluster_info[comp_idx]
        print(f"    Component {comp_idx}: {info['details']}")
    
    if len(component_cluster_info) > 10:
        print(f"    ... and {len(component_cluster_info) - 10} more components")
    
    return vertex_to_cluster_map, component_cluster_info


def export_clusters_to_json(vertex_to_cluster_map: Dict, component_cluster_info: Dict, 
                           vertex_attrs: Dict, output_path: str):
    """
    Export clustering results to JSON file.
    
    Args:
        vertex_to_cluster_map: {vertex_id: (component_id, cluster_id)}
        component_cluster_info: {component_id: {...}}
        vertex_attrs: {vertex_id: {filename, feature_idx, label, ...}}
        output_path: Path to save JSON file
    """
    print(f"\nExporting clusters to JSON: {output_path}")
    
    # Build cluster-centric view
    clusters_output = {}
    
    for comp_idx, info in component_cluster_info.items():
        comp_key = f"component_{comp_idx}"
        clusters_output[comp_key] = {
            'component_id': comp_idx,
            'method': info['method'],
            'num_clusters': info['num_clusters'],
            'modularity': info['modularity'],
            'details': info['details'],
            'clusters': {}
        }
        
        # For each cluster in this component
        for cluster_id, vertex_ids in info['clusters'].items():
            cluster_key = f"cluster_{cluster_id}"
            cluster_vertices = []
            
            for vertex_id in vertex_ids:
                vertex_info = vertex_attrs[vertex_id]
                cluster_vertices.append({
                    'vertex_id': vertex_id,
                    'filename': vertex_info['filename'],
                    'feature_idx': vertex_info['feature_idx'],
                    'label': vertex_info['label'],
                    'x_center': vertex_info.get('x_center', 0.0),
                    'y_center': vertex_info.get('y_center', 0.0),
                    'intensity': vertex_info.get('intensity', 0.0),
                    'openms_fid': vertex_info.get('openms_fid', ''),
                    'charge': vertex_info.get('charge', 0),
                    'pep_ident': vertex_info.get('pep_ident', []),
                    'prot_ident': vertex_info.get('prot_ident', []),
                    'ms2_scans': vertex_info.get('ms2_scans', '')
                })
            
            clusters_output[comp_key]['clusters'][cluster_key] = {
                'cluster_id': cluster_id,
                'num_vertices': len(vertex_ids),
                'pep_ident': aggregate_cluster_pep_idents(cluster_vertices),
                'prot_ident': aggregate_cluster_prot_idents(cluster_vertices),
                'vertices': cluster_vertices
            }
    
    # Save to JSON
    with open(output_path, 'w') as f:
        json.dump(clusters_output, f, indent=2)
    
    print(f"✓ Saved clusters to JSON: {output_path}")
    
    # Print summary
    total_clusters = sum(len(info['clusters']) for info in clusters_output.values())
    total_vertices = sum(len(cluster['vertices']) 
                        for comp in clusters_output.values() 
                        for cluster in comp['clusters'].values())
    print(f"  Summary: {len(clusters_output)} components, {total_clusters} total clusters, {total_vertices} vertices")


def _filter_component_worker(component_data: Dict) -> Dict:
    """
    Worker function to filter a single component independently.
    Runs in a separate process (ProcessPoolExecutor) to avoid GIL.
    
    Args:
        component_data: Dictionary containing:
            - comp_idx: Component index
            - comp_vertices: List of vertex IDs in component
            - comp_clusters: Dict of cluster_id -> vertex_ids
            - vertex_attrs: Vertex attributes (subset needed for this component)
            - vertex_to_cluster_map: Mapping (subset for this component)
            - info: Component info dict
            - filter_method: Filtering method
            - mega_complexity_threshold: Threshold for using heuristic
            NOTE: graph is accessed via global (inherited from parent process)
    
    Returns:
        Dictionary with results for this component
    """
    try:
        import os
        import time
        import sys
        worker_pid = os.getpid()
        
        # Log available CPUs for debugging
        available_cpus = os.sched_getaffinity(0)
        print(f"[WORKER-TASK {worker_pid}] Available CPUs: {available_cpus} (count: {len(available_cpus)})", file=sys.stderr, flush=True)
        
        # Get actual running CPU if possible (Linux-specific)
        try:
            actual_cpu_start = os.sched_getaffinity(0)
            cpu_list_start = sorted(list(actual_cpu_start))
            print(f"[WORKER-TASK {worker_pid}] Assigned CPUs at START: {cpu_list_start}", file=sys.stderr, flush=True)
        except:
            cpu_list_start = None
        
        start_time = time.time()
        
        comp_idx = component_data['comp_idx']
        comp_vertices = component_data['comp_vertices']
        comp_clusters = component_data['comp_clusters']
        vertex_attrs = component_data['vertex_attrs']
        vertex_to_cluster_map = component_data['vertex_to_cluster_map']
        info = component_data['info']
        filter_method = component_data['filter_method']
        mega_complexity_threshold = component_data.get('mega_complexity_threshold', 10000)
        
        # Skip components with single vertices
        if len(comp_vertices) <= 1:
            elapsed = time.time() - start_time
            return {
                'comp_idx': comp_idx,
                'status': 'skipped',
                'filtered_cluster_info': info,
                'deleted_vertices': [],
                'combo_estimate': 0,
                '__worker_pid': worker_pid,
                '__elapsed_ms': int(elapsed * 1000)
            }
        
        # Find clusters with duplicate files
        clusters_with_duplicates = {}
        all_non_duplicate_vertices = set()
        
        for cluster_id, vertex_ids in comp_clusters.items():
            file_vertex_map = {}
            for vertex_id in vertex_ids:
                filename = vertex_attrs[vertex_id]['filename']
                if filename not in file_vertex_map:
                    file_vertex_map[filename] = []
                file_vertex_map[filename].append(vertex_id)
            
            has_duplicates = any(len(vids) > 1 for vids in file_vertex_map.values())
            if has_duplicates:
                clusters_with_duplicates[cluster_id] = file_vertex_map
            else:
                all_non_duplicate_vertices.update(vertex_ids)
        
        # If no duplicates, keep as-is
        if not clusters_with_duplicates:
            return {
                'comp_idx': comp_idx,
                'status': 'skipped',
                'filtered_cluster_info': info,
                'deleted_vertices': [],
                'combo_estimate': 0
            }
        
        # Retrieve graph from module globals (set by parent via fork, accessible since added to global list)
        graph = globals().get('graph', None)
        
        # DEBUG: Alert if graph is None
        if graph is None:
            print(f"[COMP {comp_idx}] ⚠️  WARNING: graph is None! Heuristic and induced_subgraph will fall back", flush=True)
        
        # Build component_cluster_info for this component from work_item data
        # This is needed for heuristic calculation (internal vs external degree)
        # Even though global component_cluster_info was deleted to save memory, we can reconstruct it here
        component_cluster_info_for_heuristic = {
            comp_idx: {
                'clusters': comp_clusters,
                'vertices': comp_vertices
            }
        }
        
        # Generate candidates (returns list of all candidates + complexity_estimate + use_heuristic flag)
        candidates, combo_estimate, use_heuristic_flag = _generate_filter_candidates(comp_idx, clusters_with_duplicates, vertex_attrs, all_non_duplicate_vertices, graph, component_cluster_info_for_heuristic, mega_complexity_threshold)
        
        # Using graph.neighbors() directly for neighbor lookups
        # graph is available via global inherited from parent
        comp_vertex_set = set(comp_vertices)
        
        # Decision was made during _generate_filter_candidates based on complexity >= mega_complexity_threshold
        # If heuristic was used, candidates[0] contains the heuristic-selected vertices
        # If not, candidates contains all possible combinations to evaluate
        
        if use_heuristic_flag:
            print(f"[COMP {comp_idx}] Using heuristic (estimated complexity: {combo_estimate:.2e} >= {mega_complexity_threshold:.2e})", flush=True)
            # Heuristic candidate was already built inside _generate_filter_candidates
            # Use it directly without further evaluation
            best_candidate = candidates[0] if candidates else set(comp_vertices)
            best_modularity = None  # Heuristic doesn't compute modularity
            candidates_evaluated = 0
            deleted_vertices_list = []
        else:
            # Use exhaustive streaming evaluation
            best_candidate = None
            best_modularity = None
            candidates_evaluated = 0
            deleted_vertices_list = []
            
            # Track which scoring paths are used
            scoring_path_counts = {'induced_subgraph': 0, 'fallback_from_exception': 0, 'degree_fallback_no_graph': 0, 
                                   'modularity': 0, 'modularity_fallback': 0, 'modularity_no_graph_fallback': 0}
            
            # DEBUG: Check cluster count
            print(f"[COMP {comp_idx}] Starting streaming evaluation: {len(clusters_with_duplicates)} clusters with duplicates, combo_est={combo_estimate:.2e}", flush=True)
            
            for candidate in candidates:
                candidates_evaluated += 1
                
                # Log progress periodically
                if candidates_evaluated % 100000 == 0:
                    print(f"[COMP {comp_idx}] Evaluated {candidates_evaluated} candidates... best={best_modularity}", flush=True)
                
                # Extract vertices in this candidate that are also in component
                test_vertices = [v for v in comp_vertices if v in candidate]
                if len(test_vertices) <= 1:
                    continue
                
                # Score this ONE candidate
                if len(clusters_with_duplicates) == 1:
                    # Single cluster: score by total degree using graph
                    if graph is not None:
                        total_degree = sum(graph.degree(v) for v in test_vertices)
                    else:
                        # Fallback: use candidate size as proxy
                        total_degree = len(test_vertices)
                    score = total_degree
                else:
                    # Multiple clusters: use weighted degree or modularity depending on filter_method
                    if filter_method == 'simplified-modularity':
                        # Create induced subgraph with only the candidate vertices (match linear version)
                        scoring_path = None
                        if graph is not None:
                            try:
                                test_subgraph = graph.induced_subgraph(test_vertices)
                                scoring_path = 'induced_subgraph'
                                weighted_score = 0
                                for v_idx, v in enumerate(test_vertices):
                                    if v in vertex_to_cluster_map:
                                        v_cluster = vertex_to_cluster_map[v][1]
                                        for neighbor_idx in test_subgraph.neighbors(v_idx):
                                            neighbor_v = test_vertices[neighbor_idx]
                                            if neighbor_v in vertex_to_cluster_map:
                                                neighbor_cluster = vertex_to_cluster_map[neighbor_v][1]
                                                if neighbor_cluster == v_cluster:
                                                    weighted_score += 1  # Internal edge
                                                else:
                                                    weighted_score -= 1  # External edge
                                score = weighted_score
                            except Exception as e:
                                # Fallback if induced_subgraph fails: use simple degree
                                scoring_path = f'degree_fallback: {type(e).__name__}'
                                score = len(test_vertices)  # Simple fallback: candidate size
                        else:
                            # Graph should always be available; use simple degree fallback
                            scoring_path = 'degree_fallback_no_graph'
                            score = len(test_vertices)  # Simple fallback: candidate size
                        
                        # Track which scoring path was used
                        if scoring_path in scoring_path_counts:
                            scoring_path_counts[scoring_path] += 1
                    
                    elif filter_method == 'modularity':
                        # Proper modularity: calculate modularity while respecting cluster boundaries
                        scoring_path = None
                        if graph is not None:
                            try:
                                test_subgraph = graph.induced_subgraph(test_vertices)
                                scoring_path = 'modularity'
                                
                                # Calculate modularity with respect to cluster assignment
                                membership = []
                                cluster_map_local = {}
                                for v in test_vertices:
                                    if v in vertex_to_cluster_map:
                                        cluster_id = vertex_to_cluster_map[v][1]
                                        if cluster_id not in cluster_map_local:
                                            cluster_map_local[cluster_id] = len(cluster_map_local)
                                        membership.append(cluster_map_local[cluster_id])
                                    else:
                                        membership.append(0)
                                
                                score = test_subgraph.modularity(membership)
                            except Exception as e:
                                # Fallback to simple degree if modularity fails
                                scoring_path = f'modularity_degree_fallback: {type(e).__name__}'
                                score = len(test_vertices)  # Simple fallback: candidate size
                        else:
                            # Graph should always be available; use simple degree fallback
                            scoring_path = 'modularity_graph_neighbors_fallback'
                            score = len(test_vertices)  # Simple fallback: candidate size
                        
                        # Track which scoring path was used
                        if scoring_path in scoring_path_counts:
                            scoring_path_counts[scoring_path] += 1
                    else:
                        # Fallback to degree for unknown filter methods
                        if graph is not None:
                            score = sum(graph.degree(v) for v in test_vertices)
                        else:
                            score = len(test_vertices)
                
                # Update best if this is better
                if best_modularity is None or score > best_modularity:
                    best_modularity = score
                    best_candidate = candidate
                    if candidates_evaluated % 100000 == 0 or candidates_evaluated <= 10:
                        print(f"[COMP {comp_idx}] New best score: {best_modularity} (after {candidates_evaluated} evals)", flush=True)
                
                # Candidate is immediately garbage-collected when loop moves to next iteration
            
            print(f"[COMP {comp_idx}] ✓ Streaming complete: evaluated {candidates_evaluated} candidates, best={best_modularity}", flush=True)
            # Log scoring path breakdown
            if candidates_evaluated > 0:
                print(f"[COMP {comp_idx}]   Scoring paths used: {scoring_path_counts}", flush=True)
        
        # Record deleted vertices (same for both heuristic and exhaustive)
        filtered_clusters = {}
        
        if best_candidate:
            for vertex_id in comp_vertices:
                if vertex_id not in best_candidate:
                    comp_id, cluster_id = vertex_to_cluster_map[vertex_id]
                    reason_text = 'Duplicate file vertex - filtered via heuristic' if use_heuristic_flag else 'Duplicate file vertex - filtered for best modularity'
                    deleted_vertices_list.append({
                        'component_id': comp_id,
                        'cluster_id': cluster_id,
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
                        'best_modularity': float(best_modularity) if best_modularity else None
                    })
            
            # Rebuild cluster info
            for cluster_id, vertex_ids in comp_clusters.items():
                filtered_vertices = [v for v in vertex_ids if v in best_candidate]
                if filtered_vertices:
                    filtered_clusters[cluster_id] = filtered_vertices
            
            filtered_vertices_all = list(best_candidate)
            filtered_cluster_info = {
                'vertices': filtered_vertices_all,
                'clusters': filtered_clusters,
                'method': info['method'],
                'modularity': info.get('modularity'),
                'num_clusters': len(filtered_clusters),
                'details': f"{len(filtered_clusters)} clusters (filtered from {len(comp_clusters)})"
            }
        else:
            filtered_cluster_info = info
        
        elapsed = time.time() - start_time
        
        # VERIFY CPU AFFINITY HELD THROUGHOUT EXECUTION
        try:
            cpu_list_end = sorted(list(os.sched_getaffinity(0)))
            if cpu_list_start and cpu_list_end != cpu_list_start:
                print(f"[WORKER-TASK {worker_pid}] ⚠️  AFFINITY CHANGED during task!", file=sys.stderr, flush=True)
                print(f"[WORKER-TASK {worker_pid}]    START: {cpu_list_start}", file=sys.stderr, flush=True)
                print(f"[WORKER-TASK {worker_pid}]    END:   {cpu_list_end}", file=sys.stderr, flush=True)
            else:
                print(f"[WORKER-TASK {worker_pid}] ✓ Affinity held: {cpu_list_end} (duration: {elapsed:.2f}s)", file=sys.stderr, flush=True)
        except:
            pass
        
        return {
            'comp_idx': comp_idx,
            'status': 'processed',
            'filtered_cluster_info': filtered_cluster_info,
            'deleted_vertices': deleted_vertices_list,
            'combo_estimate': combo_estimate,
            '__worker_pid': worker_pid,
            '__elapsed_ms': int(elapsed * 1000)
        }
    
    except Exception as e:
        import traceback
        elapsed = time.time() - start_time
        # Return error status without crashing
        return {
            'comp_idx': component_data.get('comp_idx', -1),
            'status': 'error',
            'error': str(e),
            'filtered_cluster_info': component_data.get('info'),
            'deleted_vertices': [],
            'combo_estimate': -1,
            '__worker_pid': os.getpid(),
            '__elapsed_ms': int(elapsed * 1000)
        }


def _init_worker(worker_index, num_cpus):
    """
    Initializer function called once per worker process in the Pool.
    Pins this worker process to a specific CPU for true multi-core parallelism.
    
    Global data structures (edge_index, vertex_attrs, etc.) are inherited 
    from parent process via fork/copy-on-write - NO need to pass via initargs!
    
    Args:
        worker_index: Unique index for this worker (0 to num_cpus-1)
        num_cpus: Total number of CPUs to distribute workers across
    """
    import os
    import sys
    pid = os.getpid()
    
    # Check what CPUs are ACTUALLY available to this worker
    available_before = os.sched_getaffinity(0)
    print(f"[WORKER {pid}] Available CPUs BEFORE pinning: {available_before}", file=sys.stderr, flush=True)
    
    # Pin to specific CPU using worker index (0-based, no collisions)
    cpu_to_pin = worker_index % num_cpus
    
    try:
        # Try to pin to CPU
        os.sched_setaffinity(0, {cpu_to_pin})
        # Verify it worked
        affinity = os.sched_getaffinity(0)
        print(f"[WORKER {pid}] Pinned to CPU {cpu_to_pin}, affinity AFTER: {affinity}", file=sys.stderr, flush=True)
    except Exception as e:
        print(f"[WORKER {pid}] CPU pinning failed: {e}", file=sys.stderr, flush=True)


def _lazy_worker_with_self_batching(metadata_stream):
    """
    Module-level worker function for TRUE LAZY LOADING + WORKER-SIDE BATCHING.
    
    Receives an iterator of lightweight metadata dictionaries and:
    1. Builds full work_items on-demand (just-in-time)
    2. Accumulates work_items up to COMBINATIONS_BATCH_LIMIT combinations
    3. Processes batch when full
    4. Returns all results
    
    This ensures the orchestrator never loads all work_items into memory.
    Worker globals (edge_index, vertex_attrs, etc.) are set by _init_worker.
    """
    import os
    import gc
    
    # Access globals set by _init_worker in this worker's process
    global edge_index, vertex_attrs, vertex_to_cluster_map, filter_method, component_cluster_info, graph
    
    # ===== WORKER TIMING =====
    worker_timing = {}
    worker_timing['worker_start'] = time.time()
    
    results = []
    current_batch = []
    current_complexity = 0
    worker_pid = os.getpid()
    batch_count = 0
    
    def build_work_item_from_metadata(metadata):
        """Build a full work_item from lightweight metadata on-demand"""
        comp_idx = metadata['comp_idx']
        combo_est = metadata['combo_estimate']
        comp_vertices = metadata['comp_vertices']
        info = metadata['info']
        comp_clusters = info['clusters']
        
        # Removed: component_edges building
        # Using graph.neighbors() directly in worker functions
        
        # Create component-scoped copies
        comp_vertex_attrs = {v: vertex_attrs[v] for v in comp_vertices}
        comp_vertex_to_cluster_map = {v: vertex_to_cluster_map[v] for v in comp_vertices}
        
        # Estimate memory usage for this component
        # Each work_item contains: 
        #   - comp_vertices list: ~8 bytes per vertex
        #   - comp_clusters dict: ~50 bytes per entry
        #   - comp_vertex_attrs dict: ~100 bytes per entry
        comp_vertex_count = len(comp_vertices)
        comp_cluster_count = len(comp_clusters)
        estimated_item_size_bytes = (
            comp_vertex_count * 8 +
            comp_cluster_count * 50 +
            comp_vertex_count * 100
        )
        
        # Log large components
        if estimated_item_size_bytes > 50_000_000:  # >50MB
            print(f"[PID {worker_pid}] ⚠ LARGE COMPONENT DETECTED: comp_idx={comp_idx}, "
                  f"vertices={comp_vertex_count}, clusters={comp_cluster_count}, "
                  f"est_size={estimated_item_size_bytes / 1_000_000:.1f}MB", flush=True)
        
        return {
            'comp_idx': comp_idx,
            'comp_vertices': comp_vertices,
            'comp_clusters': comp_clusters,
            'vertex_attrs': comp_vertex_attrs,
            'vertex_to_cluster_map': comp_vertex_to_cluster_map,
            'info': info,
            'filter_method': filter_method,
            'combo_estimate': combo_est,
            'mega_complexity_threshold': metadata.get('mega_complexity_threshold', 10000),
            'will_use_heuristic': metadata.get('will_use_heuristic', False),
            '__estimated_size_bytes': estimated_item_size_bytes
        }
    
    # Process metadata stream with JUST-IN-TIME work_item building
    # Build work_items on-demand only when adding to current batch
    # This avoids holding all 2731 large work_item dicts in memory simultaneously
    metadata_list = list(metadata_stream)  # Convert to list to find max sizes
    print(f"[PID {worker_pid}] Processing {len(metadata_list)} metadata items", flush=True)
    
    # Find max complexity in this worker's chunk
    if metadata_list:
        max_combo = max(m.get('combo_estimate', 0) for m in metadata_list)
        print(f"[PID {worker_pid}] Max complexity in chunk: {max_combo:.2e}", flush=True)
    
    worker_timing['build_start'] = time.time()
    
    # ===== BUILD RESULTS BY BATCHING - RETURN LIST (PICKLABLE) =====
    # Process metadata in batches to control memory usage
    # This is critical for workers processing 1000+ components
    
    results = []  # Accumulate results from all batches
    batch_count = 0
    batch_times = []
    total_metadata_processed = 0
    
    current_batch = []
    current_complexity = 0
    worker_timing['processing_start'] = time.time()
    
    for idx, metadata in enumerate(metadata_list):
        # Build work_item JUST-IN-TIME (when we know it's needed for batching)
        combo_est = metadata.get('combo_estimate', 0)
        
        # Log progress for tracking
        if idx % 100 == 0 and idx > 0:
            print(f"[PID {worker_pid}] Processing work_items: {idx}/{len(metadata_list)}, "
                  f"current_batch_size={len(current_batch)}, current_complexity={current_complexity:.2e}", 
                  flush=True)
        
        # Would adding this metadata exceed the batch limit?
        if current_complexity + combo_est > COMBINATIONS_BATCH_LIMIT and current_batch:
            # Process current batch and start a new one (WITHOUT building all work_items)
            batch_count += 1
            batch_size = len(current_batch)
            batch_start = time.time()
            print(f"[PID {worker_pid}] Processing batch {batch_count}: {batch_size} components, "
                  f"complexity={current_complexity:.2e}", flush=True)
            # Build only the work_items in THIS batch, then process immediately
            batch_work_items = [build_work_item_from_metadata(m) for m in current_batch]
            batch_results = [_filter_component_worker(item) for item in batch_work_items]
            batch_elapsed = time.time() - batch_start
            batch_times.append(batch_elapsed)
            for result in batch_results:
                result['__worker_pid'] = worker_pid
                result['__batch_num'] = batch_count
                result['__batch_elapsed_ms'] = int(batch_elapsed * 1000)
                results.append(result)  # Append to result list
            print(f"[PID {worker_pid}] ✓ Batch {batch_count} complete in {batch_elapsed:.2f}s; "
                  f"collected {len(batch_results)} results", flush=True)
            # Clear batch and metadata to free memory
            current_batch = []
            batch_work_items = []  # Explicit free
            batch_results = []      # Explicit free
            current_complexity = 0
            gc.collect()
        
        # Add this metadata to current batch (store only metadata, not built work_item yet)
        current_batch.append(metadata)
        current_complexity += combo_est
        total_metadata_processed += 1
    
    # Process final batch
    if current_batch:
        batch_count += 1
        batch_size = len(current_batch)
        batch_start = time.time()
        print(f"[PID {worker_pid}] Processing final batch {batch_count}: {batch_size} components, "
              f"complexity={current_complexity:.2e}", flush=True)
        # Build only the work_items in final batch, then process
        batch_work_items = [build_work_item_from_metadata(m) for m in current_batch]
        batch_results = [_filter_component_worker(item) for item in batch_work_items]
        batch_elapsed = time.time() - batch_start
        batch_times.append(batch_elapsed)
        for result in batch_results:
            result['__worker_pid'] = worker_pid
            result['__batch_num'] = batch_count
            result['__batch_elapsed_ms'] = int(batch_elapsed * 1000)
            results.append(result)  # Append to result list
        print(f"[PID {worker_pid}] ✓ Final batch {batch_count} complete in {batch_elapsed:.2f}s; "
              f"collected {len(batch_results)} results", flush=True)
        current_batch = []
        batch_work_items = []  # Explicit free
        batch_results = []      # Explicit free
        gc.collect()
    
    worker_timing['processing_end'] = time.time()
    
    # ===== FREE LARGE GLOBALS AFTER ALL PROCESSING =====
    # NOW we can safely delete globals since build_work_item_from_metadata won't be called
    print(f"[PID {worker_pid}] Freeing global data structures to recover memory...", flush=True)
    worker_timing['free_start'] = time.time()
    try:
        del edge_index
        del vertex_attrs
        del vertex_to_cluster_map
        del component_cluster_info
    except:
        pass
    gc.collect()
    worker_timing['free_end'] = time.time()
    print(f"[PID {worker_pid}] ✓ Freed globals in {worker_timing['free_end'] - worker_timing['free_start']:.2f}s", flush=True)
    
    worker_timing['build_end'] = time.time()
    print(f"[PID {worker_pid}] ✓ Built and processed {total_metadata_processed} work_items", flush=True)
    
    # ===== WORKER TIMING SUMMARY =====
    worker_build_time = worker_timing['build_end'] - worker_timing['build_start']
    worker_free_time = worker_timing['free_end'] - worker_timing['free_start']
    worker_proc_time = worker_timing['processing_end'] - worker_timing['processing_start']
    worker_total_time = worker_timing['processing_end'] - worker_timing['worker_start']
    
    # Calculate batch time stats
    if batch_times:
        avg_batch_time = sum(batch_times) / len(batch_times)
        max_batch_time = max(batch_times)
        min_batch_time = min(batch_times)
    else:
        avg_batch_time = max_batch_time = min_batch_time = 0
    
    print(f"\n[PID {worker_pid}] {'='*60}", flush=True)
    print(f"[PID {worker_pid}] WORKER TIMING SUMMARY", flush=True)
    print(f"[PID {worker_pid}] {'='*60}", flush=True)
    print(f"[PID {worker_pid}] Building work_items........... {worker_build_time:>8.2f}s", flush=True)
    print(f"[PID {worker_pid}] Freeing globals............... {worker_free_time:>8.2f}s", flush=True)
    print(f"[PID {worker_pid}] Processing {total_metadata_processed} items in {batch_count} batches:", flush=True)
    print(f"[PID {worker_pid}]   - Avg batch time........... {avg_batch_time:>8.2f}s", flush=True)
    print(f"[PID {worker_pid}]   - Min batch time........... {min_batch_time:>8.2f}s", flush=True)
    print(f"[PID {worker_pid}]   - Max batch time........... {max_batch_time:>8.2f}s", flush=True)
    print(f"[PID {worker_pid}]   - Total batch time......... {worker_proc_time:>8.2f}s", flush=True)
    print(f"[PID {worker_pid}] TOTAL WORKER TIME............ {worker_total_time:>8.2f}s", flush=True)
    print(f"[PID {worker_pid}] {'='*60}\n", flush=True)
    
    # Return accumulated results list
    return results


# REMOVED: _evaluate_candidates_batch - dead code (never called, replaced by streaming evaluation in _filter_component_worker)



    """
    Worker function to process a BATCH of components (reduces pickling overhead).
    Called by each process in the pool - processes one batch sequentially.
    Batches are smaller (3x cores) to improve load balancing while reducing IPC.
    
    Args:
        components_batch: List of component work items
    
    Returns:
        List of result dictionaries from _filter_component_worker for each component in batch
    """
    import os
    worker_pid = os.getpid()
    
    batch_results = []
    
    for component_data in components_batch:
        result = _filter_component_worker(component_data)
        result['__worker_pid'] = worker_pid
        batch_results.append(result)
    
    return batch_results


# ============================================================================
# BATCH CONFIGURATION FOR ORCHESTRATOR-LEVEL COMPONENT BATCHING
# Components are grouped into batches where sum(complexity) <= limit
# This ensures bounded RAM usage across all workers processing work simultaneously
# 50,000 combinations per batch = ~250MB-1GB per worker (depends on candidate sizes)
COMBINATIONS_BATCH_LIMIT = 50000  # Increased to 5 million for worker-side batching (reduce IPC overhead)
# ============================================================================


def _batch_components_by_complexity_indices(complexity_estimates: List[Tuple]) -> List[List[int]]:
    """
    Group component INDICES into batches where total complexity per batch <= COMBINATIONS_BATCH_LIMIT.
    Lightweight version that works with (comp_idx, complexity) tuples instead of full work items.
    
    Uses greedy bin-packing: iterate through components, add to current batch until
    adding the next component would exceed the limit, then start a new batch.
    
    Args:
        complexity_estimates: List of (comp_idx, complexity_estimate) tuples
        
    Returns:
        List of batches, where each batch is a list of component indices
    """
    batches = []
    current_batch = []
    current_complexity = 0
    
    for comp_idx, complexity in complexity_estimates:
        # If component itself is larger than limit, give it its own batch
        if complexity > COMBINATIONS_BATCH_LIMIT:
            # First, finish current batch if not empty
            if current_batch:
                batches.append(current_batch)
                current_batch = []
                current_complexity = 0
            # Add large component alone
            batches.append([comp_idx])
        else:
            # Can this component fit in current batch?
            if current_complexity + complexity <= COMBINATIONS_BATCH_LIMIT:
                current_batch.append(comp_idx)
                current_complexity += complexity
            else:
                # Start new batch
                if current_batch:
                    batches.append(current_batch)
                current_batch = [comp_idx]
                current_complexity = complexity
    
    # Don't forget final batch
    if current_batch:
        batches.append(current_batch)
    
    return batches


def _batch_components_by_complexity(complexity_estimates: List[Tuple]) -> List[List]:
    """
    Group components into batches where total complexity per batch <= COMBINATIONS_BATCH_LIMIT.
    
    Uses greedy bin-packing: iterate through components, add to current batch until
    adding the next component would exceed the limit, then start a new batch.
    
    Args:
        complexity_estimates: List of (work_item, complexity_estimate) tuples
        
    Returns:
        List of batches, where each batch is a list of work_items
    """
    batches = []
    current_batch = []
    current_complexity = 0
    
    for work_item, complexity in complexity_estimates:
        # If component itself is larger than limit, give it its own batch
        if complexity > COMBINATIONS_BATCH_LIMIT:
            # First, finish current batch if not empty
            if current_batch:
                batches.append(current_batch)
                current_batch = []
                current_complexity = 0
            # Add large component alone
            batches.append([work_item])
        else:
            # Can this component fit in current batch?
            if current_complexity + complexity <= COMBINATIONS_BATCH_LIMIT:
                current_batch.append(work_item)
                current_complexity += complexity
            else:
                # Start new batch
                if current_batch:
                    batches.append(current_batch)
                current_batch = [work_item]
                current_complexity = complexity
    
    # Don't forget final batch
    if current_batch:
        batches.append(current_batch)
    
    return batches


def filter_duplicate_file_vertices(graph: ig.Graph, vertex_to_cluster_map: Dict, 
                                   component_cluster_info: Dict, vertex_attrs: Dict,
                                   filter_method: str = 'simplified-modularity',
                                   num_cores: int = None,
                                   use_multiprocessing: bool = True,
                                   mega_complexity_threshold: float = 10000) -> Tuple[Dict, Dict, List]:
    """
    Filter duplicate vertices from the same file within clusters (parallelized or serial).
    For each cluster, keep only one vertex per file - the one that maximizes connectivity/modularity.
    Processes independent components in parallel across multiple cores.
    
    Args:
        graph: igraph Graph object
        vertex_to_cluster_map: {vertex_id: (component_id, cluster_id)}
        component_cluster_info: {component_id: {'vertices': [...], 'clusters': {...}}}
        vertex_attrs: {vertex_id: {'filename': str, 'feature_idx': int, ...}}
        filter_method: 'simplified-modularity' (weighted internal/external edges) or 'modularity' (proper modularity)
        num_cores: Number of cores to use for parallel processing (default: use all available)
        mega_complexity_threshold: Threshold for using heuristic approach (default: 1e15). Set to 0 or negative to disable.
    
    Returns:
        Tuple of (filtered_vertex_map, filtered_cluster_info, deleted_vertices_list)
    """
    # Copy function parameters to module-level globals for worker process access via fork
    import sys
    current_module = sys.modules[__name__]
    current_module.vertex_attrs = vertex_attrs
    current_module.vertex_to_cluster_map = vertex_to_cluster_map
    current_module.filter_method = filter_method
    current_module.component_cluster_info = component_cluster_info
    
    print("\n" + "="*70)
    print("FILTERING DUPLICATE FILE VERTICES (Parallel Processing)")
    print("="*70)
    
    def _get_actual_available_cpus():
        """
        Get the ACTUAL number of CPUs available to this process,
        respecting cgroup limits (Docker/Nextflow container restrictions).
        
        Returns reported CPUs, cgroup-limited CPUs, and effective CPUs.
        """
        import os
        
        # What the host reports
        reported_cpus = os.cpu_count() or 1
        
        # What cgroups limit us to (if we're in a container)
        cgroup_cpus = None
        try:
            # Try cgroups v2 first (modern Docker)
            with open('/sys/fs/cgroup/cpuset.cpus', 'r') as f:
                cpuset_str = f.read().strip()
                if cpuset_str:
                    # Parse cpuset format (e.g., "0-3" or "0,2,4")
                    import re
                    parts = cpuset_str.split(',')
                    cpu_ids = set()
                    for part in parts:
                        if '-' in part:
                            start, end = part.split('-')
                            cpu_ids.update(range(int(start), int(end) + 1))
                        else:
                            cpu_ids.add(int(part))
                    cgroup_cpus = len(cpu_ids)
        except:
            pass
        
        if cgroup_cpus is None:
            try:
                # Try cgroups v1 (older Docker)
                with open('/sys/fs/cgroup/cpuset/cpuset.cpus', 'r') as f:
                    cpuset_str = f.read().strip()
                    if cpuset_str:
                        import re
                        parts = cpuset_str.split(',')
                        cpu_ids = set()
                        for part in parts:
                            if '-' in part:
                                start, end = part.split('-')
                                cpu_ids.update(range(int(start), int(end) + 1))
                            else:
                                cpu_ids.add(int(part))
                        cgroup_cpus = len(cpu_ids)
            except:
                pass
        
        # What sched_getaffinity reports (what we can theoretically use)
        affinity_cpus = len(os.sched_getaffinity(0))
        
        # Effective CPUs = minimum of what's actually available
        effective_cpus = min(affinity_cpus, cgroup_cpus) if cgroup_cpus else affinity_cpus
        
        return {
            'reported': reported_cpus,
            'affinity': affinity_cpus,
            'cgroup_limited': cgroup_cpus,
            'effective': effective_cpus
        }


    # Determine number of threads to use
    if num_cores is None:
        cpu_info = _get_actual_available_cpus()
        num_cores = cpu_info['effective']
        
        if cpu_info['cgroup_limited'] and cpu_info['cgroup_limited'] < cpu_info['reported']:
            print(f"[CPU INFO] Container is cgroup-limited to {cpu_info['cgroup_limited']} CPUs (host has {cpu_info['reported']})")
            print(f"[CPU INFO] Using {num_cores} cores for parallel processing")
        else:
            print(f"[CPU INFO] Using {num_cores} cores (reported by os.cpu_count())")
    else:
        num_cores = min(num_cores, os.cpu_count() or 1)
    
    # Cap at 10 cores to prevent container resource exhaustion
    num_cores = min(num_cores, 10)
    
    print(f"Using {num_cores} cores for parallel processing\n")
    
    deleted_vertices_list = []
    filtered_vertex_map = dict(vertex_to_cluster_map)
    filtered_cluster_info = {}
    
    total_components = len(component_cluster_info)
    components_with_dups = 0
    components_to_process = []
    
    # Build edge index for fast lookup (maps vertex -> neighbors)
    print(f"PHASE 0: Building edge index for {len(graph.es)} edges...", flush=True)
    import sys
    sys.stdout.flush()
    
    local_edge_index = {}
    for edge in graph.es:
        src, tgt = edge.source, edge.target
        if src not in local_edge_index:
            local_edge_index[src] = []
        if tgt not in local_edge_index:
            local_edge_index[tgt] = []
        local_edge_index[src].append(tgt)
        local_edge_index[tgt].append(src)
    
    # Assign to module globals for worker access via fork
    current_module.edge_index = local_edge_index
    current_module.graph = graph  # Make graph available to workers
    current_module.component_cluster_info = component_cluster_info  # Make component_cluster_info available to workers
    
    print(f"PHASE 0 COMPLETE: Edge index built with {len(local_edge_index)} vertices\n", flush=True)
    sys.stdout.flush()
    
    # First pass: identify components with duplicates and estimate complexity WITHOUT storing full work items
    # This prevents RAM exhaustion for large datasets
    print(f"PHASE 1: Identifying components with duplicate files ({total_components} total)...", flush=True)
    sys.stdout.flush()
    
    components_to_estimate = []  # List of (comp_idx, has_dups) - NOT full work items
    
    for pass_idx, comp_idx in enumerate(sorted(component_cluster_info.keys())):
        if pass_idx % 5000 == 0:
            print(f"  → Analyzed {pass_idx}/{total_components} components...", flush=True)
            sys.stdout.flush()
        
        info = component_cluster_info[comp_idx]
        comp_clusters = info['clusters']
        comp_vertices = info['vertices']
        
        # Skip single-vertex components
        if len(comp_vertices) <= 1:
            filtered_cluster_info[comp_idx] = info
            continue
        
        # Check for duplicates WITHOUT building full work item
        clusters_with_duplicates = {}
        for cluster_id, vertex_ids in comp_clusters.items():
            file_vertex_map = {}
            for vertex_id in vertex_ids:
                filename = vertex_attrs[vertex_id]['filename']
                if filename not in file_vertex_map:
                    file_vertex_map[filename] = []
                file_vertex_map[filename].append(vertex_id)
            
            has_duplicates = any(len(vids) > 1 for vids in file_vertex_map.values())
            if has_duplicates:
                clusters_with_duplicates[cluster_id] = file_vertex_map
        
        if not clusters_with_duplicates:
            # No duplicates in this component
            filtered_cluster_info[comp_idx] = info
            continue
        
        components_with_dups += 1
        
        # Store ONLY the component index and clusters with duplicates for estimation
        # DON'T store full work item yet - we'll build it just-in-time during batching
        components_to_estimate.append((comp_idx, clusters_with_duplicates, comp_vertices))
    
    # ====== PHASE 1.5: Prepare metadata for lazy-loading workers (TRUE lazy-loading) ======
    # Store ONLY lightweight metadata: comp_idx, combo_estimate, comp_vertices, clusters_with_dups, info
    # Full work_items are built on-demand by workers (just-in-time) as they self-batch
    # This keeps orchestrator memory usage minimal - NO full work_items held in memory
    
    if components_to_estimate:
        print(f"PHASE 1.5: Preparing metadata for {len(components_to_estimate)} components with duplicates...", flush=True)
        print(f"  (Full work_items will be built on-demand by workers - TRUE lazy-loading)", flush=True)
        print(f"  Workers will self-batch up to {COMBINATIONS_BATCH_LIMIT:,} combinations per batch", flush=True)
        
        # Prepare lightweight metadata only (not full work_items)
        components_metadata = []
        
        for item_tuple in components_to_estimate:
            comp_idx, clusters_with_dups, comp_vertices = item_tuple
            info = component_cluster_info[comp_idx]
            
            # Estimate complexity for this component (used in worker for batching decision)
            actual_combo_estimate = 1.0
            for cluster_id in clusters_with_dups:
                file_vertex_map = clusters_with_dups[cluster_id]
                num_dup_options = sum(1 for vids in file_vertex_map.values() if len(vids) > 1)
                if num_dup_options > 0:
                    avg_option_size = sum(len(vids) for vids in file_vertex_map.values() if len(vids) > 1) / max(1, num_dup_options)
                    actual_combo_estimate *= (float(max(2, int(avg_option_size))) ** num_dup_options)
                    if actual_combo_estimate > 1e15:
                        actual_combo_estimate = 1e15
                        break
            
            # If complexity >= mega_complexity_threshold, this component will use heuristic (not exhaustive)
            # Heuristic work = O(n²) greedy selection, NOT exponential combinations
            # Use realistic complexity: proportional to vertices involved, not combination count
            if actual_combo_estimate >= mega_complexity_threshold:
                # Heuristic complexity: approximate as O(vertices * clusters * avg_duplicates_per_cluster)
                # This represents the actual greedy iteration work, not the exponential search space
                num_clusters = len(clusters_with_dups)
                num_vertices = len(comp_vertices)
                # Estimate: roughly O(n * c) where n=vertices, c=clusters
                # In practice, heuristic is much faster than exhaustive would suggest
                heuristic_complexity = float(num_vertices * num_clusters) * 10  # Small coefficient for iteration
                metadata_complexity = min(heuristic_complexity, 1000)  # Cap at 1000 for small-ish units
                will_use_heuristic = True
            else:
                # Exhaustive: use actual combo estimate
                metadata_complexity = actual_combo_estimate
                will_use_heuristic = False
            
            # Store ONLY lightweight metadata (no clusters_with_dups!)
            # clusters_with_dups is NOT needed by workers - they rebuild from global data on-demand
            # This saves GIGABYTES when sending metadata to 7 workers in parallel
            metadata = {
                'comp_idx': comp_idx,
                'combo_estimate': metadata_complexity,  # Realistic complexity for load balancing
                'actual_combo_estimate': actual_combo_estimate,  # For logging
                'will_use_heuristic': will_use_heuristic,  # Hint for load balancer
                'comp_vertices': comp_vertices,
                'info': info,
                'mega_complexity_threshold': mega_complexity_threshold
                # NOTE: clusters_with_dups omitted - workers don't use it
            }
            components_metadata.append(metadata)
        
        print(f"PHASE 1.5 COMPLETE: Prepared metadata for {len(components_metadata)} components (lightweight)", flush=True)
        
        # Analyze and log component complexity distribution
        # Separate heuristic vs exhaustive components
        heuristic_components = [m for m in components_metadata if m.get('will_use_heuristic', False)]
        exhaustive_components = [m for m in components_metadata if not m.get('will_use_heuristic', False)]
        
        complexities = [m['combo_estimate'] for m in components_metadata]
        if complexities:
            min_combo = min(complexities)
            max_combo = max(complexities)
            avg_combo = sum(complexities) / len(complexities)
            
            # Count mega-components and large components (after heuristic complexity adjustment)
            mega_components = [m for m in components_metadata if m['combo_estimate'] > COMBINATIONS_BATCH_LIMIT * 10]
            large_components = [m for m in components_metadata if m['combo_estimate'] > COMBINATIONS_BATCH_LIMIT]
            
            print(f"\n  Component Complexity Statistics:", flush=True)
            print(f"    Min:  {min_combo:.2e}", flush=True)
            print(f"    Max:  {max_combo:.2e}", flush=True)
            print(f"    Avg:  {avg_combo:.2e}", flush=True)
            print(f"    ⚠ Mega-components (>10M): {len(mega_components)}", flush=True)
            if mega_components:
                for m in sorted(mega_components, key=lambda x: x['combo_estimate'], reverse=True)[:5]:
                    print(f"      - comp_idx={m['comp_idx']}: {m['combo_estimate']:.2e} combinations, "
                          f"{len(m['comp_vertices'])} vertices", flush=True)
            print(f"    ⚠ Large components (>1M): {len(large_components)}", flush=True)
            
            # Report heuristic vs exhaustive breakdown
            print(f"\n  Heuristic vs Exhaustive Breakdown:", flush=True)
            print(f"    ✓ Heuristic components: {len(heuristic_components)}", flush=True)
            if heuristic_components:
                actual_complexities = [m.get('actual_combo_estimate', 0) for m in heuristic_components]
                metadata_complexities = [m['combo_estimate'] for m in heuristic_components]
                total_actual = sum(actual_complexities)
                total_metadata = sum(metadata_complexities)
                reduction_factor = total_actual / total_metadata if total_metadata > 0 else 1.0
                print(f"      Total actual combinations (before heuristic): {total_actual:.2e}", flush=True)
                print(f"      Total metadata complexity (after adjustment): {total_metadata:.2e}", flush=True)
                print(f"      Reduction factor: {reduction_factor:.1f}x", flush=True)
                print(f"      Top heuristic components:", flush=True)
                for m in sorted(heuristic_components, key=lambda x: x.get('actual_combo_estimate', 0), reverse=True)[:3]:
                    actual = m.get('actual_combo_estimate', 0)
                    metadata = m['combo_estimate']
                    print(f"        - comp_idx={m['comp_idx']}: {actual:.2e} → {metadata:.2e} "
                          f"({actual/metadata:.0f}x reduction)", flush=True)
            print(f"    ✓ Exhaustive components: {len(exhaustive_components)}", flush=True)
        
        sys.stdout.flush()
        
        # Clean up components_to_estimate to free memory
        del components_to_estimate
    else:
        # No components to process
        components_metadata = []
    
    # Initialize worker diagnostics (before processing)
    worker_pids = {}  # pid -> component count
    task_times = []
    cores_to_use = num_cores or (1 if not hasattr(os, 'cpu_count') else (os.cpu_count() or 1))
    
    # Track wall-clock time for multiprocessing
    import time
    multiprocessing_start_time = time.time()
    
    # Process components using multiprocessing with TRUE LAZY LOADING + WORKER-SIDE BATCHING
    # Only metadata is sent to workers; each worker:
    # 1. Receives one metadata dict at a time
    # 2. Builds work_item on-demand
    # 3. Accumulates to 1M combinations
    # 4. Processes batch when full
    # Result: Orchestrator NEVER holds all work_items in memory
    
    if components_metadata:
        print(f"\n{'='*70}")
        print(f"TRUE LAZY LOADING + WORKER-SIDE BATCHING", flush=True)
        print(f"{'='*70}", flush=True)
        print(f"Total metadata items: {len(components_metadata)}", flush=True)
        print(f"Each worker will batch-accumulate up to {COMBINATIONS_BATCH_LIMIT:,} combinations\n", flush=True)
        
        # Initialize results before processing
        results = []
        
        # Check if multiprocessing is enabled
        if not use_multiprocessing:
            print(f"Processing {len(components_metadata)} components using SERIAL mode...\n", flush=True)
            results = _lazy_worker_with_self_batching(iter(components_metadata))
        else:
            # Use multiprocessing with lazy loading
            # Each worker gets its own subset of metadata and builds work_items on-demand
            # OPTIMIZATION: Use 8 workers to maximize throughput with manageable IPC overhead
            # 8 workers on 15 cores leaves headroom for system tasks
            num_workers = min(max(4, 4), cores_to_use)  # Use 8 workers if available, min 2, max available cores
            
            print(f"Multiprocessing optimization: Using {num_workers} workers (reduced from {cores_to_use - 1} to minimize IPC overhead)", flush=True)
            print(f"Each worker processes batches up to {COMBINATIONS_BATCH_LIMIT:,} combinations\n", flush=True)
            
            # Divide metadata into chunks for workers using COMPLEXITY-AWARE LOAD BALANCING
            # Sort by complexity to enable greedy bin-packing
            sorted_metadata = sorted(components_metadata, key=lambda m: m.get('combo_estimate', 0), reverse=True)
            
            # Initialize chunks and track total complexity per chunk
            chunk_complexities = [0] * num_workers
            metadata_chunks = [[] for _ in range(num_workers)]
            
            # Greedy bin-packing: assign each component to the least-loaded chunk
            for metadata in sorted_metadata:
                complexity = metadata.get('combo_estimate', 0)
                min_idx = chunk_complexities.index(min(chunk_complexities))
                metadata_chunks[min_idx].append(metadata)
                chunk_complexities[min_idx] += complexity
            
            # Report load distribution
            print(f"Distributing {len(components_metadata)} metadata items to {num_workers} workers using complexity-aware load balancing", flush=True)
            for i, (chunk, complexity) in enumerate(zip(metadata_chunks, chunk_complexities)):
                print(f"  Worker {i}: {len(chunk)} components, total_complexity={complexity:.2e}", flush=True)
            print(f"  (Load variance: min={min(chunk_complexities):.2e}, max={max(chunk_complexities):.2e})\n", flush=True)
            sys.stdout.flush()
            
            # ===== MULTIPROCESSING TIMING =====
            filter_timing = {}
            filter_timing['pool_start'] = time.time()
            
            try:
                parent_affinity = os.sched_getaffinity(0)
                parent_pid = os.getpid()
                print(f"[MAIN PROCESS PID {parent_pid}] Available CPUs: {parent_affinity} (count: {len(parent_affinity)})\n", flush=True)
                
                print(f"{'='*70}")
                print(f"LAZY LOADING WORKER PROCESSING", flush=True)
                print(f"{'='*70}\n", flush=True)
                
                pool = None
                filter_timing['pool_create_start'] = time.time()
                try:
                    worker_counter = Value('i', 0)
                    counter_lock = Lock()
                    
                    def init_worker_wrapper(counter, lock, num_cpus):
                        """Wrapper that gets worker index and calls _init_worker with only CPU pinning args"""
                        with lock:
                            worker_index = counter.value
                            counter.value += 1
                        _init_worker(worker_index, num_cpus)
                    
                    # Create Pool with ONLY CPU pinning initializer
                    # Global data structures (edge_index, vertex_attrs, etc.) are inherited via fork
                    # This saves ~2GB that would be wasted on copying in initargs!
                    # TODO refactor!
                    # TODO 
                    pool = Pool(processes=num_workers, 
                               initializer=init_worker_wrapper, 
                               initargs=(worker_counter, counter_lock, cores_to_use))
                    filter_timing['pool_create_end'] = time.time()
                    
                    # Process metadata chunks via imap_unordered with STREAMING
                    # Stream results one at a time to avoid holding all in memory
                    # This prevents OOM when collecting results from 8+ workers on 16K+ components
                    chunksize = 1  # One metadata chunk per worker task
                    filter_timing['result_collection_start'] = time.time()
                    print(f"\n{'='*70}", flush=True)
                    print(f"STREAMING RESULTS FROM {num_workers} WORKERS (no batch collection)...", flush=True)
                    print(f"{'='*70}\n", flush=True)
                    
                    results = []  # Initialize early for streaming
                    results_received = 0
                    
                    # Stream results directly without batch collection to reduce memory
                    # TODO Check if that works - maybe overloads ram - try queue?
                    for chunk_results in pool.imap_unordered(
                        _lazy_worker_with_self_batching,
                        [iter(chunk) for chunk in metadata_chunks],
                        chunksize=chunksize
                    ):
                        # Process chunk results as they arrive (streaming)
                        for result in chunk_results:
                            results.append(result)
                            results_received += 1
                            
                            # Progress update every 1000 results
                            if results_received % 1000 == 0:
                                print(f"[MAIN] Streamed {results_received} results from workers...", flush=True)
                    
                    filter_timing['result_collection_end'] = time.time()
                    
                    print(f"\n{'='*70}", flush=True)
                    print(f"RESULT STREAMING COMPLETE: Received {len(results)} total results", flush=True)
                    print(f"{'='*70}\n", flush=True)
                    
                finally:
                    if pool:
                        pool.close()
                        pool.join()
                        filter_timing['pool_end'] = time.time()
                        print(f"Pool closed after processing {len(results)} components\n", flush=True)
                
            except Exception as e:
                import traceback
                print(f"⚠ Multiprocessing failed ({type(e).__name__}: {e})", flush=True)
                print(f"Traceback: {traceback.format_exc()}", flush=True)
                print(f"Falling back to serial processing...", flush=True)
                results = _lazy_worker_with_self_batching(iter(components_metadata))
                filter_timing['pool_end'] = time.time()
        
        # Collect results
        filter_timing['result_processing_start'] = time.time()
        print(f"Processing results (streaming from worker)...", flush=True)
        
        # Free metadata after workers are done (they've all results now)
        del components_metadata
        del sorted_metadata
        del metadata_chunks
        gc.collect()
        print(f"[MAIN] Freed metadata structures to recover memory", flush=True)
        
        # Track per-worker and per-batch timing
        worker_batch_times = {}  # {pid: [(batch_num, time_ms), ...]}
        result_count = 0
        
        for result in results:
            result_count += 1
            worker_pid = result.get('__worker_pid')
            elapsed_ms = result.get('__elapsed_ms', 0)
            batch_num = result.get('__batch_num', 0)
            batch_elapsed_ms = result.get('__batch_elapsed_ms', 0)
            
            if worker_pid:
                worker_pids[worker_pid] = worker_pids.get(worker_pid, 0) + 1
            if elapsed_ms:
                task_times.append(elapsed_ms)
            
            # Track batch times per worker
            if worker_pid and batch_elapsed_ms:
                if worker_pid not in worker_batch_times:
                    worker_batch_times[worker_pid] = []
                worker_batch_times[worker_pid].append((batch_num, batch_elapsed_ms))
            
            if result.get('status') == 'error':
                print(f"⚠ Component {result['comp_idx']} processing error: {result.get('error', 'Unknown')}", flush=True)
                # Use original info on error
                filtered_cluster_info[result['comp_idx']] = result.get('filtered_cluster_info')
                continue
            
            if result.get('status') == 'skipped':
                # No filtering needed for this component
                filtered_cluster_info[result['comp_idx']] = result['filtered_cluster_info']
                continue
            
            comp_idx = result['comp_idx']
            filtered_cluster_info[comp_idx] = result['filtered_cluster_info']
            deleted_vertices_list.extend(result['deleted_vertices'])
            
            # Update filtered_vertex_map
            for deleted_vertex_info in result['deleted_vertices']:
                vertex_id = deleted_vertex_info['vertex_id']
                if vertex_id in filtered_vertex_map:
                    del filtered_vertex_map[vertex_id]
        
        filter_timing['result_processing_end'] = time.time()
        
        # Free accumulated results after processing (we've extracted what we need)
        del results
        gc.collect()
        print(f"[MAIN] Freed results list to recover memory", flush=True)
    
    # Print worker utilization diagnostics
    if worker_pids:
        elapsed_time = time.time() - multiprocessing_start_time
        
        # ===== FILTER TIMING SUMMARY =====
        if filter_timing and 'pool_start' in filter_timing:
            pool_total = filter_timing.get('pool_end', time.time()) - filter_timing['pool_start']
            pool_create = (filter_timing.get('pool_create_end', 0) - filter_timing.get('pool_create_start', 0)) if 'pool_create_start' in filter_timing else 0
            result_collection = (filter_timing.get('result_collection_end', 0) - filter_timing.get('result_collection_start', 0)) if 'result_collection_start' in filter_timing else 0
            flatten = (filter_timing.get('flatten_end', 0) - filter_timing.get('flatten_start', 0)) if 'flatten_start' in filter_timing else 0
            
            # Calculate pure worker processing time (result collection includes all work)
            pure_work = result_collection - (flatten if flatten > 0 else 0)
            
            print(f"\n{'='*70}")
            print(f"FILTERING MULTIPROCESSING BREAKDOWN")
            print(f"{'='*70}")
            print(f"  Pool creation (fork) time......... {pool_create:>8.2f}s")
            print(f"  Result collection (imap_unordered) {result_collection:>8.2f}s")
            print(f"  Result flattening time........... {flatten:>8.2f}s")
            print(f"  Pure worker processing time..... {pure_work:>8.2f}s")
            print(f"  Total multiprocessing........... {pool_total:>8.2f}s")
            print(f"{'='*70}\n", flush=True)
        
        print(f"\n{'='*70}")
        print(f"PROCESS POOL UTILIZATION DIAGNOSTICS")
        print(f"{'='*70}")
        print(f"Unique worker processes: {len(worker_pids)}")
        print(f"Max worker processes available: {cores_to_use}")
        print(f"\nProcess Task Distribution:")
        for pid, count in sorted(worker_pids.items()):
            print(f"  PID {pid}: processed {count} components")
        
        print(f"\n{'='*70}")
        print(f"PER-WORKER BATCH TIMING ANALYSIS")
        print(f"{'='*70}")
        
        # Aggregate timing by worker
        for pid in sorted(worker_batch_times.keys()):
            batch_list = worker_batch_times[pid]
            batch_times_ms = [t for _, t in batch_list]
            total_batch_time_ms = sum(batch_times_ms)
            
            if batch_times_ms:
                avg_batch_time = total_batch_time_ms / len(batch_times_ms)
                max_batch_time = max(batch_times_ms)
                min_batch_time = min(batch_times_ms)
                
                print(f"\nPID {pid}: {len(batch_list)} batches")
                print(f"  Batches: {[b for b, _ in batch_list]}")
                print(f"  Total batch time........ {total_batch_time_ms/1000:>8.2f}s")
                print(f"  Avg batch time......... {avg_batch_time/1000:>8.2f}s")
                print(f"  Min batch time......... {min_batch_time/1000:>8.2f}s")
                print(f"  Max batch time......... {max_batch_time/1000:>8.2f}s")
        
        print(f"{'='*70}\n", flush=True)
        
        # ===== IPC OVERHEAD ANALYSIS =====
        if filter_timing and 'pool_start' in filter_timing and 'result_processing_end' in filter_timing:
            ipc_pool_create = filter_timing.get('pool_create_end', 0) - filter_timing.get('pool_create_start', 0) if 'pool_create_start' in filter_timing else 0
            ipc_collection = filter_timing.get('result_collection_end', 0) - filter_timing.get('result_collection_start', 0) if 'result_collection_start' in filter_timing else 0
            ipc_flatten = filter_timing.get('flatten_end', 0) - filter_timing.get('flatten_start', 0) if 'flatten_start' in filter_timing else 0
            orchestrator_result_processing = filter_timing.get('result_processing_end', 0) - filter_timing.get('result_processing_start', 0) if 'result_processing_start' in filter_timing else 0
            
            total_ipc_overhead = ipc_pool_create + ipc_collection + ipc_flatten + orchestrator_result_processing
            multiproc_total = filter_timing.get('pool_end', time.time()) - filter_timing['pool_start']
            
            print(f"{'='*70}")
            print(f"IPC OVERHEAD ANALYSIS")
            print(f"{'='*70}")
            print(f"  Pool creation (forking)............ {ipc_pool_create:>8.2f}s")
            print(f"  Result collection (streaming)..... {ipc_collection:>8.2f}s")
            print(f"  Result flattening (orchestrator).. {ipc_flatten:>8.2f}s")
            print(f"  Result processing (orchestrator).. {orchestrator_result_processing:>8.2f}s")
            print(f"  TOTAL IPC OVERHEAD............... {total_ipc_overhead:>8.2f}s")
            print(f"  Total multiprocessing (all)....... {multiproc_total:>8.2f}s")
            
            if multiproc_total > 0:
                overhead_pct = (total_ipc_overhead / multiproc_total) * 100
                print(f"  IPC as % of total................ {overhead_pct:>8.1f}%")
            
            print(f"{'='*70}\n")
        
        print(f"\n✅ MULTIPROCESSING: True parallelism, no GIL limitation!")
        print(f"   All {len(worker_pids)} processes run independently with full CPU utilization")
        
        if task_times:
            import statistics
            avg_time = statistics.mean(task_times)
            min_time = min(task_times)
            max_time = max(task_times)
            total_task_time = sum(task_times)
            
            print(f"\nTask timing (ms):")
            print(f"  Average per task: {avg_time:.1f} ms")
            print(f"  Min: {min_time} ms")
            print(f"  Max: {max_time} ms")
            print(f"  Total if serial: {total_task_time:.0f} ms ({total_task_time/1000:.1f}s)")
            print(f"  Actual wall-clock: {elapsed_time*1000:.0f} ms ({elapsed_time:.1f}s)")
            
            # Calculate speedup
            if elapsed_time > 0:
                speedup = total_task_time / (elapsed_time * 1000)
                expected_speedup = len(worker_pids)
                efficiency = (speedup / expected_speedup) * 100 if expected_speedup > 0 else 0
                print(f"\n  Speedup: {speedup:.1f}x (expected: {expected_speedup:.1f}x)")
                print(f"  Efficiency: {efficiency:.1f}% (actual speedup / ideal speedup)")
                
                if efficiency >= 80:
                    print(f"  ✓ GOOD EFFICIENCY: Excellent process-level scaling, no GIL!")
                elif efficiency >= 50:
                    print(f"  ⚠ FAIR EFFICIENCY: Good scaling with per-process parallelism")
                else:
                    print(f"  ⚠ LOW EFFICIENCY: Overhead may still dominate (IPC/pickling)")
            
            # Detect load imbalance
            if max_time > avg_time * 5:
                print(f"\n  ⚠ WARNING: Large time variance detected (max is {max_time/avg_time:.1f}x average)")
        print(f"{'='*70}\n")
    
    print(f"\n{'='*70}")
    print(f"FILTERING SUMMARY")
    print(f"{'='*70}")
    print(f"Components processed: {components_with_dups}")
    print(f"Total vertices deleted: {len(deleted_vertices_list)}")
    print(f"Details saved in output JSON report")
    print(f"{'='*70}\n")
    
    # Sort deleted vertices by component for JSON export
    deleted_vertices_list.sort(key=lambda x: (x['component_id'], x['cluster_id'], x['vertex_id']))
    
    return filtered_vertex_map, filtered_cluster_info, deleted_vertices_list




def _generate_filter_candidates(comp_idx: int, clusters_with_duplicates: Dict, vertex_attrs: Dict, 
                               all_non_duplicate_vertices: set = None,
                               graph: ig.Graph = None, component_cluster_info: Dict = None,
                               mega_complexity_threshold: float = 10000) -> Tuple[List[set], int, bool]:
    """
    Generate all possible combinations of vertices to keep and return as list.
    
    For each duplicate set within a cluster, we must keep exactly one vertex.
    Non-duplicate vertices (alone from their file in a cluster) are always kept.
    All vertices from clusters without duplicates are also always kept.
    
    If complexity >= mega_complexity_threshold, uses heuristic instead: combined degree = internal - external.
    
    Args:
        comp_idx: Component index (for debugging)
        clusters_with_duplicates: {cluster_id: {filename: [vertex_ids]}}
        vertex_attrs: {vertex_id: {...}}
        all_non_duplicate_vertices: Set of all vertices from clusters without duplicates (always keep)
        graph: igraph Graph object (needed for heuristic combined degree calculation)
        component_cluster_info: Component clustering info (needed for heuristic combined degree calculation)
    
    Returns:
        Tuple of (candidates_list, complexity_estimate, use_heuristic_flag)
        - candidates_list: Single-element list with heuristic candidate if complexity >= mega_complexity_threshold, else all candidates
        - complexity_estimate: Estimated number of combinations (accurate, no early capping)
        - use_heuristic_flag: Boolean indicating if heuristic was applied
    """
    if all_non_duplicate_vertices is None:
        all_non_duplicate_vertices = set()
    
    # COMPLEXITY ESTIMATION for heuristic decision
    # Estimate the number of possible combinations (WITHOUT early break - let it run to completion)
    num_duplicate_clusters = len(clusters_with_duplicates)
    combo_estimate = 1
    use_heuristic = False  # Track whether heuristic should be used
    # COMPLEXITY ESTIMATION for heuristic decision
    # Estimate the number of possible combinations (WITHOUT early break - let it run to completion)
    num_duplicate_clusters = len(clusters_with_duplicates)
    combo_estimate = 1
    use_heuristic = False  # Track whether heuristic should be used
    for cluster_id in clusters_with_duplicates:
        file_vertex_map = clusters_with_duplicates[cluster_id]
        num_dup_options = [len(vids) for vids in file_vertex_map.values() if len(vids) > 1]
        for x in num_dup_options:
            combo_estimate *= x
        if combo_estimate >= mega_complexity_threshold:
            use_heuristic = True  # Will use heuristic if true at end of loop
    
    # OPTIMIZATION: If using heuristic, build single heuristic candidate instead of all combinations
    # This saves minutes of candidate cross-multiplication for large components
    if use_heuristic:
        # Use heuristic: for each duplicate set, pick the vertex with highest combined degree
        # Combined degree = internal_connections - external_connections
        # Internal = edges within the cluster, External = edges to other clusters in component
        heuristic_candidate = all_non_duplicate_vertices.copy()
        
        for cluster_id, file_vertex_map in clusters_with_duplicates.items():
            for filename, vertex_ids in file_vertex_map.items():
                if len(vertex_ids) > 1:  # Duplicate set
                    # Calculate combined degree for each candidate vertex
                    best_vertex = None
                    best_combined_degree = None
                    
                    for vertex_id in vertex_ids:
                        if graph is not None and component_cluster_info is not None:
                            # Get cluster vertices and other cluster vertices in this component
                            cluster_vertices = set(component_cluster_info[comp_idx]['clusters'][cluster_id])
                            other_cluster_vertices = set(component_cluster_info[comp_idx]['vertices']) - cluster_vertices
                            
                            # Count internal and external edges
                            internal_degree = 0
                            external_degree = 0
                            
                            if vertex_id < graph.vcount():
                                for neighbor_id in graph.neighbors(vertex_id):
                                    if neighbor_id in cluster_vertices:
                                        internal_degree += 1
                                    elif neighbor_id in other_cluster_vertices:
                                        external_degree += 1
                            
                            combined_degree = internal_degree - external_degree
                        else:
                            # Fallback to simple degree if graph not provided
                            combined_degree = vertex_attrs[vertex_id].get('degree', 0)
                        
                        if best_combined_degree is None or combined_degree > best_combined_degree:
                            best_combined_degree = combined_degree
                            best_vertex = vertex_id
                    
                    if best_vertex is not None:
                        heuristic_candidate.add(best_vertex)
                else:
                    # Single vertex from this file
                    heuristic_candidate.add(vertex_ids[0])
        
        return [heuristic_candidate], combo_estimate, use_heuristic
    
    # Process each cluster independently to maintain cluster structure
    cluster_candidates = {}  # cluster_id -> list of candidate sets for that cluster
    
    for cluster_id in sorted(clusters_with_duplicates.keys()):
        file_vertex_map = clusters_with_duplicates[cluster_id]
        
        # For this cluster, find duplicate sets and non-duplicates
        cluster_choice_options = []  # Options for this cluster's duplicates
        cluster_all_verts = set()    # All vertices in this cluster
        cluster_all_dups = set()     # Vertices from files with >1 vert in this cluster
        
        for filename in sorted(file_vertex_map.keys()):
            vertex_ids = file_vertex_map[filename]
            cluster_all_verts.update(vertex_ids)
            if len(vertex_ids) > 1:  # Duplicate set
                cluster_choice_options.append(vertex_ids)
                cluster_all_dups.update(vertex_ids)
        
        # Vertices to always keep in this cluster (singletons from non-duplicate files)
        cluster_non_dups = cluster_all_verts - cluster_all_dups
        
        # Generate candidates for this cluster
        if not cluster_choice_options:
            # No duplicates in this cluster - keep all vertices
            cluster_candidates[cluster_id] = [cluster_all_verts]
        else:
            # Generate combinations for this cluster only
            cluster_cands = []
            
            def gen_cluster_combos(opt_idx, current):
                if opt_idx == len(cluster_choice_options):
                    # Add non-duplicates to complete this candidate
                    complete = set(current).union(cluster_non_dups)
                    cluster_cands.append(complete)
                    return
                
                for vertex_id in cluster_choice_options[opt_idx]:
                    gen_cluster_combos(opt_idx + 1, current + [vertex_id])
            
            gen_cluster_combos(0, [])
            cluster_candidates[cluster_id] = cluster_cands
    
    # Now cross-multiply candidates from all clusters using GENERATOR (on-demand, not all upfront)
    if not cluster_candidates:
        # Edge case: no clusters with duplicates
        return iter([all_non_duplicate_vertices]), combo_estimate, use_heuristic
    
    # Get list of cluster IDs in order
    cluster_ids = sorted(cluster_candidates.keys())
    
    # Cross-product generator: yield candidates one-at-a-time to avoid holding millions in RAM
    def cross_multiply_generator(cluster_idx, accumulated):
        """Generator that yields final candidate sets on-demand (not all upfront)"""
        if cluster_idx == len(cluster_ids):
            # accumulated is a list of candidate sets, one per cluster
            # Combine them into a single candidate set
            combined = all_non_duplicate_vertices.copy()  # Start with non-duplicate clusters' vertices
            for cluster_cand_set in accumulated:
                combined.update(cluster_cand_set)  # Add vertices from this cluster's candidate
            yield combined
            return
        
        cid = cluster_ids[cluster_idx]
        for candidate_set in cluster_candidates[cid]:
            yield from cross_multiply_generator(cluster_idx + 1, accumulated + [candidate_set])
    
    # Return generator (not list) - candidates produced on-demand during evaluation
    return cross_multiply_generator(0, []), combo_estimate, use_heuristic


def export_filtered_clusters_to_json(vertex_to_cluster_map: Dict, component_cluster_info: Dict, 
                                     vertex_attrs: Dict, output_path: str):
    """Export filtered clusters to JSON file with '_filtered' suffix."""
    # Insert '_filtered' before .json extension
    if output_path.endswith('.json'):
        output_path = output_path.replace('.json', '_filtered.json')
    
    clusters_output = {}
    
    for comp_idx, info in component_cluster_info.items():
        comp_key = f"component_{comp_idx}"
        clusters_output[comp_key] = {
            'component_id': comp_idx,
            'method': info['method'],
            'num_clusters': info['num_clusters'],
            'modularity': info.get('modularity'),
            'details': info['details'],
            'clusters': {}
        }
        
        # For each cluster in this component
        for cluster_id, vertex_ids in info['clusters'].items():
            cluster_key = f"cluster_{cluster_id}"
            cluster_vertices = []
            
            for vertex_id in vertex_ids:
                vertex_info = vertex_attrs[vertex_id]
                cluster_vertices.append({
                    'vertex_id': vertex_id,
                    'filename': vertex_info['filename'],
                    'feature_idx': vertex_info['feature_idx'],
                    'label': vertex_info['label'],
                    'pep_ident': vertex_info.get('pep_ident', []),
                    'prot_ident': vertex_info.get('prot_ident', []),
                    'x_center': vertex_info.get('x_center', 0.0),
                    'y_center': vertex_info.get('y_center', 0.0),
                    'intensity': vertex_info.get('intensity', 0.0),
                    'openms_fid': vertex_info.get('openms_fid', ''),
                    'ms2_scans': vertex_info.get('ms2_scans', ''),
                    'charge': vertex_info.get('charge', 0)
                })
            
            # Aggregate pep_idents and prot_idents for the cluster
            cluster_pep_ident = aggregate_cluster_pep_idents(cluster_vertices)
            cluster_prot_ident = aggregate_cluster_prot_idents(cluster_vertices)
            
            clusters_output[comp_key]['clusters'][cluster_key] = {
                'cluster_id': cluster_id,
                'num_vertices': len(vertex_ids),
                'pep_ident': cluster_pep_ident,
                'prot_ident': cluster_prot_ident,
                'vertices': cluster_vertices
            }
    
    # Save to JSON
    with open(output_path, 'w') as f:
        json.dump(clusters_output, f, indent=2)
    
    print(f"✓ Saved filtered clusters to JSON: {output_path}")
    
    # Print summary
    total_clusters = sum(len(info['clusters']) for info in clusters_output.values())
    total_vertices = sum(len(cluster['vertices']) 
                        for comp in clusters_output.values() 
                        for cluster in comp['clusters'].values())
    print(f"  Summary: {len(clusters_output)} components, {total_clusters} total clusters, {total_vertices} vertices")


def aggregate_cluster_pep_idents(cluster_vertices: List) -> any:
    """
    Aggregate pep_idents from all vertices in a cluster.
    
    Returns:
    - null if no vertices have any pep_idents
    - List of {sequence, percentage, count} for each unique pep_ident found
    """
    pep_ident_counts = {}
    total_with_ident = 0
    
    for vertex in cluster_vertices:
        vertex_pep_idents = vertex.get('pep_ident', [])
        if vertex_pep_idents and len(vertex_pep_idents) > 0:
            for pep_seq in vertex_pep_idents:
                pep_ident_counts[pep_seq] = pep_ident_counts.get(pep_seq, 0) + 1
            total_with_ident += 1
    
    if not pep_ident_counts:
        return None
    
    # Build result with percentages
    total_vertices = len(cluster_vertices)
    result = []
    for pep_seq in sorted(pep_ident_counts.keys()):
        count = pep_ident_counts[pep_seq]
        percentage = (count / total_vertices) * 100
        result.append({
            'sequence': pep_seq,
            'count': count,
            'percentage': round(percentage, 1)
        })
    
    return result


def aggregate_cluster_prot_idents(cluster_vertices: List) -> any:
    """
    Aggregate prot_idents from all vertices in a cluster.
    
    Returns:
    - null if no vertices have any prot_idents
    - List of {id, percentage, count} for each unique prot_ident found
    """
    prot_ident_counts = {}
    total_with_ident = 0
    
    for vertex in cluster_vertices:
        vertex_prot_idents = vertex.get('prot_ident', [])
        if vertex_prot_idents and len(vertex_prot_idents) > 0:
            for prot_id in vertex_prot_idents:
                prot_ident_counts[prot_id] = prot_ident_counts.get(prot_id, 0) + 1
            total_with_ident += 1
    
    if not prot_ident_counts:
        return None
    
    # Build result with percentages
    total_vertices = len(cluster_vertices)
    result = []
    for prot_id in sorted(prot_ident_counts.keys()):
        count = prot_ident_counts[prot_id]
        percentage = (count / total_vertices) * 100
        result.append({
            'id': prot_id,
            'count': count,
            'percentage': round(percentage, 1)
        })
    
    return result


def calculate_file_linkage_scores(recovered_clusters_json: Dict, total_files: int) -> Dict[str, float]:
    """
    Calculate file linkage scores based on cluster composition.
    
    For each file, calculates the average coverage of its features:
    - coverage(feature) = (cluster_size - 1) / (total_files - 1)
    - file_score = average(coverage) across all features in that file
    
    A score of 1.0 indicates all features matched with all other files.
    A score of 0.0 indicates all features are isolated (no matches across files).
    
    Args:
        recovered_clusters_json: Final cluster structure with all vertices
        total_files: Total number of files in the dataset
    
    Returns:
        Dict mapping filename -> linkage_score (0.0 to 1.0)
    """
    # Collect all features per file with their cluster sizes
    file_features = {}  # {filename: [cluster_sizes]}
    
    for component_key, component_data in recovered_clusters_json.items():
        clusters = component_data.get('clusters', {})
        
        for cluster_key, cluster_data in clusters.items():
            cluster_size = cluster_data.get('num_vertices', 0)
            vertices = cluster_data.get('vertices', [])
            
            # Track which files appear in this cluster
            files_in_cluster = set()
            for vertex in vertices:
                filename = vertex.get('filename', 'unknown')
                files_in_cluster.add(filename)
                
                # Add this cluster size to the file's feature list
                if filename not in file_features:
                    file_features[filename] = []
                file_features[filename].append(cluster_size)
    
    # Calculate average coverage per file
    file_linkage_scores = {}
    
    if total_files <= 1:
        # Can't calculate coverage if only 1 file
        for filename in file_features.keys():
            file_linkage_scores[filename] = 0.0
        return file_linkage_scores
    
    for filename, cluster_sizes in file_features.items():
        if not cluster_sizes:
            file_linkage_scores[filename] = 0.0
            continue
        
        # Calculate coverage for each feature
        coverages = []
        for cluster_size in cluster_sizes:
            # Cluster size tells us how many files have this feature
            # If cluster_size=6 and we have 6 files, then feature matches all files
            # If cluster_size=3 and we have 6 files, then feature matches only 3 files
            files_matched = cluster_size - 1  # Don't count the original file
            coverage = files_matched / (total_files - 1)
            coverage = max(0.0, min(1.0, coverage))  # Clamp to [0, 1]
            coverages.append(coverage)
        
        # File score is average coverage of all its features
        file_linkage_scores[filename] = sum(coverages) / len(coverages) if coverages else 0.0
    
    return file_linkage_scores


def save_deleted_vertices_report(deleted_vertices_list: List, output_dir: str = None, filter_method: str = None, num_clusters_with_recoveries: int = 0, total_recovered: int = 0, recovery_tracking: Dict = None):
    """Save deletion report as JSON file."""
    if not deleted_vertices_list:
        print("\nNo vertices were deleted")
        return
    
    # Determine output directory
    if output_dir is None:
        output_dir = '.'  # Save to current directory by default
    elif output_dir.endswith('.json'):
        # If a full path with .json extension is provided, use it directly
        output_path = output_dir
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
        
        # Count recovered vertices
        num_recovered = sum(len(d.get('recovered_vertices', [])) for d in deleted_vertices_list)
        
        # Build vertex ID to filename mapping
        vertex_to_filename = {}
        for deletion in deleted_vertices_list:
            vertex_to_filename[deletion['vertex_id']] = deletion['filename']
        
        # Build recovered vertices list for reporting
        recovered_vertices_list = []
        if recovery_tracking:
            for (comp_idx, cluster_id), vertex_ids in recovery_tracking.items():
                for vertex_id in vertex_ids:
                    filename = vertex_to_filename.get(vertex_id, 'unknown')
                    recovered_vertices_list.append({
                        'vertex_id': vertex_id,
                        'filename': filename,
                        'target_component': comp_idx,
                        'target_cluster': cluster_id
                    })
        
        # Organize by component for readability
        report = {
            'summary': {
                'total_deleted': len(deleted_vertices_list),
                'total_recovered': total_recovered,
                'clusters_with_recovered_vertices': num_clusters_with_recoveries if num_clusters_with_recoveries else 'N/A',
                'timestamp': time.strftime("%Y%m%d_%H%M%S"),
                'deletion_reason': 'Duplicate file vertices removed to optimize clustering modularity',
                'filter_method': filter_method or 'unknown'
            },
            'deleted_by_component': {},
            'recovered_vertices': recovered_vertices_list
        }
        
        # Group by component
        for deletion in deleted_vertices_list:
            comp_id = deletion['component_id']
            if comp_id not in report['deleted_by_component']:
                report['deleted_by_component'][comp_id] = {
                    'component_id': comp_id,
                    'vertices': []
                }
            report['deleted_by_component'][comp_id]['vertices'].append(deletion)
        
        # Save report
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"✓ Saved deletion report: {output_path}")
        return
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Create filename with timestamp
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    output_path = os.path.join(output_dir, f'deleted_vertices_{timestamp}.json')
    
    # Count recovered vertices
    num_recovered = sum(len(d.get('recovered_vertices', [])) for d in deleted_vertices_list)
    
    # Build recovered vertices list for reporting
    recovered_vertices_list = []
    if recovery_tracking:
        for (comp_idx, cluster_id), vertex_ids in recovery_tracking.items():
            for vertex_id in vertex_ids:
                recovered_vertices_list.append({
                    'vertex_id': vertex_id,
                    'target_component': comp_idx,
                    'target_cluster': cluster_id
                })
    
    # Organize by component for readability
    report = {
        'summary': {
            'total_deleted': len(deleted_vertices_list),
            'total_recovered': total_recovered,
            'clusters_with_recovered_vertices': num_clusters_with_recoveries if num_clusters_with_recoveries else 'N/A',
            'timestamp': timestamp,
            'deletion_reason': 'Duplicate file vertices removed to optimize clustering modularity',
            'filter_method': filter_method or 'unknown'
        },
        'deleted_by_component': {},
        'recovered_vertices': recovered_vertices_list
    }
    
    # Group by component
    for deletion in deleted_vertices_list:
        comp_id = deletion['component_id']
        if comp_id not in report['deleted_by_component']:
            report['deleted_by_component'][comp_id] = {
                'component_id': comp_id,
                'vertices': []
            }
        report['deleted_by_component'][comp_id]['vertices'].append(deletion)
    
    # Save report
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"✓ Saved deletion report: {output_path}")


def save_graph(graph: ig.Graph, output_path: str, format: str = 'graphml'):
    """Save graph in various formats."""
    if format == 'graphml':
        graph.write_graphml(output_path)
        print(f"✓ Saved graph (GraphML): {output_path}")
    elif format == 'gml':
        graph.write_gml(output_path)
        print(f"✓ Saved graph (GML): {output_path}")
    elif format == 'edgelist':
        graph.write_edgelist(output_path)
        print(f"✓ Saved graph (EdgeList): {output_path}")
    else:
        raise ValueError(f"Unsupported format: {format}")



def print_analysis(analysis: Dict):
    """Pretty print analysis results."""
    if not analysis:
        print("No graph data to analyze")
        return
    
    print(f"\n{'='*70}")
    print("Network Graph Analysis")
    print(f"{'='*70}")
    
    print(f"\nBasic Statistics:")
    print(f"  Vertices: {analysis.get('num_vertices', 'N/A')}")
    print(f"  Edges: {analysis.get('num_edges', 'N/A')}")
    print(f"  Density: {analysis.get('density', 'N/A'):.4f}")
    print(f"  Average Degree: {analysis.get('avg_degree', 'N/A'):.2f} ± {analysis.get('std_degree', 'N/A'):.2f}")
    print(f"  Max Degree: {analysis.get('max_degree', 'N/A')}")
    print(f"  Min Degree: {analysis.get('min_degree', 'N/A')}")
    
    if 'components' in analysis:
        print(f"\nComponents:")
        print(f"  Number of Components: {analysis['components']['num_components']}")
        print(f"  Largest Component: {analysis['components']['largest_component_size']} vertices")
        if len(analysis['components']['component_sizes']) <= 10:
            print(f"  Component Sizes: {analysis['components']['component_sizes']}")
    
    if 'clustering_coefficient' in analysis:
        print(f"\nClustering Coefficient: {analysis['clustering_coefficient']}")
    
    if 'diameter' in analysis:
        print(f"Diameter: {analysis['diameter']}")
    
    if 'degree_by_file' in analysis:
        print(f"\nDegree Statistics by File:")
        for filename, stats in analysis['degree_by_file'].items():
            print(f"  {Path(filename).stem}:")
            print(f"    Vertices: {stats['num_vertices']}")
            print(f"    Avg Degree: {stats['avg_degree']:.2f}")
            print(f"    Max Degree: {stats['max_degree']}")
            print(f"    Min Degree: {stats['min_degree']}")

def save_deleted_vertices(graph: ig.Graph, deleted_vertices_list: List[Dict], 
                         component_cluster_info: Dict, vertex_attrs: Dict,
                         vertex_to_cluster_map: Dict) -> Tuple[List[Dict], int]:
    """
    Recover deleted vertices by adding them back to adjacent clusters if:
    1. The deleted vertex had edges to that cluster
    2. The cluster doesn't already have a vertex from that file
    
    For each deleted vertex with multiple cluster connections, add to the cluster
    with the strongest connection (most edges).
    
    Args:
        graph: igraph Graph object
        deleted_vertices_list: List of deleted vertex dictionaries
        component_cluster_info: Component clustering information
        vertex_attrs: Vertex attributes dictionary
        vertex_to_cluster_map: Current vertex to cluster mapping
    
    Returns:
        Tuple of (updated_cluster_info, num_clusters_recovered_to)
    """
    if not deleted_vertices_list:
        return component_cluster_info, 0
    
    print(f"\n{'='*70}")
    print("RECOVERING DELETED VERTICES")
    print(f"{'='*70}")
    print(f"Checking {len(deleted_vertices_list)} deleted vertices for recovery candidates...\n")
    
    recovered_cluster_info = {}  # Copy of component_cluster_info to modify
    recovered_clusters_json = {}  # JSON structure for export with full vertex metadata
    cluster_recovery_tracking = {}  # Track recovered vertices per cluster: (comp_idx, cluster_id) -> [vertices]
    
    for comp_idx, info in component_cluster_info.items():
        recovered_cluster_info[comp_idx] = {
            'vertices': list(info['vertices']),
            'clusters': {cid: list(vlist) for cid, vlist in info['clusters'].items()},
            'method': info['method'],
            'modularity': info['modularity'],
            'num_clusters': info['num_clusters'],
            'details': info['details'],
            'recovered_vertices': []  # Track recovered vertices
        }
        
        # Build JSON structure with full vertex metadata
        comp_key = f"component_{comp_idx}"
        recovered_clusters_json[comp_key] = {
            'component_id': comp_idx,
            'method': info['method'],
            'num_clusters': info['num_clusters'],
            'modularity': info.get('modularity'),
            'details': info['details'],
            'clusters': {}
        }
        
        # Add filtered clusters with full vertex info
        for cluster_id, vertex_ids in info['clusters'].items():
            cluster_key = f"cluster_{cluster_id}"
            cluster_vertices = []
            
            for vertex_id in vertex_ids:
                vertex_info = vertex_attrs[vertex_id]
                cluster_vertices.append({
                    'vertex_id': vertex_id,
                    'filename': vertex_info['filename'],
                    'feature_idx': vertex_info['feature_idx'],
                    'label': vertex_info['label'],
                    'pep_ident': vertex_info.get('pep_ident', []),
                    'prot_ident': vertex_info.get('prot_ident', []),
                    'x_center': vertex_info.get('x_center', 0.0),
                    'y_center': vertex_info.get('y_center', 0.0),
                    'intensity': vertex_info.get('intensity', 0.0),
                    'openms_fid': vertex_info.get('openms_fid', ''),
                    'ms2_scans': vertex_info.get('ms2_scans', ''),
                    'charge': vertex_info.get('charge', 0)
                })
            
            # Aggregate pep_idents and prot_idents for the cluster
            cluster_pep_ident = aggregate_cluster_pep_idents(cluster_vertices)
            cluster_prot_ident = aggregate_cluster_prot_idents(cluster_vertices)
            
            recovered_clusters_json[comp_key]['clusters'][cluster_key] = {
                'cluster_id': cluster_id,
                'num_vertices': len(vertex_ids),
                'pep_ident': cluster_pep_ident,
                'prot_ident': cluster_prot_ident,
                'vertices': cluster_vertices
            }
    
    total_recovered = 0
    clusters_with_recoveries = set()
    
    for deleted_vertex_info in deleted_vertices_list:
        deleted_vertex_id = deleted_vertex_info['vertex_id']
        deleted_filename = deleted_vertex_info['filename']
        comp_idx = deleted_vertex_info['component_id']
        orig_cluster_id = deleted_vertex_info['cluster_id']
        
        # Find edges from this deleted vertex to other clusters
        # Get all neighbors of deleted vertex
        neighbors = graph.neighbors(deleted_vertex_id)
        
        # Build map: cluster_id -> weighted_sum (edge_count * edge_weight)
        cluster_weighted_edges = {}
        
        for neighbor_id in neighbors:
            # Check if neighbor is still in the graph (wasn't deleted)
            if neighbor_id not in vertex_to_cluster_map:
                continue
            
            neighbor_comp_idx, neighbor_cluster_id = vertex_to_cluster_map[neighbor_id]
            
            # Only for different clusters
            if neighbor_comp_idx != comp_idx or neighbor_cluster_id == orig_cluster_id:
                continue
            
            # Get edge weight
            edge_id = graph.get_eid(deleted_vertex_id, neighbor_id, directed=False, error=False)
            if edge_id == -1:
                # Edge doesn't exist
                continue
            
            edge_weight = graph.es[edge_id]['weight'] if 'weight' in graph.es[edge_id].attributes() else 1.0
            
            # Sum weighted edges to this cluster
            cluster_key = (neighbor_comp_idx, neighbor_cluster_id)
            cluster_weighted_edges[cluster_key] = cluster_weighted_edges.get(cluster_key, 0.0) + edge_weight
        
        if not cluster_weighted_edges:
            # No edges to other clusters
            continue
        
        # Find the cluster with the strongest weighted connection
        best_cluster_key = max(cluster_weighted_edges, key=cluster_weighted_edges.get)
        target_comp_idx, target_cluster_id = best_cluster_key
        
        # Check if target cluster already has a vertex from this file
        target_cluster_vertices = recovered_cluster_info[target_comp_idx]['clusters'][target_cluster_id]
        file_already_present = any(
            vertex_attrs[v]['filename'] == deleted_filename 
            for v in target_cluster_vertices
        )
        
        if file_already_present:
            # Skip recovery
            continue
        
        # Recover the vertex!
        recovered_cluster_info[target_comp_idx]['clusters'][target_cluster_id].append(deleted_vertex_id)
        recovered_cluster_info[target_comp_idx]['vertices'].append(deleted_vertex_id)
        recovered_cluster_info[target_comp_idx]['recovered_vertices'].append({
            'vertex_id': deleted_vertex_id,
            'filename': deleted_filename,
            'feature_idx': deleted_vertex_info['feature_idx'],
            'label': deleted_vertex_info['label'],
            'target_cluster_id': target_cluster_id,
            'weighted_edge_sum': cluster_weighted_edges[best_cluster_key]
        })
        
        # Also add to JSON structure for recovered clusters export
        comp_key = f"component_{target_comp_idx}"
        cluster_key = f"cluster_{target_cluster_id}"
        recovered_vertex_json = {
            'vertex_id': deleted_vertex_id,
            'filename': deleted_filename,
            'feature_idx': deleted_vertex_info['feature_idx'],
            'label': deleted_vertex_info['label'],
            'pep_ident': deleted_vertex_info.get('pep_ident', []),
            'prot_ident': deleted_vertex_info.get('prot_ident', []),
            'x_center': deleted_vertex_info.get('x_center', 0.0),
            'y_center': deleted_vertex_info.get('y_center', 0.0),
            'intensity': deleted_vertex_info.get('intensity', 0.0),
            'openms_fid': deleted_vertex_info.get('openms_fid', ''),
            'ms2_scans': deleted_vertex_info.get('ms2_scans', ''),
            'charge': deleted_vertex_info.get('charge', 0)
        }
        recovered_clusters_json[comp_key]['clusters'][cluster_key]['vertices'].append(recovered_vertex_json)
        recovered_clusters_json[comp_key]['clusters'][cluster_key]['num_vertices'] += 1
        
        # Recompute cluster pep_ident and prot_ident aggregates after adding the vertex
        cluster_pep_ident = aggregate_cluster_pep_idents(recovered_clusters_json[comp_key]['clusters'][cluster_key]['vertices'])
        recovered_clusters_json[comp_key]['clusters'][cluster_key]['pep_ident'] = cluster_pep_ident
        cluster_prot_ident = aggregate_cluster_prot_idents(recovered_clusters_json[comp_key]['clusters'][cluster_key]['vertices'])
        recovered_clusters_json[comp_key]['clusters'][cluster_key]['prot_ident'] = cluster_prot_ident
        
        # Track recovery per cluster for reporting
        cluster_key = (target_comp_idx, target_cluster_id)
        if cluster_key not in cluster_recovery_tracking:
            cluster_recovery_tracking[cluster_key] = []
        cluster_recovery_tracking[cluster_key].append(deleted_vertex_id)
        
        # Mark vertex as recovered in vertex map
        vertex_to_cluster_map[deleted_vertex_id] = (target_comp_idx, target_cluster_id)
        
        # Mark vertex attributes as recovered
        vertex_attrs[deleted_vertex_id]['recovered'] = True
        vertex_attrs[deleted_vertex_id]['original_cluster'] = orig_cluster_id
        
        total_recovered += 1
        clusters_with_recoveries.add(best_cluster_key)
    
    print(f"✓ Recovered {total_recovered} vertices into {len(clusters_with_recoveries)} cluster(s)")
    for cluster_key in sorted(clusters_with_recoveries):
        comp_idx, target_cluster_id = cluster_key
        recovered_count = len(cluster_recovery_tracking.get(cluster_key, []))
        print(f"  Component {comp_idx}, Cluster {target_cluster_id}: {recovered_count} recovered")
    print()
    
    return recovered_cluster_info, len(clusters_with_recoveries), total_recovered, cluster_recovery_tracking, recovered_clusters_json


def create_singleton_clusters_for_unrecovered_deleted_vertices(
    deleted_vertices_list: List[Dict],
    recovered_vertex_ids: set,
    component_cluster_info: Dict,
    vertex_to_cluster_map: Dict,
    vertex_attrs: Dict
) -> Tuple[Dict, Dict, int]:
    """
    Create singleton clusters for deleted vertices that couldn't be recovered.
    Each unrecovered deleted vertex becomes its own cluster in its own component.
    
    Args:
        deleted_vertices_list: List of all deleted vertex dictionaries
        recovered_vertex_ids: Set of vertex IDs that were successfully recovered
        component_cluster_info: Component clustering structure to update
        vertex_to_cluster_map: Vertex-to-cluster mapping to update
        vertex_attrs: Vertex attributes to update
    
    Returns:
        Tuple of (updated_component_cluster_info, singleton_clusters_json, num_singletons_created)
    """
    if not deleted_vertices_list:
        return component_cluster_info, {}, 0
    
    # Find unrecovered deleted vertices
    unrecovered_deleted = [
        v for v in deleted_vertices_list
        if v['vertex_id'] not in recovered_vertex_ids
    ]
    
    if not unrecovered_deleted:
        print(f"\n✓ All {len(deleted_vertices_list)} deleted vertices were recovered!")
        return component_cluster_info, {}, 0
    
    print(f"\n{'='*70}")
    print("CREATING SINGLETON CLUSTERS FOR UNRECOVERED DELETED VERTICES")
    print(f"{'='*70}")
    print(f"Creating {len(unrecovered_deleted)} singleton cluster(s) for deleted vertices with no recovery options...\n")
    
    # Find the highest existing component and vertex IDs
    next_component_id = max(component_cluster_info.keys()) + 1 if component_cluster_info else 0
    next_vertex_id = max(vertex_to_cluster_map.keys()) + 1 if vertex_to_cluster_map else 0
    
    singleton_clusters_json = {}
    singleton_count = 0
    
    for deleted_vertex_info in unrecovered_deleted:
        # Create new component for this singleton
        component_id = next_component_id
        cluster_id = 0  # Singleton clusters always have cluster_id = 0
        vertex_id = next_vertex_id
        
        # Extract vertex information from deleted vertex record
        filename = deleted_vertex_info['filename']
        feature_idx = deleted_vertex_info['feature_idx']
        label = deleted_vertex_info.get('label', f"deleted_vertex_{vertex_id}")
        pep_ident = deleted_vertex_info.get('pep_ident', [])
        prot_ident = deleted_vertex_info.get('prot_ident', [])
        x_center = deleted_vertex_info.get('x_center', 0.0)
        y_center = deleted_vertex_info.get('y_center', 0.0)
        intensity = deleted_vertex_info.get('intensity', 0.0)
        openms_fid = deleted_vertex_info.get('openms_fid', '')
        ms2_scans = deleted_vertex_info.get('ms2_scans', '')
        charge = deleted_vertex_info.get('charge', 0)
        
        # Build vertex attributes
        vertex_attrs[vertex_id] = {
            'filename': filename,
            'feature_idx': feature_idx,
            'label': label,
            'pep_ident': pep_ident,
            'prot_ident': prot_ident,
            'x_center': x_center,
            'y_center': y_center,
            'intensity': intensity,
            'openms_fid': openms_fid,
            'ms2_scans': ms2_scans,
            'charge': charge,
            'color': '#FF6B6B',  # Red for unrecovered deleted vertices (singleton)
            'recovery_type': 'deleted_singleton',  # Mark as deleted singleton for tracking
            'recovered': False,
            'original_component': deleted_vertex_info['component_id'],
            'original_cluster': deleted_vertex_info['cluster_id'],
            'deletion_reason': deleted_vertex_info.get('reason', 'Duplicate file vertex'),
            'degree': 0  # No edges in singleton
        }
        
        # Add to vertex-to-cluster map
        vertex_to_cluster_map[vertex_id] = (component_id, cluster_id)
        
        # Create component and cluster in component_cluster_info
        component_cluster_info[component_id] = {
            'vertices': [vertex_id],
            'clusters': {cluster_id: [vertex_id]},
            'method': 'singleton',
            'modularity': None,
            'num_clusters': 1,
            'details': f"Singleton cluster for deleted vertex {filename}_{feature_idx}",
            'is_singleton': True,
            'recovery_type': 'deleted_singleton',
            'original_component': deleted_vertex_info['component_id'],
            'original_cluster': deleted_vertex_info['cluster_id'],
            'deletion_reason': deleted_vertex_info.get('reason', 'Duplicate file vertex')
        }
        
        # Add to JSON export structure
        comp_key = f"component_{component_id}"
        cluster_key = f"cluster_{cluster_id}"
        
        singleton_clusters_json[comp_key] = {
            'component_id': component_id,
            'method': 'singleton',
            'num_clusters': 1,
            'modularity': None,
            'details': f"Singleton cluster for deleted vertex {filename}_{feature_idx}",
            'is_singleton': True,
            'recovery_type': 'deleted_singleton',
            'original_component': deleted_vertex_info['component_id'],
            'original_cluster': deleted_vertex_info['cluster_id'],
            'deletion_reason': deleted_vertex_info.get('reason', 'Duplicate file vertex'),
            'clusters': {
                cluster_key: {
                    'cluster_id': cluster_id,
                    'num_vertices': 1,
                    'pep_ident': pep_ident if pep_ident else None,
                    'prot_ident': prot_ident if prot_ident else None,
                    'vertices': [{
                        'vertex_id': vertex_id,
                        'filename': filename,
                        'feature_idx': feature_idx,
                        'label': label,
                        'pep_ident': pep_ident,
                        'prot_ident': prot_ident,
                        'x_center': x_center,
                        'y_center': y_center,
                        'intensity': intensity,
                        'openms_fid': openms_fid,
                        'ms2_scans': ms2_scans,
                        'charge': charge
                    }]
                }
            }
        }
        
        next_vertex_id += 1
        next_component_id += 1
        singleton_count += 1
    
    print(f"✓ Created {singleton_count} singleton cluster(s) for unrecovered deleted vertices:")
    for deleted_vertex_info in unrecovered_deleted:
        print(f"  - {deleted_vertex_info['filename']}_{deleted_vertex_info['feature_idx']} (from comp {deleted_vertex_info['component_id']}, cluster {deleted_vertex_info['cluster_id']})")
    print()
    
    return component_cluster_info, singleton_clusters_json, singleton_count


def add_cutoff_suffix_to_filename(filename: str, edge_cutoff: float = float('inf'), 
                                   mz_cutoff: float = None, rt_cutoff: float = None,
                                   fraction: float = None) -> str:
    """
    Add cutoff and fraction information to output filename before the extension.
    
    Args:
        filename: Original output filename
        edge_cutoff: Edge cutoff value used
        mz_cutoff: m/z cutoff value used
        rt_cutoff: RT cutoff value used
        fraction: Fraction of edges used (for sampling)
    
    Returns:
        Modified filename with cutoff/fraction info inserted before extension
    """
    if filename is None:
        return None
    
    # Start with fraction if specified
    suffix_parts = []
    if fraction is not None:
        suffix_parts.append(f"fraction{fraction}")
    
    # Determine which cutoff was used
    if mz_cutoff is not None and rt_cutoff is not None:
        # Coordinate-based cutoff
        suffix_parts.append(f"mz{mz_cutoff}_rt{rt_cutoff}")
    elif edge_cutoff != float('inf'):
        # Distance-based cutoff
        suffix_parts.append(f"dist{edge_cutoff}")
    
    # Combine parts with underscores
    if suffix_parts:
        cutoff_str = "_" + "_".join(suffix_parts)
    else:
        cutoff_str = ""
    
    if cutoff_str:
        # Insert cutoff info before file extension
        path = Path(filename)
        return str(path.parent / (path.stem + cutoff_str + path.suffix))
    
    return filename


def load_unpaired_vertices_json(unpaired_json_path: str) -> List[Dict]:
    """
    Load unpaired vertices from JSON file generated by pair_features.py.
    
    Args:
        unpaired_json_path: Path to unpaired_vertices JSON file
    
    Returns:
        List of unpaired vertex dictionaries
    """
    if not unpaired_json_path or not os.path.exists(unpaired_json_path):
        return []
    
    try:
        with open(unpaired_json_path, 'r') as f:
            data = json.load(f)
        unpaired_list = data.get('unpaired_vertices', [])
        print(f"\n[PROGRESS] Loaded {len(unpaired_list)} unpaired vertices from: {os.path.basename(unpaired_json_path)}")
        return unpaired_list
    except Exception as e:
        print(f"✗ Warning: Could not load unpaired vertices from {unpaired_json_path}: {e}")
        return []


def create_singleton_clusters_for_unpaired(
    unpaired_vertices: List[Dict],
    component_cluster_info: Dict,
    vertex_to_cluster_map: Dict,
    vertex_attrs: Dict,
    next_vertex_id: int = 0,
    next_component_id: int = 0
) -> Tuple[int, int]:
    """
    Create singleton clusters for unpaired vertices.
    Each unpaired vertex becomes its own cluster in its own component.
    
    Args:
        unpaired_vertices: List of unpaired vertex dictionaries from JSON
        component_cluster_info: Internal clustering structure to update
        vertex_to_cluster_map: Vertex-to-cluster mapping to update
        vertex_attrs: Vertex attributes to update
        next_vertex_id: Starting vertex ID (default: max existing + 1)
        next_component_id: Starting component ID (default: max existing + 1)
    
    Returns:
        Tuple of (next_available_vertex_id, next_available_component_id)
    """
    if not unpaired_vertices:
        return next_vertex_id, next_component_id
    
    print(f"\n[PROGRESS] Creating {len(unpaired_vertices)} singleton clusters for unpaired vertices...")
    
    # Find the highest existing IDs if not provided
    if next_vertex_id == 0 and vertex_to_cluster_map:
        next_vertex_id = max(vertex_to_cluster_map.keys()) + 1
    if next_component_id == 0 and component_cluster_info:
        next_component_id = max(component_cluster_info.keys()) + 1
    
    singleton_count = 0
    for unpaired in unpaired_vertices:
        vertex_id = next_vertex_id
        component_id = next_component_id
        cluster_id = 0  # Singleton clusters always have cluster_id = 0
        
        # Build vertex attributes
        vertex_attrs[vertex_id] = {
            'filename': unpaired.get('filename', 'unknown'),
            'feature_idx': unpaired.get('feature_idx', -1),
            'label': unpaired.get('vertex_id', f'unpaired_{vertex_id}'),
            'pep_ident': unpaired.get('pep_ident', []),
            'prot_ident': unpaired.get('prot_ident', []),
            'x_center': unpaired.get('x_center', 0.0),
            'y_center': unpaired.get('y_center', 0.0),
            'intensity': unpaired.get('intensity', 0.0),
            'openms_fid': unpaired.get('openms_fid', ''),
            'ms2_scans': unpaired.get('ms2_scans', ''),
            'charge': unpaired.get('charge', 0),
            'color': '#808080',  # Gray for unpaired
            'recovery_type': 'unpaired'  # Mark as unpaired for tracking
        }
        
        # Add to vertex-to-cluster map
        vertex_to_cluster_map[vertex_id] = (component_id, cluster_id)
        
        # Create component and cluster in component_cluster_info
        component_cluster_info[component_id] = {
            'vertices': [vertex_id],
            'clusters': {cluster_id: [vertex_id]},
            'method': 'singleton',
            'modularity': None,
            'num_clusters': 1,
            'details': f"Singleton cluster for unpaired vertex {unpaired.get('filename')}_{unpaired.get('feature_idx')}",
            'is_singleton': True,
            'recovery_type': 'unpaired'
        }
        
        next_vertex_id += 1
        next_component_id += 1
        singleton_count += 1
    
    print(f"  ✓ Created {singleton_count} singleton clusters for unpaired vertices")
    return next_vertex_id, next_component_id


def main():
    parser = argparse.ArgumentParser(
        description="Build and analyze network graph from edge dictionary using igraph. "
                    "Input can be edges dictionary {file_idx: [edges]} or paired features file with network_edges."
    )
    parser.add_argument("--input_pkl", required=True, 
                       help="Input pickle file: edges dictionary {file_idx: [edges]} or paired features file")
    parser.add_argument("--edge_cutoff", type=float, default=float('inf'), 
                       help="Maximum euclidean distance for edge inclusion (default: no cutoff)")
    parser.add_argument("--mz_cutoff", type=float, default=None,
                       help="Maximum m/z (x-axis) coordinate difference for edge filtering (overrides --edge_cutoff if set with --rt_cutoff)")
    parser.add_argument("--rt_cutoff", type=float, default=None,
                       help="Maximum RT (y-axis) coordinate difference for edge filtering (overrides --edge_cutoff if set with --mz_cutoff)")
    parser.add_argument("--output_graphml", help="Output GraphML file path")
    parser.add_argument("--output_gml", help="Output GML file path")
    parser.add_argument("--output_edgelist", help="Output edge list file path")
    parser.add_argument("--output_analysis", help="Output analysis JSON file path")
    parser.add_argument("--output_composition", help="Output graph composition JSON file path")
    parser.add_argument("--skip-analysis", action='store_true', 
                       help="Skip graph analysis computation and statistics (speeds up execution)")
    
    # Clustering arguments
    parser.add_argument("--enable-clustering", action='store_true',
                       help="Enable clustering on connected components (detects sub-groups in graphs)")
    parser.add_argument("--clustering-method", type=str, default='louvain',
                       choices=['louvain', 'leiden', 'walktrap', 'label_propagation', 'edge_betweenness', 'spinglass', 'cluster_optimal', 'cluster_faster_greedy'],
                       help="Clustering algorithm to use (default: louvain)")
    parser.add_argument("--clustering-use-weights", action='store_true', default=True,
                       help="Use edge distances as weights in clustering (default: enabled)")
    parser.add_argument("--clustering-no-weights", action='store_false', dest='clustering_use_weights',
                       help="Disable weighted clustering and use unweighted algorithms. If algorithms supoort weighhts, they often pull them directly from the input data. Activating this option on weighted igraph function can lead to better results")
    parser.add_argument("--clustering-weight-mode", type=str, default='inverse',
                       choices=['inverse', 'distance'],
                       help="Weight mode for clustering: inverse (1/(1+distance)) or distance (default: inverse)")
    parser.add_argument("--resolution-parameter", type=float, default=0.25,
                       help="Resolution parameter for leiden clustering (default: 0.5, lower=more clusters, higher=fewer clusters)")
    parser.add_argument("--objective-function", type=str, default='CPM',
                       choices=['CPM', 'modularity'],
                       help="Objective function for leiden clustering (default: CPM)")
    parser.add_argument("--output-clusters-json", help="Output JSON file with cluster information")
    
    # Resolution optimization
    parser.add_argument("--auto-select-resolution", action='store_true',
                       help="Automatically select optimal resolution parameter using resolution optimization")
    parser.add_argument("--resolution-optimization-method", type=str, default='gradient-descent',
                       choices=['gradient-descent', 'grid-search'],
                       help="Resolution optimization method: 'gradient-descent' (adaptive stepping, default) or 'grid-search' (exhaustive search over range)")
    parser.add_argument("--resolution-initial", type=float, default=0.1,
                       help="Initial resolution parameter for gradient descent search (default: 0.1)")
    parser.add_argument("--resolution-min", type=float, default=0.01,
                       help="Minimum resolution value for grid search (default: 0.01, only used with --resolution-optimization-method grid-search)")
    parser.add_argument("--resolution-max", type=float, default=2.0,
                       help="Maximum resolution value for grid search (default: 2.0, only used with --resolution-optimization-method grid-search)")
    parser.add_argument("--resolution-num-points", type=int, default=20,
                       help="Number of evenly spaced resolution points to evaluate in grid search (default: 20, only used with --resolution-optimization-method grid-search)")
    parser.add_argument("--resolution-metric", type=str, default='combined',
                       choices=['combined', 'modularity', 'avg_cluster_size'],
                       help="Optimization metric: 'combined' (maximize cluster size, minimize deletions), 'modularity' (maximize modularity directly), or 'avg_cluster_size' (maximize overall average cluster size of all clusters after filtering & recovery) (default: combined)")
    parser.add_argument("--resolution-lambda", type=float, default=0.1,
                       help="Lambda parameter controlling weight of vertices_deleted penalty in optimization score (default: 0.1, only used with 'combined' metric)")
    parser.add_argument("--resolution-max-iterations", type=int, default=15,
                       help="Maximum number of iterations for gradient descent optimization (default: 15, only used with --resolution-optimization-method gradient-descent)")
    parser.add_argument("--resolution-tolerance", type=float, default=1e-5,
                       help="Convergence tolerance for resolution parameter in scipy Brent's method (default: 1e-5, lower=more precise but slower)")
    parser.add_argument("--output-resolution-iterations", type=str, default=None,
                       help="Output JSON file with resolution optimization iteration history (default: resolution_iterations_TIMESTAMP.json in current directory)")
    
    # Filtering
    parser.add_argument("--filter-duplicate-file-vertices", action='store_true',
                       help="Filter duplicate vertices from same file within clusters, keeping the vertex that maximizes component modularity")
    parser.add_argument("--filter-method", type=str, default='simplified-modularity',
                       choices=['simplified-modularity', 'modularity'],
                       help="Method for selecting which duplicate vertices to keep: 'simplified-modularity' (weighted internal/external edges) or 'modularity' (proper modularity calculation)")
    parser.add_argument("--mega_complexity_threshold", type=float, default=5000,
                       help="Complexity threshold for using heuristic approach instead of exhaustive search. For components with combo_estimate > threshold, uses clusterwise modularity optimization (fast, greedy) instead of evaluating all combinations. Set to 0 or negative to disable heuristic. (default: 1e15)")
    parser.add_argument("--output-deleted-vertices", type=str, default=None,
                       help="Output path for deleted vertices JSON report (default: deleted_vertices_TIMESTAMP.json in current directory)")
    parser.add_argument("--save-deleted-vertices", action='store_true',
                       help="Attempt to recover deleted vertices by adding them to adjacent clusters if they have edges to those clusters and the cluster doesn't have another vertex from the same file")
    parser.add_argument("--disable-multiprocessing", action='store_true',
                       help="Disable multiprocessing for filtering duplicate vertices - use serial processing instead (useful for debugging or when container CPU allocation is limited)")
    
    # Random seed for reproducibility
    parser.add_argument("--random-seed", type=int, default=None,
                       help="Random seed for reproducible results (default: None = no seed, non-deterministic)")
    
    # Combined ident lookup for cluster metadata (pep_ident + prot_ident)
    parser.add_argument("--ident-lookup", type=str, default=None,
                       help="JSON file with combined ident lookup: {filename_feat_idx: {pep_ident: [...], prot_ident: [...]}} (from pair_features.py)")
    parser.add_argument("--unpaired-json", type=str, default=None,
                       help="JSON file with unpaired vertices (features with no edges after pairing) from pair_features.py")
    
    args = parser.parse_args()
    
    # Determine number of cores for parallel processing
    num_cores = os.cpu_count() or 1
    
    # Set random seed if provided
    if args.random_seed is not None:
        import random
        random.seed(args.random_seed)
        np.random.seed(args.random_seed)
        # Ensure igraph uses the seeded Python random module
        ig.set_random_number_generator(random)
        print(f"Random seed set to: {args.random_seed}")
    
    print(f"\n{'='*70}")
    print("Building Network Graph")
    print(f"{'='*70}")
    
    # ===== TIMING INSTRUMENTATION =====
    timing_marks = {}
    timing_marks['start'] = time.time()
    
    # Load edges and cutoff parameters in a SINGLE file operation (avoids double-loading)
    print(f"\nLoading edges and parameters from: {args.input_pkl}")
    timing_marks['load_edges_start'] = time.time()
    edges_full, cutoff_params = load_edges_and_cutoff_params_from_paired_features(args.input_pkl)
    timing_marks['load_edges_end'] = time.time()
    if cutoff_params:
        # Use loaded parameters as defaults if not specified on command line
        if args.edge_cutoff == float('inf') and args.mz_cutoff is None and args.rt_cutoff is None:
            # No command-line cutoffs specified, use loaded values
            if cutoff_params.get('mz_cutoff') is not None and cutoff_params.get('rt_cutoff') is not None:
                args.mz_cutoff = cutoff_params['mz_cutoff']
                args.rt_cutoff = cutoff_params['rt_cutoff']
                print(f"  ✓ Using m/z/RT cutoffs from file: mz={args.mz_cutoff}, rt={args.rt_cutoff}")
            elif cutoff_params.get('edges_cutoff') is not None and cutoff_params['edges_cutoff'] != float('inf'):
                args.edge_cutoff = cutoff_params['edges_cutoff']
                print(f"  ✓ Using edge cutoff from file: {args.edge_cutoff}")
        elif args.mz_cutoff is not None or args.rt_cutoff is not None or args.edge_cutoff != float('inf'):
            # Command-line cutoffs specified, print what's being used
            if args.mz_cutoff is not None and args.rt_cutoff is not None:
                print(f"  ⓘ Using command-line m/z/RT cutoffs: mz={args.mz_cutoff}, rt={args.rt_cutoff}")
            elif args.edge_cutoff != float('inf'):
                print(f"  ⓘ Using command-line edge cutoff: {args.edge_cutoff}")
    
    if not edges_full:
        print("✗ No edges found in input file!")
        return
    
    print(f"  Loaded {len(edges_full)} total edges")
    
    # Use all edges (no sampling in refactored version)
    edges = edges_full
    
    # Build graph
    cutoff_info = f"edge_cutoff={args.edge_cutoff}" if args.mz_cutoff is None else f"mz_cutoff={args.mz_cutoff}, rt_cutoff={args.rt_cutoff}"
    print(f"\nBuilding graph with {cutoff_info}...")
    timing_marks['build_graph_start'] = time.time()
    
    # Load combined ident lookup if provided
    ident_lookup = {}
    if args.ident_lookup and os.path.exists(args.ident_lookup):
        try:
            with open(args.ident_lookup, 'r') as f:
                ident_lookup = json.load(f)
            print(f"Loaded ident_lookup from {os.path.basename(args.ident_lookup)}")
            print(f"  Items in ident_lookup: {len(ident_lookup)}")
            # Debug: Show first few entries
            for i, (key, value) in enumerate(list(ident_lookup.items())[:3]):
                pep_ident = value.get('pep_ident', [])
                print(f"    {key}: pep_ident={pep_ident}")
        except Exception as e:
            print(f"Warning: Could not load ident_lookup: {e}")
    else:
        print(f"No ident_lookup file provided or file not found")
    
    graph, vertex_id_map, vertex_attrs = build_network_graph(edges, edge_cutoff=args.edge_cutoff,
                                                              mz_cutoff=args.mz_cutoff, rt_cutoff=args.rt_cutoff,
                                                              ident_lookup=ident_lookup)
    timing_marks['build_graph_end'] = time.time()
    
    # ===== MEMORY OPTIMIZATION: Delete large temporary data structures =====
    # These are no longer needed after graph construction
    print("[MEMORY] Freeing temporary data structures (ident_lookup, edges)...", flush=True)
    try:
        del ident_lookup       # 2-5 GB - data merged into vertex_attrs, no longer needed
        del edges             # Variable reassignment from sampling, safe to delete
        del edges_full        # 5-8 GB - no longer needed after sampling
        gc.collect()          # Force garbage collection
        print("[MEMORY] ✓ Freed temporary structures", flush=True)
    except NameError:
        pass  # Variables might not exist if not loaded (e.g., no ident_lookup file)
    
    if graph is None:
        return
    
    # Analyze graph (if not skipped)
    analysis = None
    graph_composition = None
    if not args.skip_analysis:
        print("\nAnalyzing graph...")
        timing_marks['analyze_start'] = time.time()
        analysis = analyze_graph(graph, vertex_attrs)
        print_analysis(analysis)
        
        # Also compute graph composition statistics
        print("\nAnalyzing graph composition...")
        graph_composition = analyze_graph_composition(graph, vertex_attrs)
        timing_marks['analyze_end'] = time.time()
        
        # Print composition summary
        comp_summary = graph_composition.get('composition_summary', {})
        if comp_summary:
            print(f"  Total connected components: {comp_summary.get('total_components')}")
            print(f"  Total features in components: {comp_summary.get('total_features')}")
            print(f"  Average features per component: {comp_summary.get('avg_features_per_component', 0):.2f}")
            print(f"  Average files per component: {comp_summary.get('avg_files_per_component', 0):.2f}")
            print(f"  Components by size: {comp_summary.get('components_by_size', {})}")
    else:
        print("Skipping graph analysis (--skip-analysis enabled)")
    
    # Perform clustering (if enabled)
    vertex_to_cluster_map = {}
    component_cluster_info = {}
    optimization_result = {}  # Track optimization results (empty if no optimization)
    optimization_complete = False  # Track if optimization was completed
    
    # Load unpaired vertices (needed for optimization scoring if available)
    unpaired_vertices = load_unpaired_vertices_json(args.unpaired_json)
    
    if args.enable_clustering:
        print("\n" + "="*70)
        print("Clustering Analysis")
        print("="*70)
        timing_marks['clustering_start'] = time.time()
        
        # Initialize clustering data (will be populated if optimization is performed)
        vertex_to_cluster_map = {}
        component_cluster_info = {}
        initial_cluster_distribution_final = {}
        
        # Auto-select resolution if enabled
        if args.auto_select_resolution and args.clustering_method.lower() == 'leiden':
            # Prepare output path for iteration history
            iterations_output_path = args.output_resolution_iterations
            if not iterations_output_path:
                timestamp = time.strftime("%Y%m%d_%H%M%S")
                iterations_output_path = f"resolution_iterations_{timestamp}.json"
            
            # ===== PHASE 2: Setup iteration output directory =====
            # Create directory to store per-iteration clustering results
            iteration_output_dir = None
            try:
                iteration_output_dir = os.path.dirname(iterations_output_path) or '.'
                iterations_subdir = os.path.join(iteration_output_dir, 'iterations')
                os.makedirs(iterations_subdir, exist_ok=True)
                iteration_output_dir = iterations_subdir
                print(f"\n✓ Iteration output directory: {iteration_output_dir}")
            except Exception as e:
                print(f"\n✗ Failed to create iteration output directory: {str(e)}")
                raise
            
            # ===== PHASE 2: Prepare pre-computed analysis for all iterations =====
            # (Compute analysis once, use same for all iterations since graph doesn't change)
            pre_computed_analysis = analysis if analysis is not None else None
            pre_computed_composition = graph_composition if graph_composition is not None else None
            
            # Run optimization - choose method based on parameter
            if args.resolution_optimization_method.lower() == 'grid-search':
                optimal_resolution, best_metrics, iteration_history, _, best_iteration_num = grid_search_resolution(
                    graph,
                    vertex_attrs,
                    method=args.clustering_method.lower(),
                    use_weights=args.clustering_use_weights,
                    weight_mode=args.clustering_weight_mode,
                    objective_function=args.objective_function,
                    apply_filtering=args.filter_duplicate_file_vertices,
                    filter_method=args.filter_method,
                    total_features_after_filter=graph.vcount(),
                    min_resolution=args.resolution_min,
                    max_resolution=args.resolution_max,
                    num_points=args.resolution_num_points,
                    lambda_param=args.resolution_lambda,
                    optimization_metric=args.resolution_metric,
                    output_iterations_path=iterations_output_path,
                    unpaired_vertices=unpaired_vertices,
                    initial_seed=args.random_seed if args.random_seed is not None else 42,
                    mega_complexity_threshold=args.mega_complexity_threshold,
                    iteration_output_dir=iteration_output_dir,
                    pre_computed_analysis=pre_computed_analysis,
                    pre_computed_composition=pre_computed_composition
                )
            else:  # gradient-descent (default)
                optimal_resolution, best_metrics, iteration_history, _, best_iteration_num = select_optimal_resolution(
                    graph,
                    vertex_attrs,
                    method=args.clustering_method.lower(),
                    use_weights=args.clustering_use_weights,
                    weight_mode=args.clustering_weight_mode,
                    objective_function=args.objective_function,
                    apply_filtering=args.filter_duplicate_file_vertices,
                    filter_method=args.filter_method,
                    total_features_after_filter=graph.vcount(),  # Approximate: actual value will be from report generation
                    initial_resolution=args.resolution_initial,
                    lambda_param=args.resolution_lambda,
                    max_iterations=args.resolution_max_iterations,
                    resolution_tolerance=args.resolution_tolerance,
                    optimization_metric=args.resolution_metric,
                    output_iterations_path=iterations_output_path,
                    unpaired_vertices=unpaired_vertices,
                    initial_seed=args.random_seed if args.random_seed is not None else 42,
                    mega_complexity_threshold=args.mega_complexity_threshold,
                    iteration_output_dir=iteration_output_dir,
                    pre_computed_analysis=pre_computed_analysis,
                    pre_computed_composition=pre_computed_composition
                )
            
            # Store optimization result for later reporting
            optimization_result = {
                'optimization_method': args.resolution_optimization_method.lower(),
                'initial_resolution': args.resolution_initial,
                'optimal_resolution': optimal_resolution,
                'optimization_metric': args.resolution_metric,
                'best_score': best_metrics['score'],
                'avg_modularity': best_metrics.get('avg_modularity'),
                'avg_cluster_size': best_metrics['avg_cluster_size'],
                'vertices_deleted': best_metrics['vertices_deleted'],
                'num_clusters': best_metrics['num_clusters'],
                'best_iteration_num': best_iteration_num,
                'iterations': iteration_history,
                'iterations_output_path': iterations_output_path
            }
            
            # Use optimized resolution for final clustering
            final_resolution = optimal_resolution
            print(f"\n  ✓ Using optimized resolution: {optimal_resolution:.4f}")
            print(f"  ✓ Best iteration found: iteration_{best_iteration_num}")
            
            # PHASE 5: Copy clustering results from best iteration (all JSON files already saved in Phase 2)
            best_iteration_dir = os.path.join(iteration_output_dir, f"iteration_{best_iteration_num}")
            print(f"  ✓ Copying clustering results from best iteration directory")
            
            # Define the 5 JSON files to copy
            json_files_to_copy = [
                'clusters_fraction1.0.json',
                'clusters_filtered.json',
                'clusters_recovered.json',
                'deleted_vertices.json',
                'graph_analysis_fraction1.0.json',
                'graph_composition_fraction1.0.json'
            ]
            
            # Get output directory (same as where resolution_iterations.json is)
            output_dir = os.path.dirname(iterations_output_path) or '.'
            
            # Copy each file
            for json_file in json_files_to_copy:
                src = os.path.join(best_iteration_dir, json_file)
                dst = os.path.join(output_dir, json_file)
                try:
                    if os.path.exists(src):
                        shutil.copy(src, dst)
                        print(f"    ✓ Copied {json_file}")
                    else:
                        # Optional file (e.g., clusters_recovered.json might not exist in all iterations)
                        print(f"    ⊘ Skipped {json_file} (not found)")
                except Exception as e:
                    print(f"    ✗ ERROR copying {json_file}: {str(e)}")
                    raise
            
            print(f"\n  ✓ Phase 5 complete: All clustering results from best iteration copied to {output_dir}")
            print(f"    (Original iteration files preserved in {iteration_output_dir} for audit trail)")
            
            # Mark that we have completed the optimization - skip the final re-export
            optimization_complete = True
        else:
            # Non-optimization path: Use provided resolution parameter
            final_resolution = args.resolution_parameter
            optimization_complete = False  # Not using optimization results
            if args.auto_select_resolution and args.clustering_method.lower() != 'leiden':
                print("  ⚠ Warning: --auto-select-resolution requires --clustering-method leiden")
                print(f"  Using provided resolution: {args.resolution_parameter}\n")
            
            # Prepare output path for single-iteration history (for consistency with optimization path)
            iterations_output_path = args.output_resolution_iterations
            if not iterations_output_path:
                timestamp = time.strftime("%Y%m%d_%H%M%S")
                iterations_output_path = f"resolution_iterations_{timestamp}.json"
            
            # Perform clustering with provided resolution
            # Compute seed based on final_resolution for determinism
            resolution_seed = int(final_resolution * 10000) % 10000
            if resolution_seed == 0:
                resolution_seed = 1
            vertex_to_cluster_map, component_cluster_info = cluster_graph_independent_components(
                graph,
                vertex_attrs,
                method=args.clustering_method.lower(),
                use_weights=args.clustering_use_weights,
                weight_mode=args.clustering_weight_mode,
                resolution_parameter=final_resolution,
                objective_function=args.objective_function,
                random_seed=resolution_seed
            )
            
            # Calculate initial (pre-filtering) cluster distribution for this final clustering
            initial_cluster_sizes_final = []
            for comp_id, comp_info in component_cluster_info.items():
                for cluster_id, vertex_list in comp_info.get('clusters', {}).items():
                    initial_cluster_sizes_final.append(len(vertex_list))
            initial_cluster_distribution_final = create_cluster_size_distribution(initial_cluster_sizes_final) if initial_cluster_sizes_final else {}
        timing_marks['clustering_end'] = time.time()
    else:
        print("\nClustering disabled (use --enable-clustering to enable)")
        timing_marks['clustering_start'] = None
        timing_marks['clustering_end'] = None
    
    # Add cutoff suffix to output filenames
    output_graphml = add_cutoff_suffix_to_filename(args.output_graphml, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, None)
    output_gml = add_cutoff_suffix_to_filename(args.output_gml, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, None)
    output_edgelist = add_cutoff_suffix_to_filename(args.output_edgelist, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, None)
    output_analysis = add_cutoff_suffix_to_filename(args.output_analysis, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, None)
    output_composition = add_cutoff_suffix_to_filename(args.output_composition, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, None)
    
    # Save graph
    if output_graphml:
        save_graph(graph, output_graphml, format='graphml')
    
    if output_gml:
        save_graph(graph, output_gml, format='gml')
    
    if output_edgelist:
        save_graph(graph, output_edgelist, format='edgelist')
    
    timing_marks['export_end'] = time.time()
    
    # ===== PRINT TIMING SUMMARY =====
    print(f"\n{'='*70}")
    print("TIMING SUMMARY")
    print(f"{'='*70}")
    
    timing_phases = [
        ('Loading edges', 'load_edges_start', 'load_edges_end'),
        ('Building graph', 'build_graph_start', 'build_graph_end'),
        ('Graph analysis', 'analyze_start', 'analyze_end'),
        ('Clustering', 'clustering_start', 'clustering_end'),
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
        elif phase_name == 'Clustering' and (not timing_marks.get('clustering_start') or not timing_marks.get('clustering_end')):
            print(f"  {phase_name:.<40} {'SKIPPED':>8s}")
        elif phase_name == 'Filtering duplicates' and ('filter_start' not in timing_marks or 'filter_end' not in timing_marks):
            print(f"  {phase_name:.<40} {'SKIPPED':>8s}")
        elif phase_name == 'Graph analysis' and ('analyze_start' not in timing_marks or 'analyze_end' not in timing_marks):
            print(f"  {phase_name:.<40} {'SKIPPED':>8s}")
    
    print(f"{'-'*70}")
    print(f"  {'TOTAL':.<40} {total_elapsed:>8.2f}s")
    print(f"{'='*70}\n")
    
    # Find and report longest phase
    if phase_times:
        longest_phase = max(phase_times, key=phase_times.get)
        longest_time = phase_times[longest_phase]
        pct = (longest_time / total_elapsed * 100) if total_elapsed > 0 else 0
        print(f"⏱️  BOTTLENECK: '{longest_phase}' takes {longest_time:.2f}s ({pct:.1f}% of total)\n")
    
    # Export degree distribution histogram (if not skipped)
    timing_marks['export_start'] = time.time()
    if output_analysis and analysis is not None:
        with open(output_analysis, 'w') as f:
            json.dump(analysis, f, indent=2)
        print(f"✓ Saved analysis: {output_analysis}")
    elif output_analysis and analysis is None:
        print("✗ Cannot save analysis (analysis was skipped)")
    
    # Save graph composition
    if output_composition and graph_composition is not None:
        with open(output_composition, 'w') as f:
            json.dump(graph_composition, f, indent=2)
        print(f"✓ Saved graph composition: {output_composition}")
    elif output_composition and graph_composition is None:
        print("✗ Cannot save graph composition (analysis was skipped)")
    
    # Save clusters (if clustering was enabled)
    if args.enable_clustering and vertex_to_cluster_map:
        # PHASE 2: Create singleton clusters for unpaired vertices BEFORE export
        # (unpaired_vertices already loaded earlier for optimization scoring)
        if not optimization_complete:
            next_vid, next_cid = create_singleton_clusters_for_unpaired(
                unpaired_vertices,
                component_cluster_info,
                vertex_to_cluster_map,
                vertex_attrs,
                next_vertex_id=max(vertex_to_cluster_map.keys()) + 1 if vertex_to_cluster_map else 0,
                next_component_id=max(component_cluster_info.keys()) + 1 if component_cluster_info else 0
            )
        
        # NOTE: We do NOT calculate or update final distributions here
        # The iteration data in resolution_iterations.json represents what the optimization algorithm actually computed
        # Unpaired singletons are ALREADY INCLUDED in the optimization JSON files for metric consistency
        # This ensures data consistency: iteration_history records the true optimization state with singletons included
        
        if not optimization_complete:
            print(f"✓ Added {len(unpaired_vertices) if unpaired_vertices else 0} unpaired singletons to cluster structure")
            print(f"  ✓ resolution_iterations.json remains unchanged (represents pure optimization result)")
            print(f"  ✓ Output files (clusters.json) include unpaired singletons as separate components")
        else:
            print(f"✓ Using pre-computed clustering from best optimization iteration")
            print(f"  ✓ All clustering results including unpaired singletons already in copied files")
        
        # PHASE 5 SIMPLIFIED: Only export if NOT using optimization results
        # (If optimization was used, files were already exported and copied in Phase 2+5)
        if not optimization_complete and args.output_clusters_json:
            output_clusters = add_cutoff_suffix_to_filename(args.output_clusters_json, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff)
            export_clusters_to_json(vertex_to_cluster_map, component_cluster_info, vertex_attrs, output_clusters)
            
            # Filter duplicate file vertices if enabled
            num_clusters_with_recoveries = 0
            filtered_cluster_sizes_final = []
            filtered_cluster_distribution_final = {}
            
            if args.filter_duplicate_file_vertices:
                filtered_vertex_map, filtered_cluster_info, deleted_vertices_list = filter_duplicate_file_vertices(
                    graph, vertex_to_cluster_map, component_cluster_info, vertex_attrs, args.filter_method, num_cores, not args.disable_multiprocessing
                )
                
                # Calculate filtered cluster distribution (after filtering, before recovery)
                filtered_cluster_sizes_final = []
                for comp_id, comp_info in filtered_cluster_info.items():
                    for cluster_id, vertex_list in comp_info.get('clusters', {}).items():
                        filtered_cluster_sizes_final.append(len(vertex_list))
                filtered_cluster_distribution_final = create_cluster_size_distribution(filtered_cluster_sizes_final) if filtered_cluster_sizes_final else {}
                
                # Attempt to recover deleted vertices if enabled
                if args.save_deleted_vertices and deleted_vertices_list:
                    filtered_cluster_info, num_clusters_with_recoveries, total_recovered, recovery_tracking, recovered_clusters_json = save_deleted_vertices(
                        graph, deleted_vertices_list, filtered_cluster_info, vertex_attrs, filtered_vertex_map
                    )
                    
                    # Create singleton clusters for unrecovered deleted vertices
                    recovered_vertex_ids = set()
                    for vertex_list in recovery_tracking.values():
                        recovered_vertex_ids.update(vertex_list)
                    
                    # Get original vertex IDs from deleted_vertices_list
                    deleted_vertex_ids = {v['vertex_id'] for v in deleted_vertices_list}
                    
                    # Create singletons for any deleted vertices not recovered
                    filtered_cluster_info, singleton_clusters_json, num_singletons = create_singleton_clusters_for_unrecovered_deleted_vertices(
                        deleted_vertices_list,
                        recovered_vertex_ids,
                        filtered_cluster_info,
                        filtered_vertex_map,
                        vertex_attrs
                    )
                    
                    # Merge singleton clusters into recovered clusters JSON
                    if num_singletons > 0 and singleton_clusters_json:
                        recovered_clusters_json.update(singleton_clusters_json)
                        
                        # Export singleton clusters JSON separately for reference
                        singletons_output_path = output_clusters.replace('.json', '_singletons.json')
                        with open(singletons_output_path, 'w') as f:
                            json.dump(singleton_clusters_json, f, indent=2)
                        print(f"✓ Saved singleton clusters for unrecovered deleted vertices to JSON: {singletons_output_path}")
                else:
                    total_recovered = 0
                    recovery_tracking = {}
                    recovered_clusters_json = {}
                
                # Calculate file linkage scores (if recovered clusters exist)
                if recovered_clusters_json:
                    # Count unique files in the dataset
                    all_files = set()
                    for component_data in recovered_clusters_json.values():
                        for cluster_data in component_data.get('clusters', {}).values():
                            for vertex in cluster_data.get('vertices', []):
                                all_files.add(vertex.get('filename'))
                    
                    num_files = len(all_files)
                    if num_files > 0:
                        file_linkage_scores = calculate_file_linkage_scores(recovered_clusters_json, num_files)
                        
                        # Add scores to the JSON structure at root level
                        recovered_clusters_json['file_linkage_scores'] = file_linkage_scores
                        print(f"\n✓ Calculated file linkage scores for {num_files} files:")
                        for filename in sorted(file_linkage_scores.keys()):
                            score = file_linkage_scores[filename]
                            print(f"  {filename}: {score:.4f}")
                
                # Export filtered clusters
                if deleted_vertices_list:  # Only export if something was deleted
                    export_filtered_clusters_to_json(filtered_vertex_map, filtered_cluster_info, vertex_attrs, output_clusters)
                    
                    # Save recovered clusters JSON (already built by save_deleted_vertices)
                    if args.save_deleted_vertices and recovered_clusters_json:
                        recovered_output_path = output_clusters.replace('.json', '_recovered.json')
                        with open(recovered_output_path, 'w') as f:
                            json.dump(recovered_clusters_json, f, indent=2)
                        print(f"✓ Saved recovered clusters to JSON: {recovered_output_path}")
                        
                        # Calculate recovered cluster distribution from final state
                        recovered_cluster_sizes_final = []
                        for comp_key, comp_data in recovered_clusters_json.items():
                            # Skip metadata like 'file_linkage_scores'
                            if not isinstance(comp_data, dict) or 'clusters' not in comp_data:
                                continue
                            for cluster_key, cluster_data in comp_data.get('clusters', {}).items():
                                num_vertices = cluster_data.get('num_vertices', 0)
                                recovered_cluster_sizes_final.append(num_vertices)
                        
                        recovered_cluster_distribution_final = create_cluster_size_distribution(recovered_cluster_sizes_final) if recovered_cluster_sizes_final else {}
                        
                        # DO NOT update resolution_iterations.json 
                        # The JSON file represents pure optimization results and should not be modified by post-processing
                        print(f"✓ Cluster filtering complete - {len(recovered_cluster_sizes_final)} clusters after recovery")
                    
                    # Save deletion report (use provided path, or default to clusters directory)
                    deletion_report_path = args.output_deleted_vertices or (os.path.dirname(output_clusters) or 'results')
                    save_deleted_vertices_report(deleted_vertices_list, deletion_report_path, args.filter_method, num_clusters_with_recoveries, total_recovered, recovery_tracking)
            else:
                # Filtering disabled - use initial state as final state
                recovered_cluster_sizes_final = initial_cluster_sizes_final if 'initial_cluster_sizes_final' in locals() else []
                recovered_cluster_distribution_final = initial_cluster_distribution_final if 'initial_cluster_distribution_final' in locals() else {}
                filtered_cluster_distribution_final = initial_cluster_distribution_final if 'initial_cluster_distribution_final' in locals() else {}
        else:
            # No clustering output specified - still update resolution_iterations if possible
            recovered_cluster_sizes_final = initial_cluster_sizes_final if 'initial_cluster_sizes_final' in locals() else []
            recovered_cluster_distribution_final = initial_cluster_distribution_final if 'initial_cluster_distribution_final' in locals() else {}
            filtered_cluster_distribution_final = initial_cluster_distribution_final if 'initial_cluster_distribution_final' in locals() else {}
        
        
        # ===== SAVE SINGLE-ITERATION RESOLUTION DATA (NON-OPTIMIZATION PATH) =====
        # When optimization is OFF, create resolution_iterations.json with single entry
        # This allows reports to be generated with normal quality even without optimization
        if not optimization_complete and 'iterations_output_path' in locals():
            # Calculate metrics for this single resolution iteration
            
            # Count deleted vertices
            vertices_deleted = len(deleted_vertices_list) if 'deleted_vertices_list' in locals() and deleted_vertices_list else 0
            
            # Count recovered vertices  
            vertices_recovered = total_recovered if 'total_recovered' in locals() else 0
            
            # Calculate statistics from the recovered cluster size histogram (same as optimization path)
            if 'recovered_cluster_distribution_final' in locals() and recovered_cluster_distribution_final:
                # recovered_cluster_distribution_final is a dict: {cluster_size: count}
                total_singletons = recovered_cluster_distribution_final.get(1, 0)
                total_clusters = sum(recovered_cluster_distribution_final.values())
                total_vertices = sum(size * count for size, count in recovered_cluster_distribution_final.items())
                avg_cluster_size = total_vertices / total_clusters if total_clusters > 0 else 0.0
                num_clusters = total_clusters
            else:
                # Fallback to initial cluster sizes if recovery data not available
                avg_cluster_size = sum(initial_cluster_sizes_final) / len(initial_cluster_sizes_final) if initial_cluster_sizes_final else 0.0
                num_clusters = len(initial_cluster_sizes_final) if initial_cluster_sizes_final else 0
                total_singletons = sum(1 for size in initial_cluster_sizes_final if size == 1)
            
            # Calculate average modularity across components
            modularity_values = []
            for comp_id, comp_info in component_cluster_info.items():
                modularity = comp_info.get('modularity')
                if modularity is not None:
                    modularity_values.append(modularity)
            avg_modularity = sum(modularity_values) / len(modularity_values) if modularity_values else None
            
            # Count pep_ident conflicts and duplicates
            clusters_with_conflicting_pep_idents = 0
            components_with_cluster_pep_ident_duplicates = 0
            for comp_id, comp_info in component_cluster_info.items():
                # Check for inter-cluster pep_ident duplicates
                if comp_info.get('num_clusters', 0) > 1:
                    # Need to check for duplicates across clusters in this component
                    pep_idents_by_cluster = {}
                    for cluster_id, vertex_list in comp_info.get('clusters', {}).items():
                        pep_idents = set()
                        for v_id in vertex_list:
                            pep_ident = vertex_attrs.get(v_id, {}).get('pep_ident')
                            if pep_ident:
                                # Handle both list and single value
                                if isinstance(pep_ident, list):
                                    pep_idents.update(pep_ident)
                                else:
                                    pep_idents.add(pep_ident)
                        pep_idents_by_cluster[cluster_id] = pep_idents
                    
                    # Check for duplicates across clusters
                    all_pep_idents = set()
                    has_duplicate = False
                    for pep_idents in pep_idents_by_cluster.values():
                        if all_pep_idents & pep_idents:  # Intersection
                            has_duplicate = True
                            break
                        all_pep_idents.update(pep_idents)
                    
                    if has_duplicate:
                        components_with_cluster_pep_ident_duplicates += 1
                
                # Check for intra-cluster conflicts (multiple different pep_idents in same cluster)
                for cluster_id, vertex_list in comp_info.get('clusters', {}).items():
                    pep_idents = set()
                    for v_id in vertex_list:
                        pep_ident = vertex_attrs.get(v_id, {}).get('pep_ident')
                        if pep_ident:
                            # Handle both list and single value
                            if isinstance(pep_ident, list):
                                pep_idents.update(pep_ident)
                            else:
                                pep_idents.add(pep_ident)
                    if len(pep_idents) > 1:
                        clusters_with_conflicting_pep_idents += 1
            
            # Calculate score (using same formula as optimization)
            vertices_deleted_ratio = vertices_deleted / graph.vcount() if graph.vcount() > 0 else 0.0
            score = avg_cluster_size - args.resolution_lambda * vertices_deleted_ratio
            
            # Build single-iteration history entry
            iteration_data = {
                'iteration': 1,
                'resolution': final_resolution,
                'vertices_deleted': vertices_deleted,
                'vertices_recovered': vertices_recovered,
                'avg_cluster_size': avg_cluster_size,
                'avg_modularity': avg_modularity,
                'score': score,
                'num_clusters': num_clusters,
                'clusters_with_conflicting_pep_idents': clusters_with_conflicting_pep_idents,
                'components_with_cluster_pep_ident_duplicates': components_with_cluster_pep_ident_duplicates,
                'initial_cluster_distribution': initial_cluster_distribution_final if 'initial_cluster_distribution_final' in locals() else {},
                'filtered_cluster_distribution': filtered_cluster_distribution_final if 'filtered_cluster_distribution_final' in locals() else {},
                'cluster_recovered_final': recovered_cluster_distribution_final if 'recovered_cluster_distribution_final' in locals() else {}
            }
            
            # Create output structure matching optimization format
            iterations_output = {
                'optimization_config': {
                    'method': 'fixed_resolution',
                    'clustering_method': args.clustering_method.lower(),
                    'resolution': final_resolution,
                    'optimization_metric': 'N/A',
                    'lambda': args.resolution_lambda,
                    'filter_enabled': args.filter_duplicate_file_vertices,
                    'filter_method': args.filter_method if args.filter_duplicate_file_vertices else None
                },
                'iterations': [iteration_data]
            }
            
            # Save to file
            try:
                with open(iterations_output_path, 'w') as f:
                    json.dump(iterations_output, f, indent=2)
                print(f"✓ Saved single-iteration resolution data: {iterations_output_path}")
                print(f"  (Fixed resolution={final_resolution:.4f}, avg_cluster_size={avg_cluster_size:.3f}, score={score:.4f})")
            except Exception as e:
                print(f"✗ Failed to save resolution iterations JSON: {str(e)}")
        
        # resolution_iterations.json is NOT modified after optimization
        # It represents the pure optimization result without post-processing
        
        if not args.output_clusters_json:
            print("Clustering performed but --output-clusters-json not specified. Cluster data not saved.")
    
    # Return graph for interactive use
    return graph, vertex_attrs, analysis


if __name__ == "__main__":
    main()

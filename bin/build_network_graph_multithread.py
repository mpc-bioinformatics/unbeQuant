#!/usr/bin/env python3
"""
Build a network graph from paired feature edges using igraph.
Visualize and analyze feature connections across files.
"""

import argparse
import pickle
import json
import numpy as np
import os
import shutil
import sys
import math
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
    import graphviz
    GRAPHVIZ_AVAILABLE = True
except ImportError:
    GRAPHVIZ_AVAILABLE = False

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


def load_edges_from_paired_features(paired_features_file: str) -> List[Dict]:
    """Load edges from paired features pickle or JSON file."""
    if paired_features_file.endswith('.pkl'):
        with open(paired_features_file, 'rb') as f:
            data = pickle.load(f)
    elif paired_features_file.endswith('.json'):
        with open(paired_features_file, 'r') as f:
            data = json.load(f)
    else:
        raise ValueError("File must be .pkl or .json")
    
    # Handle multiple edge formats
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
    
    return edges


def load_cutoff_params_from_paired_features(pkl_path: str) -> Dict:
    """Load cutoff parameters from paired features pickle file.
    Returns dict with keys: mz_cutoff, rt_cutoff, edges_cutoff (or empty dict if not found).
    """
    try:
        if pkl_path.endswith('.pkl'):
            with open(pkl_path, 'rb') as f:
                data = pickle.load(f)
        elif pkl_path.endswith('.json'):
            with open(pkl_path, 'r') as f:
                data = json.load(f)
        else:
            return {}
        
        if isinstance(data, dict) and 'cutoff_params' in data:
            return data['cutoff_params']
        
        return {}
    
    except Exception as e:
        print(f"  ⓘ Note: Could not load cutoff parameters from {pkl_path}: {e}")
        return {}


def extract_vertex_coordinates_from_edges(edges: List[Dict], vertex_id_map: Dict) -> Dict[int, Tuple[float, float]]:
    """
    Extract original (m/z, RT) coordinates for each vertex from edge data.
    
    Args:
        edges: List of edge dictionaries with coordinate information
        vertex_id_map: Mapping of (filename, feature_idx) -> vertex_id
    
    Returns:
        Dictionary mapping vertex_id -> (x_coord, y_coord) from original heatmap space
    """
    vertex_coords = {}
    
    # Extract coordinates from edges
    # Each edge has file1/file2 or current_file/matched_file with feature references
    # We need to find any edge that references each vertex to get its coordinates
    
    # Build a coordinate lookup from edges
    # Edges don't directly contain coordinates, but we can infer relative positions from distances
    # For now, return empty dict - coordinates need to be passed separately
    
    print("  Note: Coordinate extraction from edges requires feature data with x_center, y_center fields")
    return vertex_coords


def load_feature_coordinates_from_jsons(feature_json_files: List[str]) -> Dict[Tuple[str, int], Tuple[float, float]]:
    """
    Load feature coordinates (x_center, y_center) from feature data JSON files.
    
    Args:
        feature_json_files: List of paths to feature data JSON files
    
    Returns:
        Dictionary mapping (filename, feature_idx) -> (x_center, y_center)
    """
    print(f"\nLoading feature coordinates from {len(feature_json_files)} JSON files...")
    
    coord_map = {}
    
    for json_file in feature_json_files:
        try:
            with open(json_file, 'r') as f:
                features = json.load(f)
            
            if not isinstance(features, list):
                print(f"  ⚠ Warning: {json_file} does not contain a feature list")
                continue
            
            # Extract coordinates for each feature
            for feature in features:
                if 'filename' not in feature or '_feat_idx' not in feature:
                    continue
                if 'x_center' not in feature or 'y_center' not in feature:
                    continue
                
                key = (feature['filename'], feature['_feat_idx'])
                coord_map[key] = (float(feature['x_center']), float(feature['y_center']))
            
            print(f"  ✓ Loaded {len(features)} features from {Path(json_file).name}")
        
        except Exception as e:
            print(f"  ✗ Error loading {json_file}: {e}")
            continue
    
    print(f"  ✓ Total features with coordinates: {len(coord_map)}")
    return coord_map


def map_coordinates_to_vertices(coord_map: Dict[Tuple[str, int], Tuple[float, float]], 
                                 vertex_id_map: Dict[Tuple[str, int], int],
                                 vertex_attrs: Dict) -> Dict[int, Tuple[float, float]]:
    """
    Map feature coordinates to vertex IDs and normalize for visualization.
    
    Args:
        coord_map: Dictionary mapping (filename, feature_idx) -> (x, y)
        vertex_id_map: Dictionary mapping (filename, feature_idx) -> vertex_id
        vertex_attrs: Dictionary of vertex attributes
    
    Returns:
        Dictionary mapping vertex_id -> (x_scaled, y_scaled) for Graphviz
    """
    print("\nMapping coordinates to vertices...")
    
    vertex_coords = {}
    
    # Map coordinates to vertex IDs
    for vertex_key, vertex_id in vertex_id_map.items():
        if vertex_key in coord_map:
            vertex_coords[vertex_id] = coord_map[vertex_key]
    
    if not vertex_coords:
        print("  ✗ No coordinate matches found")
        return {}
    
    print(f"  ✓ Matched {len(vertex_coords)}/{len(vertex_id_map)} vertices with coordinates")
    
    # Normalize coordinates for visualization
    # Find coordinate ranges
    x_coords = [x for x, y in vertex_coords.values()]
    y_coords = [y for x, y in vertex_coords.values()]
    
    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)
    
    x_range = max_x - min_x if max_x > min_x else 1.0
    y_range = max_y - min_y if max_y > min_y else 1.0
    
    # Scale to a reasonable coordinate range for Graphviz (e.g., 0-100)
    scale_factor = 100.0
    
    normalized_coords = {}
    for vertex_id, (x, y) in vertex_coords.items():
        # Normalize to 0-1 range, then scale
        x_norm = ((x - min_x) / x_range) * scale_factor
        y_norm = ((y - min_y) / y_range) * scale_factor
        normalized_coords[vertex_id] = (x_norm, y_norm)
    
    print(f"  Coordinate ranges: x=[{min_x:.2f}, {max_x:.2f}], y=[{min_y:.2f}, {max_y:.2f}]")
    print(f"  Scaled to: x=[0, {scale_factor}], y=[0, {scale_factor}]")
    
    return normalized_coords


def get_file_colors(filenames: List[str]) -> Dict[str, str]:
    """
    Generate distinct colors for each unique filename.
    Uses a palette of distinct, visually appealing colors.
    
    Args:
        filenames: List of unique filenames
    
    Returns:
        Dictionary mapping filename -> hex color code
    """
    # Define a palette of distinct, nice colors (hex codes)
    palette = [
        '#FF6B6B',  # Red
        '#4ECDC4',  # Teal
        '#45B7D1',  # Blue
        '#FFA07A',  # Light Salmon
        '#98D8C8',  # Mint
        '#F7DC6F',  # Yellow
        '#BB8FCE',  # Purple
        '#85C1E2',  # Sky Blue
        '#F8B88B',  # Peach
        '#AED6F1',  # Light Blue
        '#F5B041',  # Orange
        '#82E0AA',  # Light Green
        '#F1948A',  # Light Red
        '#85C1E2',  # Steel Blue
        '#D7BDE2',  # Lavender
    ]
    
    # Create mapping with cycling through palette if more files than colors
    colors = {}
    for idx, filename in enumerate(sorted(set(filenames))):
        color_idx = idx % len(palette)
        colors[filename] = palette[color_idx]
    
    return colors


def compute_bidirectional_map(edges: List[Dict]) -> Dict[Tuple, bool]:
    """
    Pre-compute which edges are bidirectional in the FULL dataset.
    Returns a map: (src_key, tgt_key) -> is_bidirectional
    
    This must be done on the complete dataset before sampling,
    since bidirectional pairs are scattered throughout the data.
    
    Args:
        edges: Complete list of all edges
    
    Returns:
        Dictionary mapping (src, tgt) -> True if bidirectional exists
    """
    # Build a set of all directed edges for quick lookup
    edge_set = set()
    for edge in edges:
        if 'file1' in edge:
            src_key = (edge['file1']['filename'], edge['file1']['feature_idx'])
            tgt_key = (edge['file2']['filename'], edge['file2']['feature_idx'])
        else:
            src_key = (edge['current_file']['filename'], edge['current_file']['feature_idx'])
            tgt_key = (edge['matched_file']['filename'], edge['matched_file']['feature_idx'])
        edge_set.add((src_key, tgt_key))
    
    # Now build map of which edges are bidirectional
    bidirectional_map = {}
    for (src, tgt) in edge_set:
        # Edge is bidirectional if the reverse edge also exists
        is_bidir = (tgt, src) in edge_set
        bidirectional_map[(src, tgt)] = is_bidir
    
    bidir_count = sum(1 for v in bidirectional_map.values() if v)
    print(f"  Pre-computed bidirectional edges: {bidir_count}/{len(bidirectional_map)} ({100*bidir_count/len(bidirectional_map):.1f}%)")
    
    return bidirectional_map


def sample_edges_for_testing(edges: List[Dict], fraction: float) -> List[Dict]:
    """
    Sample complete feature groups (connected components) for testing.
    Instead of randomly picking edges, this traces all connections within
    feature groups and includes complete groups until edge budget is full.
    
    Args:
        edges: List of edge dictionaries
        fraction: Fraction of edges to keep (0.0 to 1.0)
    
    Returns:
        Sampled list of edges (complete feature groups)
    """
    if fraction <= 0.0 or fraction >= 1.0:
        if fraction >= 1.0:
            return edges
        else:
            return []
    
    target_edges = max(10, int(len(edges) * fraction))
    
    # Build adjacency map: vertex -> list of edge indices
    vertex_to_edges = {}
    
    for edge_idx, edge in enumerate(edges):
        if 'file1' in edge:
            v1 = (edge['file1']['filename'], edge['file1']['feature_idx'])
            v2 = (edge['file2']['filename'], edge['file2']['feature_idx'])
        else:
            v1 = (edge['current_file']['filename'], edge['current_file']['feature_idx'])
            v2 = (edge['matched_file']['filename'], edge['matched_file']['feature_idx'])
        
        if v1 not in vertex_to_edges:
            vertex_to_edges[v1] = []
        if v2 not in vertex_to_edges:
            vertex_to_edges[v2] = []
        
        vertex_to_edges[v1].append(edge_idx)
        vertex_to_edges[v2].append(edge_idx)
    
    # Find connected components using BFS
    visited_vertices = set()
    components = []  # Each component is a set of edge indices
    
    for start_vertex in vertex_to_edges:
        if start_vertex in visited_vertices:
            continue
        
        # BFS to find all vertices and edges in this component
        component_edges = set()
        queue = [start_vertex]
        visited_vertices.add(start_vertex)
        
        while queue:
            vertex = queue.pop(0)
            
            # Add all edges connected to this vertex
            for edge_idx in vertex_to_edges[vertex]:
                component_edges.add(edge_idx)
                
                # Find the other vertex of this edge
                edge = edges[edge_idx]
                if 'file1' in edge:
                    v1 = (edge['file1']['filename'], edge['file1']['feature_idx'])
                    v2 = (edge['file2']['filename'], edge['file2']['feature_idx'])
                else:
                    v1 = (edge['current_file']['filename'], edge['current_file']['feature_idx'])
                    v2 = (edge['matched_file']['filename'], edge['matched_file']['feature_idx'])
                
                other_vertex = v2 if v1 == vertex else v1
                
                if other_vertex not in visited_vertices:
                    visited_vertices.add(other_vertex)
                    queue.append(other_vertex)
        
        components.append(component_edges)
    
    # Sort components by edge count (largest first) for better representation
    components.sort(key=len, reverse=True)
    
    # Greedily select complete components until edge budget is full
    selected_edges = set()
    selected_components = 0
    
    for component in components:
        if len(selected_edges) + len(component) <= target_edges:
            # Add complete component
            selected_edges.update(component)
            selected_components += 1
        elif len(selected_edges) < target_edges:
            # Last component: only add if it would get us closer to the target
            remaining = target_edges - len(selected_edges)
            component_list = list(component)
            if len(component_list) <= remaining:
                # Component fits completely
                selected_edges.update(component)
                selected_components += 1
            else:
                # Partial fill only at the very end
                np.random.seed(42)
                partial = np.random.choice(len(component_list), size=remaining, replace=False)
                selected_edges.update([component_list[i] for i in partial])
                selected_components += 1
            break
        else:
            break
    
    sampled_edges = [edges[i] for i in sorted(selected_edges)]
    
    # Count unique vertices
    vertices = set()
    for edge in sampled_edges:
        if 'file1' in edge:
            v1 = (edge['file1']['filename'], edge['file1']['feature_idx'])
            v2 = (edge['file2']['filename'], edge['file2']['feature_idx'])
        else:
            v1 = (edge['current_file']['filename'], edge['current_file']['feature_idx'])
            v2 = (edge['matched_file']['filename'], edge['matched_file']['feature_idx'])
        vertices.add(v1)
        vertices.add(v2)
    
    print(f"  Smart sampling: {len(sampled_edges)} edges from {selected_components} complete feature groups ({len(vertices)} unique vertices)")
    
    return sampled_edges


def build_network_graph(edges: List[Dict], edge_cutoff: float = float('inf'), 
                       mz_cutoff: float = None, rt_cutoff: float = None,
                       bidirectional_map: Dict = None, ident_lookup: Dict = None) -> Tuple[ig.Graph, Dict]:
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
        bidirectional_map: Pre-computed map of (src, tgt) -> is_bidirectional from full dataset
    
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
    
    # Debug: Compare generated keys with ident_lookup keys
    if ident_lookup:
        lookup_keys = set(ident_lookup.keys())
        generated_keys = set()
        matched_keys = set()
        unmatched_keys = set()
        
        for idx, attrs in vertex_attrs.items():
            lookup_key = f"{attrs['filename']}_{attrs['feature_idx']}"
            generated_keys.add(lookup_key)
            if lookup_key in lookup_keys:
                matched_keys.add(lookup_key)
            else:
                unmatched_keys.add(lookup_key)
        
        print(f"\n  Lookup key matching:")
        print(f"    Generated keys: {len(generated_keys)}")
        print(f"    Ident_lookup keys: {len(lookup_keys)}")
        print(f"    Matched keys: {len(matched_keys)}")
        print(f"    Unmatched keys: {len(unmatched_keys)}")
        sys.stdout.flush()
        
        # Show sample unmatched keys to debug mismatch
        if unmatched_keys and len(unmatched_keys) <= 5:
            print(f"    Unmatched samples: {list(unmatched_keys)[:3]}")
        if matched_keys:
            sample_matched = list(matched_keys)[0]
            print(f"    Sample matched: {sample_matched}")
            print(f"      pep_ident in lookup: {ident_lookup[sample_matched].get('pep_ident', [])}")
        sys.stdout.flush()
    
    # Generate colors for each source file
    unique_filenames = sorted(set(filename for filename, _ in vertices))
    file_colors = get_file_colors(unique_filenames)
    
    # Assign color to each vertex based on its source file
    for vertex_id, attrs in vertex_attrs.items():
        filename = attrs['filename']
        attrs['color'] = file_colors[filename]
    
    # Print file color mapping
    print("\nFile color mapping:")
    for filename, color in file_colors.items():
        vertex_count = sum(1 for v_filename, _ in vertices if v_filename == filename)
        print(f"  {Path(filename).stem}: {color} ({vertex_count} vertices)")
    
    # Create graph with vertices
    graph = ig.Graph(len(vertices))
    
    # Add vertex attributes
    for vertex_id, attrs in vertex_attrs.items():
        for attr_name, attr_value in attrs.items():
            graph.vs[vertex_id][attr_name] = attr_value
    
    # Add edges with distance weights
    # Edge direction: current_file (source) → matched_file (target)
    edge_list = []
    edge_weights = []
    edge_sources = []  # Track edge source vertex IDs for validation
    edge_targets = []  # Track edge target vertex IDs for validation
    edge_bidirectional = []  # Track if edge is bidirectional
    
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
        edge_sources.append(vertex1_id)
        edge_targets.append(vertex2_id)
        
        # Check if edge is bidirectional (using pre-computed map if available)
        is_bidir = False
        if bidirectional_map:
            is_bidir = bidirectional_map.get((vertex1_key, vertex2_key), False)
        edge_bidirectional.append(is_bidir)
    
    # Add all edges to graph at once
    graph.add_edges(edge_list)
    
    # Add edge distances and track source/target for visualization
    for edge_id, weight in enumerate(edge_weights):
        graph.es[edge_id]['distance'] = weight
        graph.es[edge_id]['source_vertex'] = edge_sources[edge_id]
        graph.es[edge_id]['target_vertex'] = edge_targets[edge_id]
        graph.es[edge_id]['is_bidirectional'] = edge_bidirectional[edge_id]
    
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


def export_degree_distribution_histogram(graph: ig.Graph, vertex_attrs: Dict, output_path: str, graph_composition: Dict = None):
    """Export degree distribution and component size distribution as histogram images."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  ⓘ Note: matplotlib not available for histogram export")
        return False
    
    degrees = graph.degree()
    
    # Create subplot grid: 2 rows, 3 columns
    # Top: full degree dist (all), full degree dist (by file), component size (linear Y-axis)
    # Bottom: zoomed 95% quartile (all), zoomed 95% quartile (by file), component size (log Y-axis)
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    ax1, ax2, ax3_linear = axes[0]  # Top row
    ax1_zoom, ax2_zoom, ax3 = axes[1]  # Bottom row
    
    # Smart binning for degree distribution
    max_degree = max(degrees) if degrees else 1
    mean_degree = np.mean(degrees) if degrees else 0
    std_degree = np.std(degrees) if degrees else 1
    
    # Determine appropriate number of bins
    num_bins = min(50, max(int(max_degree / max(1, std_degree / 10)), 5))
    degree_bins = np.linspace(min(degrees), max_degree, num_bins + 1)
    
    # Calculate 95th percentile for zoom
    percentile_95 = np.percentile(degrees, 95)
    
    # --- PLOT 1: Overall degree distribution ---
    ax1.hist(degrees, bins=degree_bins, edgecolor='black', color='steelblue', alpha=0.7)
    ax1.set_xlabel('Vertex Degree', fontsize=11)
    ax1.set_ylabel('Number of Vertices', fontsize=11)
    ax1.set_title('Degree Distribution (All Vertices)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axvline(np.mean(degrees), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(degrees):.2f}')
    ax1.axvline(np.median(degrees), color='green', linestyle='--', linewidth=2, label=f'Median: {np.median(degrees):.2f}')
    ax1.legend(fontsize=9)
    
    # --- PLOT 1 ZOOM: 95% quartile of degree distribution ---
    degrees_95 = [d for d in degrees if d <= percentile_95]
    zoom_bins = np.linspace(min(degrees_95), max(degrees_95), min(30, max(len(set(degrees_95)), 5)) + 1)
    ax1_zoom.hist(degrees_95, bins=zoom_bins, edgecolor='black', color='steelblue', alpha=0.7)
    ax1_zoom.set_xlabel('Vertex Degree', fontsize=11)
    ax1_zoom.set_ylabel('Number of Vertices', fontsize=11)
    ax1_zoom.set_title(f'Degree Distribution - 95% Quartile (≤{percentile_95:.1f})', fontsize=12, fontweight='bold')
    ax1_zoom.grid(True, alpha=0.3)
    ax1_zoom.axvline(np.mean(degrees_95), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(degrees_95):.2f}')
    ax1_zoom.axvline(np.median(degrees_95), color='green', linestyle='--', linewidth=2, label=f'Median: {np.median(degrees_95):.2f}')
    ax1_zoom.legend(fontsize=9)
    
    # --- PLOT 2: Degree distribution by file ---
    unique_filenames = sorted(set(vertex_attrs[v]['filename'] for v in range(graph.vcount())))
    palette = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F', '#BB8FCE', '#85C1E2', '#F8B88B', '#AED6F1']
    
    for idx, filename in enumerate(unique_filenames):
        file_degrees = [graph.degree(v) for v in range(graph.vcount()) if vertex_attrs[v]['filename'] == filename]
        color = palette[idx % len(palette)]
        from pathlib import Path
        file_label = Path(filename).stem
        ax2.hist(file_degrees, bins=degree_bins, alpha=0.5, label=file_label, color=color, edgecolor='black')
    
    ax2.set_xlabel('Vertex Degree', fontsize=11)
    ax2.set_ylabel('Number of Vertices', fontsize=11)
    ax2.set_title('Degree Distribution (by Source File)', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=8, loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # --- PLOT 2 ZOOM: 95% quartile by file ---
    for idx, filename in enumerate(unique_filenames):
        file_degrees = [graph.degree(v) for v in range(graph.vcount()) if vertex_attrs[v]['filename'] == filename]
        file_degrees_95 = [d for d in file_degrees if d <= percentile_95]
        if file_degrees_95:
            color = palette[idx % len(palette)]
            file_label = Path(filename).stem
            ax2_zoom.hist(file_degrees_95, bins=zoom_bins, alpha=0.5, label=file_label, color=color, edgecolor='black')
    
    ax2_zoom.set_xlabel('Vertex Degree', fontsize=11)
    ax2_zoom.set_ylabel('Number of Vertices', fontsize=11)
    ax2_zoom.set_title(f'Degree Distribution by File - 95% Quartile (≤{percentile_95:.1f})', fontsize=12, fontweight='bold')
    ax2_zoom.legend(fontsize=8, loc='upper right')
    ax2_zoom.grid(True, alpha=0.3)
    
    # --- PLOT 3: Components by size distribution (with log Y-axis, minimal labels) ---
    if graph_composition and 'composition_summary' in graph_composition:
        comp_by_size = graph_composition['composition_summary'].get('components_by_size', {})
        if comp_by_size:
            # Sort by size numerically
            sizes = sorted([int(k) for k in comp_by_size.keys()])
            counts = [comp_by_size[str(s)] for s in sizes]
            
            ax3.bar(range(len(sizes)), counts, color='darkorange', alpha=0.7, edgecolor='black')
            ax3.set_xlabel('Component Size (# vertices)', fontsize=11)
            ax3.set_ylabel('Number of Components', fontsize=11)
            ax3.set_title('Components by Size Distribution (Log Scale)', fontsize=12, fontweight='bold')
            ax3.set_yscale('log')  # Log scale for Y-axis
            ax3.grid(True, alpha=0.3, axis='y')
            
            # Reduce X-axis label clutter: show every nth label based on number of sizes
            tick_interval = max(1, len(sizes) // 10)  # Show approx 10 labels max
            ax3.set_xticks(range(0, len(sizes), tick_interval))
            ax3.set_xticklabels([str(sizes[i]) for i in range(0, len(sizes), tick_interval)], fontsize=9)
    
    # --- PLOT 6: Components by size distribution (with linear Y-axis) ---
    if graph_composition and 'composition_summary' in graph_composition:
        comp_by_size = graph_composition['composition_summary'].get('components_by_size', {})
        if comp_by_size:
            # Sort by size numerically
            sizes = sorted([int(k) for k in comp_by_size.keys()])
            counts = [comp_by_size[str(s)] for s in sizes]
            
            ax3_linear.bar(range(len(sizes)), counts, color='darkorange', alpha=0.7, edgecolor='black')
            ax3_linear.set_xlabel('Component Size (# vertices)', fontsize=11)
            ax3_linear.set_ylabel('Number of Components', fontsize=11)
            ax3_linear.set_title('Components by Size Distribution (Linear Scale)', fontsize=12, fontweight='bold')
            ax3_linear.grid(True, alpha=0.3, axis='y')
            
            # Reduce X-axis label clutter: show every nth label based on number of sizes
            tick_interval = max(1, len(sizes) // 10)  # Show approx 10 labels max
            ax3_linear.set_xticks(range(0, len(sizes), tick_interval))
            ax3_linear.set_xticklabels([str(sizes[i]) for i in range(0, len(sizes), tick_interval)], fontsize=9)
    
    fig.suptitle(f'Network Graph Analysis ({graph.vcount()} vertices, {graph.ecount()} edges)', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved network analysis histograms: {output_path}")
    return True


def analyze_graph_composition(graph: ig.Graph, vertex_attrs: Dict) -> Dict:
    """
    Analyze graph composition: how many connected components have how many features and files.
    
    Returns a dictionary with statistics like:
    {
        "component_composition": {
            "3_features_2_files": 5,  # 5 components with 3 features from 2 files
            "4_features_3_files": 2,  # 2 components with 4 features from 3 files
            ...
        },
        "composition_summary": {
            "total_components": 10,
            "total_features": 35,
            "avg_features_per_component": 3.5,
            "avg_files_per_component": 2.1,
            "max_features_in_component": 10,
            "max_files_in_component": 3,
            "components_by_size": {
                "1": 2,  # 2 components with 1 feature
                "2": 5,  # 5 components with 2 features
                "3": 3,  # 3 components with 3 features
            }
        }
    }
    """
    if graph is None or graph.vcount() == 0:
        return {'component_composition': {}, 'composition_summary': {}}
    
    components = graph.components()
    composition_map = {}  # Key: f"{num_features}_{num_files}", Value: count
    features_per_component = []
    files_per_component = []
    components_by_size = {}  # Key: num_features, Value: count
    
    for component_vertices in components:
        # Count unique files in this component
        unique_files = set()
        for vertex_id in component_vertices:
            filename = vertex_attrs[vertex_id]['filename']
            unique_files.add(filename)
        
        num_features = len(component_vertices)
        num_files = len(unique_files)
        
        # Create composition key
        comp_key = f"{num_features}_{num_files}"
        composition_map[comp_key] = composition_map.get(comp_key, 0) + 1
        
        features_per_component.append(num_features)
        files_per_component.append(num_files)
        
        # Track components by size
        size_key = str(num_features)
        components_by_size[size_key] = components_by_size.get(size_key, 0) + 1
    
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
        'composition_summary': summary
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
                       unpaired_singletons_count: int = 0,
                       num_cores: int = None,
                       use_multiprocessing: bool = True,
                       initial_seed: int = 42,
                       mega_complexity_threshold: float = 10000) -> Dict:
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
        unpaired_singletons_count: Number of unpaired singleton vertices to include in avg_cluster_size calculation
    
    Returns:
        Dictionary with metrics: {
            'resolution': float,
            'vertices_deleted': int,
            'avg_cluster_size': float,
            'score': float
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
    
    # Include unpaired singletons AND unrecovered deleted vertices in average cluster size calculation
    # Both become singleton (size 1) clusters
    total_singletons = unpaired_singletons_count + unrecovered_deleted
    total_size_sum = sum(all_cluster_sizes) + total_singletons  
    total_cluster_count = len(all_cluster_sizes) + total_singletons 
    
    # DEBUG: Print singleton and cluster count metrics
    print(f"[DEBUG] total_recovered={total_recovered}, deleted_vertices={deleted_vertices}, unrecovered_deleted={unrecovered_deleted}, total_singletons={total_singletons}, total_size_sum={total_size_sum}, total_cluster_count={total_cluster_count}")
    
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
    
    # Add singleton clusters (unpaired + unrecovered) to the distribution
    if total_singletons > 0:
        cluster_recovered_final[1] = cluster_recovered_final.get(1, 0) + total_singletons
    
    # Note: lambda will be passed in main function for scoring
    # For now, just return the metrics
    
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
        # Include clustering info for reuse in final processing (avoid re-clustering)
        '_clustering_data': {
            'vertex_to_cluster_map': vertex_to_cluster_map,
            'component_cluster_info': component_cluster_info,
            'filtered_vertex_map': filtered_vertex_map if filtered_vertex_map else {},
            'filtered_cluster_info': filtered_cluster_info if filtered_cluster_info else {}
        }
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
                             unpaired_singletons_count: int = 0,
                             num_cores: int = None,
                             initial_seed: int = 42,
                             mega_complexity_threshold: float = 10000) -> Tuple[float, Dict, List]:
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
    
    # Storage for iteration history
    iteration_history = []
    call_count = [0]  # Use list to allow modification in nested function
    
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
            unpaired_singletons_count=unpaired_singletons_count,
            num_cores=num_cores,
            initial_seed=initial_seed,
            mega_complexity_threshold=mega_complexity_threshold
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
            'cluster_recovered_final': metrics.get('cluster_recovered_final', {}),
            '_clustering_data': metrics.get('_clustering_data', {})  # Store clustering for reuse
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
    best_clustering_data = best_iteration.get('_clustering_data', {})
    
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
    
    return best_resolution, best_metrics, iteration_history, best_clustering_data


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
                          unpaired_singletons_count: int = 0,
                          num_cores: int = None,
                          initial_seed: int = 42,
                          mega_complexity_threshold: float = 10000) -> Tuple[float, Dict, List]:
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
    best_clustering_data = None  # Store clustering from best iteration to avoid re-clustering
    
    for iteration, current_resolution in enumerate(resolution_values):
        print(f"  [{iteration+1:2d}/{num_points}] Evaluating resolution={current_resolution:.4f}...", end='', flush=True)
        
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
            unpaired_singletons_count=unpaired_singletons_count,
            num_cores=num_cores,
            initial_seed=initial_seed,
            mega_complexity_threshold=mega_complexity_threshold
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
            'iteration': iteration + 1,
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
            best_clustering_data = metrics.get('_clustering_data', {})  # Save clustering from best iteration
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
    
    return best_resolution, best_metrics, iteration_history, best_clustering_data


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
                    'charge': vertex_info.get('charge', 0)
                })
            
            clusters_output[comp_key]['clusters'][cluster_key] = {
                'cluster_id': cluster_id,
                'num_vertices': len(vertex_ids),
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



def _clusterwise_modularity_heuristic(comp_idx: int, clusters_with_duplicates: Dict, 
                                       comp_clusters: Dict, comp_vertices: List,
                                       vertex_attrs: Dict, vertex_to_cluster_map: Dict,
                                       graph_edges: List) -> tuple:
    """
    Clusterwise modularity optimization heuristic for mega-components.
    
    For each cluster, keeps only ONE vertex per file (the one with highest modularity score).
    Score = internal_degree - external_degree
      - internal_degree: edges to vertices in same cluster
      - external_degree: edges to vertices in different clusters (same component)
    
    Returns:
        (vertices_to_keep: frozenset, deleted_vertices_list: list)
    """
    print(f"[COMP {comp_idx}] ⚠  USING HEURISTIC (complexity too high for exhaustive search)", flush=True)
    
    # Build neighbor map for degree calculations
    neighbor_map = {}
    for src, tgt in graph_edges:
        if src not in neighbor_map:
            neighbor_map[src] = []
        if tgt not in neighbor_map:
            neighbor_map[tgt] = []
        neighbor_map[src].append(tgt)
        neighbor_map[tgt].append(src)
    
    vertices_to_keep = set(comp_vertices)  # Start with all, then remove selected ones
    deleted_vertices_list = []
    
    # For each cluster with duplicates
    for cluster_id, file_vertex_map in clusters_with_duplicates.items():
        # file_vertex_map: {filename: [vertex_ids]}
        
        for filename, vertex_ids in file_vertex_map.items():
            if len(vertex_ids) <= 1:
                continue  # No duplicates in this file within this cluster
            
            # Score each vertex in this file group
            scores = {}
            for vertex_id in vertex_ids:
                internal_degree = 0
                external_degree = 0
                
                for neighbor in neighbor_map.get(vertex_id, []):
                    neighbor_cluster = vertex_to_cluster_map[neighbor][1]
                    if neighbor_cluster == cluster_id:
                        internal_degree += 1
                    else:
                        external_degree += 1
                
                score = internal_degree - external_degree
                scores[vertex_id] = score
            
            # Keep only the highest-scoring vertex for this file in this cluster
            best_vertex = max(vertex_ids, key=lambda v: scores[v])
            
            # Remove all others from vertices_to_keep
            for vertex_id in vertex_ids:
                if vertex_id != best_vertex:
                    vertices_to_keep.discard(vertex_id)
                    comp_id, clust_id = vertex_to_cluster_map[vertex_id]
                    deleted_vertices_list.append({
                        'component_id': comp_id,
                        'cluster_id': clust_id,
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
                        'reason': 'Duplicate file vertex - filtered via clusterwise heuristic',
                        'heuristic_score': float(scores[vertex_id])
                    })
    
    print(f"[COMP {comp_idx}] Heuristic deleted {len(deleted_vertices_list)} duplicates, keeping {len(vertices_to_keep)} vertices", flush=True)
    return frozenset(vertices_to_keep), deleted_vertices_list


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
            - graph_edges: List of edges in this component's subgraph
            - mega_complexity_threshold: Threshold for using heuristic
    
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
        graph_edges = component_data['graph_edges']
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
        
        # Reconstruct graph edges for this component for degree calculations
        # graph_edges is list of (src, tgt) tuples for edges in this component's subgraph
        neighbor_map = {}  # vertex_id -> list of neighbors in this component
        
        for src, tgt in graph_edges:
            if src not in neighbor_map:
                neighbor_map[src] = []
            if tgt not in neighbor_map:
                neighbor_map[tgt] = []
            neighbor_map[src].append(tgt)
            neighbor_map[tgt].append(src)  # Undirected
        
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
            scoring_path_counts = {'induced_subgraph': 0, 'fallback_from_exception': 0, 'neighbor_map_no_graph': 0, 
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
                    # Single cluster: score by total degree
                    total_degree = sum(len(neighbor_map.get(v, [])) for v in test_vertices)
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
                                # Fallback if induced_subgraph fails - use neighbor_map approach
                                scoring_path = f'fallback_from_exception: {type(e).__name__}'
                                weighted_score = 0
                                for v in test_vertices:
                                    v_cluster = vertex_to_cluster_map[v][1]
                                    for neighbor in neighbor_map.get(v, []):
                                        if neighbor in candidate:
                                            neighbor_cluster = vertex_to_cluster_map[neighbor][1]
                                            if neighbor_cluster == v_cluster:
                                                weighted_score += 1
                                            else:
                                                weighted_score -= 1
                                score = weighted_score
                        else:
                            # No graph available, use neighbor_map
                            scoring_path = 'neighbor_map_no_graph'
                            weighted_score = 0
                            for v in test_vertices:
                                v_cluster = vertex_to_cluster_map[v][1]
                                for neighbor in neighbor_map.get(v, []):
                                    if neighbor in candidate:
                                        neighbor_cluster = vertex_to_cluster_map[neighbor][1]
                                        if neighbor_cluster == v_cluster:
                                            weighted_score += 1
                                        else:
                                            weighted_score -= 1
                            score = weighted_score
                        
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
                                # Fallback to weighted degree if modularity fails
                                scoring_path = f'modularity_fallback: {type(e).__name__}'
                                weighted_score = 0
                                for v in test_vertices:
                                    v_cluster = vertex_to_cluster_map[v][1]
                                    for neighbor in neighbor_map.get(v, []):
                                        if neighbor in candidate:
                                            neighbor_cluster = vertex_to_cluster_map[neighbor][1]
                                            if neighbor_cluster == v_cluster:
                                                weighted_score += 1
                                            else:
                                                weighted_score -= 1
                                score = weighted_score
                        else:
                            # No graph available, use neighbor_map
                            scoring_path = 'modularity_no_graph_fallback'
                            weighted_score = 0
                            for v in test_vertices:
                                v_cluster = vertex_to_cluster_map[v][1]
                                for neighbor in neighbor_map.get(v, []):
                                    if neighbor in candidate:
                                        neighbor_cluster = vertex_to_cluster_map[neighbor][1]
                                        if neighbor_cluster == v_cluster:
                                            weighted_score += 1
                                        else:
                                            weighted_score -= 1
                            score = weighted_score
                        
                        # Track which scoring path was used
                        if scoring_path in scoring_path_counts:
                            scoring_path_counts[scoring_path] += 1
                    else:
                        # Fallback to degree for unknown filter methods
                        score = sum(len(neighbor_map.get(v, [])) for v in test_vertices)
                
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
        
        # Extract edges using pre-built edge index
        component_vertex_set = set(comp_vertices)
        component_edges = []
        for vertex_id in comp_vertices:
            for neighbor in edge_index.get(vertex_id, []):
                if neighbor > vertex_id and neighbor in component_vertex_set:
                    component_edges.append((vertex_id, neighbor))
        
        # Create component-scoped copies
        comp_vertex_attrs = {v: vertex_attrs[v] for v in comp_vertices}
        comp_vertex_to_cluster_map = {v: vertex_to_cluster_map[v] for v in comp_vertices}
        
        # Estimate memory usage for this component
        # Each work_item contains: 
        #   - comp_vertices list: ~8 bytes per vertex
        #   - comp_clusters dict: ~50 bytes per entry
        #   - comp_vertex_attrs dict: ~100 bytes per entry
        #   - component_edges list: ~16 bytes per edge
        comp_vertex_count = len(comp_vertices)
        comp_cluster_count = len(comp_clusters)
        comp_edge_count = len(component_edges)
        estimated_item_size_bytes = (
            comp_vertex_count * 8 +
            comp_cluster_count * 50 +
            comp_vertex_count * 100 +
            comp_edge_count * 16
        )
        
        # Log large components
        if estimated_item_size_bytes > 50_000_000:  # >50MB
            print(f"[PID {worker_pid}] ⚠ LARGE COMPONENT DETECTED: comp_idx={comp_idx}, "
                  f"vertices={comp_vertex_count}, clusters={comp_cluster_count}, edges={comp_edge_count}, "
                  f"est_size={estimated_item_size_bytes / 1_000_000:.1f}MB", flush=True)
        
        return {
            'comp_idx': comp_idx,
            'comp_vertices': comp_vertices,
            'comp_clusters': comp_clusters,
            'vertex_attrs': comp_vertex_attrs,
            'vertex_to_cluster_map': comp_vertex_to_cluster_map,
            'info': info,
            'filter_method': filter_method,
            'graph_edges': component_edges,
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
    
    # Return generator that yields results as they're produced (not accumulated)
    return result_generator()


def _evaluate_candidates_batch(candidates_batch, comp_idx, batch_num, comp_vertices, 
                               neighbor_map, clusters_with_duplicates, filter_method, 
                               vertex_to_cluster_map):
    """
    Evaluate a batch of candidates and return the best one in this batch.
    
    Args:
        candidates_batch: List of candidate sets to evaluate
        comp_idx: Component index (for logging)
        batch_num: Batch number (for logging)
        comp_vertices: All vertices in component
        neighbor_map: vertex_id -> list of neighbors
        clusters_with_duplicates: {cluster_id: {filename: [vertex_ids]}}
        filter_method: Filtering method name
        vertex_to_cluster_map: {vertex_id: (comp_idx, cluster_id)}
    
    Returns:
        Tuple of (best_candidate, best_score) for this batch
    """
    batch_best_candidate = None
    batch_best_score = None
    
    for idx, candidate in enumerate(candidates_batch):
        test_vertices = [v for v in comp_vertices if v in candidate]
        if len(test_vertices) <= 1:
            continue
        
        # Score this candidate
        if len(clusters_with_duplicates) == 1:
            # Single cluster: score by total degree
            total_degree = sum(len(neighbor_map.get(v, [])) for v in test_vertices)
            score = total_degree
        else:
            # Multiple clusters: use weighted degree
            if filter_method == 'simplified-modularity':
                weighted_score = 0
                for v in test_vertices:
                    v_cluster = vertex_to_cluster_map[v][1]
                    for neighbor in neighbor_map.get(v, []):
                        if neighbor in candidate:
                            neighbor_cluster = vertex_to_cluster_map[neighbor][1]
                            if neighbor_cluster == v_cluster:
                                weighted_score += 1  # Internal edge weight - match linear version (was 2)
                            else:
                                weighted_score -= 1  # External edge weight
                score = weighted_score
            else:
                # Fallback to degree
                score = sum(len(neighbor_map.get(v, [])) for v in test_vertices)
        
        # Track batch best
        if batch_best_score is None or score > batch_best_score:
            batch_best_score = score
            batch_best_candidate = candidate
    
    if batch_best_candidate is not None:
        print(f"[COMP {comp_idx}] Batch {batch_num}: best score = {batch_best_score}", flush=True)
    
    return (batch_best_candidate, batch_best_score)



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
COMBINATIONS_BATCH_LIMIT = 100_000  # Increased to 5 million for worker-side batching (reduce IPC overhead)
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
            num_workers = min(max(2, 6), cores_to_use)  # Use 8 workers if available, min 2, max available cores
            
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
    for cluster_id in clusters_with_duplicates:
        file_vertex_map = clusters_with_duplicates[cluster_id]
        num_dup_options = sum(1 for vids in file_vertex_map.values() if len(vids) > 1)
        if num_dup_options > 0:
            # Estimate average options per duplicate set
            avg_option_size = sum(len(vids) for vids in file_vertex_map.values() if len(vids) > 1) / max(1, num_dup_options)
            combo_estimate *= (max(2, int(avg_option_size)) ** num_dup_options)
            # Track heuristic decision WITHOUT breaking - let complexity calculation complete
            if combo_estimate >= mega_complexity_threshold:
                use_heuristic = True  # Will use heuristic if true at end of loop
            # Cap for exponential safety (avoid inf/overflow) but don't break early
            if combo_estimate > 1e20:
                combo_estimate = 1e20
    
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


def get_cluster_colors(num_clusters: int) -> List[str]:
    """Generate distinct colors for clusters."""
    # Color palette for clusters
    palette = [
        '#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8',
        '#F7DC6F', '#BB8FCE', '#85C1E2', '#F8B88B', '#AED6F1',
        '#F5B041', '#82E0AA', '#F1948A', '#85C1E2', '#D7BDE2',
        '#F8B739', '#6C5B7B', '#FA7921', '#87E8D5', '#FF4757',
        '#1ABC9C', '#3498DB', '#9B59B6', '#E74C3C', '#E67E22',
    ]
    
    colors = []
    for i in range(num_clusters):
        colors.append(palette[i % len(palette)])
    
    return colors


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


def visualize_graph(graph: ig.Graph, vertex_attrs: Dict, output_image: str,
                   cluster_map: Dict = None, component_cluster_info: Dict = None,
                   layout_engine: str = 'dot',
                   use_cluster_subgraphs: bool = False,
                   layout_use_weights: bool = False,
                   use_coordinate_layout: bool = False,
                   vertex_coordinates: Dict = None,
                   deleted_vertex_ids: set = None):
    """Generate and save network graph visualization using Graphviz (with optimizations).
    
    Args:
        use_coordinate_layout: Use original feature coordinates as initial positions
        vertex_coordinates: Dict mapping vertex_id -> (x, y) coordinates from heatmap
    """
    if not GRAPHVIZ_AVAILABLE:
        print("✗ graphviz Python package not available for visualization")
        print("  Install with: pip install graphviz")
        return False
    
    try:
        start_time = time.time()
        
        # Determine output format from file extension
        output_format = 'svg'  # SVG is much faster to render than PNG
        if output_image.endswith('.pdf'):
            output_format = 'pdf'
        elif output_image.endswith('.png'):
            output_format = 'png'
        
        # Use a temporary path without dots (to avoid Graphviz confusion with filename containing dots)
        # Create temp path in same directory as output_image
        output_dir = str(Path(output_image).parent)
        temp_name = f"_temp_graph_{int(time.time() * 1000000) % 1000000}"
        if output_dir and output_dir != '.':
            temp_path = str(Path(output_dir) / temp_name)
        else:
            temp_path = temp_name

        # Calculate degree-based node sizes
        degrees = graph.degree()
        min_degree = min(degrees) if degrees else 1
        max_degree = max(degrees) if degrees else 1
        degree_range = max_degree - min_degree if max_degree > min_degree else 1
        
        num_vertices = graph.vcount()
        num_edges = graph.ecount()
        print(f"\nRendering {num_vertices} vertices and {num_edges} edges...")
        print(f"  Output format: {output_format} (SVG is fastest, ~2x faster than PNG)")
        print(f"  Layout engine: {layout_engine}")
        
        # Estimate rendering time based on graph size
        estimated_format_time = (num_vertices * 0.0001 + num_edges * 0.00005) / 60  # in minutes
        estimated_render_time = (num_edges * 0.003) / 60  # in minutes
        estimated_total = estimated_format_time + estimated_render_time
        
        print(f"  Estimated format time: ~{estimated_format_time:.1f} min")
        print(f"  Estimated render time: ~{estimated_render_time:.1f} min")
        print(f"  Estimated total time: ~{estimated_total:.1f} min")
        print(f"\n  Note: Graphviz dot layout is single-threaded. For faster rendering:")
        print(f"    - Use SVG format (--output_image graph.svg)")
        print(f"    - Filter by edge weight (--edge_cutoff <distance>)")
        print(f"    - Use simpler layout: graph.gml or graph.graphml format instead")
        
        # Validate layout engine compatibility for weighted layouts
        if layout_use_weights and layout_engine not in ['neato', 'fdp', 'sfdp']:
            print(f"✗ Layout weights are only supported for neato, fdp, sfdp (got: {layout_engine})")
            return False
        
        # Validate layout engine compatibility for coordinate-based layouts
        if use_coordinate_layout and layout_engine not in ['neato', 'fdp']:
            print(f"✗ Coordinate-based layout is only supported for neato, fdp (got: {layout_engine})")
            return False
        
        if use_coordinate_layout and not vertex_coordinates:
            print(f"✗ Coordinate-based layout requested but no coordinates provided")
            return False

        # Create graph
        dot = graphviz.Digraph(comment='Feature Network', format=output_format, engine=layout_engine)
        
        # Optimize DOT attributes for speed and label placement
        if layout_engine in ['sfdp', 'fdp', 'neato']:
            dot.attr(overlap='false', splines='true', sep='+1.0', pad='0.3')
        else:
            dot.attr(rankdir='LR', splines='polyline', overlap='compress', sep='+0.5', pad='0.2')
        dot.attr('node', shape='circle', style='filled', fontsize='1.5', margin='0.01', width='0.1', height='0.1')
        dot.attr('edge', fontsize='1')
        
        print(f"\n  [1/3] Adding {num_vertices} vertices...", end='', flush=True)
        vertex_time_start = time.time()
        
        # Precompute cluster colors if clustering is enabled
        cluster_color_map = {}
        if cluster_map and component_cluster_info:
            cluster_keys = sorted(set(cluster_map.values()))
            cluster_colors = get_cluster_colors(len(cluster_keys))
            cluster_color_map = {key: cluster_colors[i] for i, key in enumerate(cluster_keys)}

        # Build Graphviz cluster subgraphs when requested
        cluster_subgraphs = {}
        if use_cluster_subgraphs and cluster_map:
            for cluster_key in sorted(set(cluster_map.values())):
                cluster_name = f"cluster_{cluster_key[0]}_{cluster_key[1]}"
                subgraph = graphviz.Digraph(name=cluster_name)
                subgraph.attr(label=f"C{cluster_key[0]}:{cluster_key[1]}",
                              style='dashed', color='#999999')
                cluster_subgraphs[cluster_key] = subgraph

        # Add vertices (nodes) with small labels
        for vertex_id in range(num_vertices):
            # Base color from vertex attributes (assigned based on source file)
            if vertex_attrs and vertex_id in vertex_attrs and 'color' in vertex_attrs[vertex_id]:
                base_color = vertex_attrs[vertex_id]['color']
            else:
                # Fallback color (should not happen)
                base_color = '#87CEEB'

            # Cluster overlay color (fill) if available
            fill_color = base_color
            if cluster_map and vertex_id in cluster_map:
                cluster_key = cluster_map[vertex_id]
                fill_color = cluster_color_map.get(cluster_key, base_color)
            
            # Get the label from vertex_attrs
            label = vertex_attrs[vertex_id]['label'] if vertex_id in vertex_attrs else str(vertex_id)
            
            # Build node attributes
            node_attrs = {'label': label, 'color': base_color, 'fillcolor': fill_color}
            
            # Check if vertex is recovered (was deleted but re-added to another cluster)
            is_recovered = vertex_attrs[vertex_id].get('recovered', False) if vertex_id in vertex_attrs else False
            
            # Apply special styling for recovered vertices
            if is_recovered:
                node_attrs['penwidth'] = '4.0'  # Thick border
                node_attrs['color'] = '#FF1493'  # Deep pink border
                node_attrs['style'] = 'filled'
            # Apply reduced opacity if vertex was deleted during filtering and not recovered
            elif deleted_vertex_ids and vertex_id in deleted_vertex_ids:
                node_attrs['penwidth'] = '0.5'
                # Make the color lighter/more transparent by blending with white
                # Convert hex color to RGB, blend with white, convert back
                try:
                    # Parse hex color
                    if fill_color.startswith('#'):
                        hex_color = fill_color[1:]
                    else:
                        hex_color = fill_color
                    
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    # Blend with white (50% opacity = blend 50% with white)
                    r = int(r * 0.5 + 255 * 0.5)
                    g = int(g * 0.5 + 255 * 0.5)
                    b = int(b * 0.5 + 255 * 0.5)
                    lighter_color = f"#{r:02x}{g:02x}{b:02x}"
                    node_attrs['fillcolor'] = lighter_color
                except:
                    # Fallback: use a light gray
                    node_attrs['fillcolor'] = '#CCCCCC'
                
                # Also set a dashed style to make it visually distinct
                node_attrs['style'] = 'dashed,filled'
            
            # Add initial position if coordinate layout is enabled
            if use_coordinate_layout and vertex_coordinates and vertex_id in vertex_coordinates:
                x, y = vertex_coordinates[vertex_id]
                # Graphviz pos format: "x,y!" (exclamation mark pins the position)
                # Without !, it's used as a starting point but can be adjusted
                node_attrs['pos'] = f"{x},{y}"
            
            # Add node with label visible and assigned color
            if use_cluster_subgraphs and cluster_map and vertex_id in cluster_map:
                cluster_key = cluster_map[vertex_id]
                subgraph = cluster_subgraphs.get(cluster_key, dot)
                subgraph.node(str(vertex_id), **node_attrs)
            else:
                dot.node(str(vertex_id), **node_attrs)
            
            # Progress every 50k vertices
            if (vertex_id + 1) % 50000 == 0:
                elapsed = time.time() - vertex_time_start
                rate = (vertex_id + 1) / elapsed
                remaining = (num_vertices - vertex_id - 1) / rate / 60
                print(f"\n    {vertex_id + 1}/{num_vertices} ({100*(vertex_id+1)/num_vertices:.1f}%) - ETA: {remaining:.1f} min", end='', flush=True)
        
        # Attach cluster subgraphs to the main graph after nodes are added
        if cluster_subgraphs:
            for subgraph in cluster_subgraphs.values():
                dot.subgraph(subgraph)

        vertex_elapsed = time.time() - vertex_time_start
        print(f" ✓ ({vertex_elapsed:.1f}s)")
        
        # Add all edges with weight labels
        edge_time_start = time.time()
        print(f"  [2/3] Adding {num_edges} edges...", end='', flush=True)
        
        # Validate a sample of edges to confirm correct direction
        sample_edges = min(5, num_edges)
        print(f"\n    Validating edge directions (sample of {sample_edges} edges):")
        for i in range(sample_edges):
            edge = graph.es[i]
            src_vertex = graph.vs[edge.source]
            tgt_vertex = graph.vs[edge.target]
            src_label = src_vertex['label'] if src_vertex else str(edge.source)
            tgt_label = tgt_vertex['label'] if tgt_vertex else str(edge.target)
            print(f"      Edge {i}: {src_label} → {tgt_label} (distance: {edge['distance']:.3f})")
        print()
        
        # Build edge pair map for bidirectional detection and deduplication
        edge_pairs = {}  # (min, max) -> edges in that pair
        directed_edge_count = {}  # (src, tgt) -> count of edges with same direction
        
        for edge in graph.es:
            src, tgt = edge.source, edge.target
            pair_key = tuple(sorted([src, tgt]))
            if pair_key not in edge_pairs:
                edge_pairs[pair_key] = {'forward': None, 'backward': None, 'forward_list': [], 'backward_list': []}
            
            # Track edge direction
            if src == pair_key[0] and tgt == pair_key[1]:
                edge_pairs[pair_key]['forward'] = edge
                edge_pairs[pair_key]['forward_list'].append(edge)
            else:
                edge_pairs[pair_key]['backward'] = edge
                edge_pairs[pair_key]['backward_list'].append(edge)
            
            # Count directed edge duplicates
            directed_key = (src, tgt)
            directed_edge_count[directed_key] = directed_edge_count.get(directed_key, 0) + 1
        
        # Count bidirectional edges and find examples
        bidirectional_count = 0
        bidirectional_examples = []
        duplicate_directed_examples = []
        
        for edge_dict in edge_pairs.values():
            if edge_dict['forward'] is not None and edge_dict['backward'] is not None:
                bidirectional_count += 1
                if len(bidirectional_examples) < 3:  # Collect first 3 examples
                    fwd = edge_dict['forward']
                    bwd = edge_dict['backward']
                    src_v = graph.vs[fwd.source]
                    tgt_v = graph.vs[fwd.target]
                    bidirectional_examples.append({
                        'fwd_src': src_v['label'],
                        'fwd_tgt': tgt_v['label'],
                        'fwd_dist': fwd['distance'],
                        'bwd_dist': bwd['distance'],
                        'fwd_count': len(edge_dict['forward_list']),
                        'bwd_count': len(edge_dict['backward_list'])
                    })
        
        # Find duplicate directed edges (A→B appears multiple times)
        for (src, tgt), count in directed_edge_count.items():
            if count > 1 and len(duplicate_directed_examples) < 3:
                src_v = graph.vs[src]
                tgt_v = graph.vs[tgt]
                # Get the distances
                edges_with_direction = [e for e in graph.es if e.source == src and e.target == tgt]
                distances = [e['distance'] for e in edges_with_direction]
                duplicate_directed_examples.append({
                    'src': src_v['label'],
                    'tgt': tgt_v['label'],
                    'count': count,
                    'distances': distances
                })
        
        print(f"  Bidirectional edges detected: {bidirectional_count}/{len(edge_pairs)}")
        if bidirectional_examples:
            print(f"  Examples of bidirectional edge pairs:")
            for i, ex in enumerate(bidirectional_examples):
                print(f"    {i+1}. {ex['fwd_src']} ↔ {ex['fwd_tgt']}")
                print(f"       → distance: {ex['fwd_dist']:.3f}, ← distance: {ex['bwd_dist']:.3f}")
        
        if duplicate_directed_examples:
            print(f"  Examples of duplicate directed edges (same A→B appears multiple times):")
            for i, ex in enumerate(duplicate_directed_examples):
                print(f"    {i+1}. {ex['src']} → {ex['tgt']} ({ex['count']} copies)")
                print(f"       distances: {', '.join(f'{d:.3f}' for d in ex['distances'])}")
        
        edges_added = 0
        bidirectional_rendered = 0
        
        for pair_key, edge_dict in edge_pairs.items():
            # Use forward edge if available, otherwise backward
            edge_to_use = edge_dict['forward'] if edge_dict['forward'] else edge_dict['backward']
            
            src, tgt = edge_to_use.source, edge_to_use.target
            distance = edge_to_use['distance'] if 'distance' in edge_to_use.attributes() else 1.0
            
            # Use the pre-computed bidirectional flag from the full dataset
            attrs = edge_to_use.attributes()
            is_bidirectional = attrs.get('is_bidirectional', False) if attrs else False
            if is_bidirectional:
                bidirectional_rendered += 1
            
            # Format distance for label
            distance_label = f"{distance:.3f}"
            
            # Color based on distance: closer to 0 = darker/stronger connection
            alpha = int(255 * distance / (1.0 + distance))
            alpha = max(0, min(255, alpha))

            # Use a darker outline to keep very light edges visible
            outline_alpha = max(0, alpha - 40)
            
            # If bidirectional, use bidirectional arrow; otherwise use single direction
            arrow_dir = 'both' if is_bidirectional else 'forward'
            
            # Add edge with weight label
            edge_attrs = {
                'label': distance_label,
                'color': f'#{outline_alpha:02x}{outline_alpha:02x}{outline_alpha:02x}',
                'fillcolor': f'#{alpha:02x}{alpha:02x}{alpha:02x}',
                'style': 'filled',
                'penwidth': '0.5',
                'arrowsize': '0.3',
                'dir': arrow_dir
            }
            if layout_use_weights:
                edge_attrs['len'] = f"{max(0.1, float(distance)):.3f}"
            dot.edge(str(src), str(tgt), **edge_attrs)
            edges_added += 1
            
            # Progress every 10k edge pairs
            if (edges_added + 1) % 10000 == 0:
                elapsed = time.time() - edge_time_start
                rate = edges_added / elapsed
                remaining = (len(edge_pairs) - edges_added) / rate / 60 if rate > 0 else 0
                print(f"\n    {edges_added}/{len(edge_pairs)} ({100*edges_added/len(edge_pairs):.1f}%) - ETA: {remaining:.1f} min", end='', flush=True)
        
        edge_elapsed = time.time() - edge_time_start
        duplicates_removed = num_edges - edges_added
        print(f" ✓ ({edge_elapsed:.1f}s, {duplicates_removed} duplicates removed, {bidirectional_rendered} bidirectional arrows)")
        
        # Add graph title
        dot.attr('graph', label=f'Feature Network ({num_vertices} vertices, {edges_added} unique edge pairs, {bidirectional_rendered}↔)')
        
        # Render and save
        print(f"  [3/3] Rendering with {layout_engine} (using {output_format} format)...")
        render_time_start = time.time()
        
        # Render to temporary file using clean name (avoids Graphviz dot-in-filename issues)
        rendered_temp_file = dot.render(temp_path, cleanup=True)
        
        # Rename the rendered file to the final output name with suffix
        if rendered_temp_file:
            try:
                shutil.move(rendered_temp_file, output_image)
            except (OSError, IOError) as e:
                print(f"  Warning: Could not rename file to {output_image}: {e}")
                print(f"  File created as: {rendered_temp_file}")
        
        render_elapsed = time.time() - render_time_start
        
        total_elapsed = time.time() - start_time
        print(f"✓ Saved graph visualization: {output_image}")
        print(f"  Bidirectional edges: {bidirectional_rendered}/{edges_added} ({100*bidirectional_rendered/edges_added:.1f}%)")
        print(f"\nTiming summary:")
        print(f"  Format vertices: {vertex_elapsed:.1f}s")
        print(f"  Format edges: {edge_elapsed:.1f}s")
        print(f"  Render: {render_elapsed:.1f}s ({render_elapsed/60:.1f} min)")
        print(f"  Total: {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")
        
        return True
        
    except Exception as e:
        print(f"✗ Error generating graph visualization: {e}")
        import traceback
        traceback.print_exc()
        return False


def save_cluster_visualization_report(graph: ig.Graph, vertex_attrs: Dict, 
                                     vertex_to_cluster_map: Dict, component_cluster_info: Dict,
                                     output_path: str):
    """
    Create an SVG visualization report showing cluster assignments.
    This generates a legend-based report showing vertices grouped by cluster.
    
    Args:
        graph: igraph Graph object
        vertex_attrs: Dictionary of vertex attributes
        vertex_to_cluster_map: {vertex_id: (component_id, cluster_id)}
        component_cluster_info: {component_id: {...clusters...}}
        output_path: Output SVG file path
    """
    print(f"\nGenerating cluster visualization report: {output_path}")
    
    if not vertex_to_cluster_map:
        print("✗ No cluster data available")
        return False
    
    try:
        # Build HTML/SVG report
        html_content = """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Cluster Analysis Report</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        h1 {
            color: #333;
            border-bottom: 2px solid #007bff;
            padding-bottom: 10px;
        }
        h2 {
            color: #555;
            margin-top: 30px;
        }
        .component {
            background-color: white;
            border-left: 4px solid #007bff;
            padding: 15px;
            margin: 15px 0;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .cluster {
            background-color: #f9f9f9;
            border-left: 3px solid #17a2b8;
            padding: 10px;
            margin: 10px 0 10px 20px;
            border-radius: 3px;
        }
        .cluster-header {
            font-weight: bold;
            color: #17a2b8;
            margin-bottom: 5px;
        }
        .vertex-list {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(250px, 1fr));
            gap: 10px;
            margin-top: 8px;
        }
        .vertex {
            background-color: white;
            padding: 8px;
            border-radius: 3px;
            border-left: 3px solid #ddd;
            font-size: 0.9em;
            font-family: monospace;
        }
        .vertex-id {
            color: #666;
            font-size: 0.85em;
        }
        .stats {
            background-color: #e7f3ff;
            padding: 15px;
            border-radius: 5px;
            margin: 15px 0;
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 10px;
            margin-top: 10px;
        }
        .stat-item {
            background-color: white;
            padding: 10px;
            border-radius: 3px;
            text-align: center;
        }
        .stat-value {
            font-size: 1.5em;
            font-weight: bold;
            color: #007bff;
        }
        .stat-label {
            font-size: 0.9em;
            color: #666;
            margin-top: 5px;
        }
    </style>
</head>
<body>
    <h1>Cluster Analysis Report</h1>
"""
        
        # Overall statistics
        num_components = len(component_cluster_info)
        total_clusters = sum(info['num_clusters'] for info in component_cluster_info.values())
        total_vertices = len(vertex_attrs)
        
        html_content += f"""
    <div class="stats">
        <h2>Overall Statistics</h2>
        <div class="stats-grid">
            <div class="stat-item">
                <div class="stat-value">{num_components}</div>
                <div class="stat-label">Connected Components</div>
            </div>
            <div class="stat-item">
                <div class="stat-value">{total_clusters}</div>
                <div class="stat-label">Total Clusters</div>
            </div>
            <div class="stat-item">
                <div class="stat-value">{total_vertices}</div>
                <div class="stat-label">Total Features</div>
            </div>
            <div class="stat-item">
                <div class="stat-value">{total_clusters / num_components:.2f}</div>
                <div class="stat-label">Avg Clusters per Component</div>
            </div>
        </div>
    </div>
"""
        
        # Component-by-component report
        html_content += "<h2>Component Details</h2>\n"
        
        for comp_idx in sorted(component_cluster_info.keys()):
            info = component_cluster_info[comp_idx]
            modularity_str = f"{info['modularity']:.3f}" if info['modularity'] is not None else 'N/A'
            
            html_content += f"""
    <div class="component">
        <h3>Component {comp_idx}</h3>
        <p><strong>Method:</strong> {info['method']} | <strong>Clusters:</strong> {info['num_clusters']} | <strong>Modularity:</strong> {modularity_str}</p>
        <p>{info['details']}</p>
"""
            
            # List clusters in this component
            for cluster_id in sorted(info['clusters'].keys()):
                vertex_ids = info['clusters'][cluster_id]
                
                html_content += f"""
        <div class="cluster">
            <div class="cluster-header">Cluster {cluster_id} ({len(vertex_ids)} vertices)</div>
            <div class="vertex-list">
"""
                
                for vertex_id in sorted(vertex_ids):
                    vertex_info = vertex_attrs[vertex_id]
                    label = vertex_info['label']
                    filename = Path(vertex_info['filename']).stem
                    
                    html_content += f"""
                <div class="vertex">
                    {label}
                    <div class="vertex-id">({filename})</div>
                </div>
"""
                
                html_content += """
            </div>
        </div>
"""
            
            html_content += "    </div>\n"
        
        html_content += """
</body>
</html>
"""
        
        # Save as HTML (can be opened in browser)
        with open(output_path.replace('.svg', '.html'), 'w') as f:
            f.write(html_content)
        
        print(f"✓ Saved cluster visualization report (HTML): {output_path.replace('.svg', '.html')}")
        
        return True
        
    except Exception as e:
        print(f"✗ Error generating cluster visualization: {e}")
        import traceback
        traceback.print_exc()
        return False


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

def visualize_component_on_demand(graph: ig.Graph, vertex_attrs: Dict, components: List, 
                                  component_cluster_info: Dict = None, vertex_to_cluster_map: Dict = None,
                                  output_image_template: str = None,
                                  layout_engine: str = 'dot',
                                  layout_use_weights: bool = False,
                                  clustering_method: str = None,
                                  clustering_weight_mode: str = 'inverse',
                                  deleted_vertices_list: List = None,
                                  filter_method: str = None):
    """
    Interactive visualization mode: compute everything, then await user input to visualize specific components.
    
    Args:
        graph: Full igraph Graph object
        vertex_attrs: Dictionary of vertex attributes
        components: List of connected components (vertex ID lists)
        component_cluster_info: Clustering info per component (optional)
        vertex_to_cluster_map: Vertex to cluster mapping (optional)
        output_image_template: Template for output filenames (e.g., "output_component_{n}.svg")
        layout_engine: Graphviz layout engine to use
        layout_use_weights: Whether to use edge weights in layout
        clustering_method: Clustering method used (for filename annotation)
        clustering_weight_mode: Weight mode used (for filename annotation)
        deleted_vertices_list: List of deleted vertices (optional, for on-demand display)
        filter_method: Filter method used for duplicate removal (for filename annotation)
    """
    if not output_image_template:
        output_image_template = "component_{n}.svg"
    
    print(f"\n{'='*70}")
    print("Visualization On Demand Mode")
    print(f"{'='*70}")
    print(f"\nAvailable components: {len(components)}")
    
    # Debug: Check what's in vertex_attrs for first few vertices
    print("\nDEBUG: Checking first 3 vertex attributes:")
    for i in range(min(3, len(vertex_attrs))):
        if i in vertex_attrs:
            attrs = vertex_attrs[i]
            pep_ident = attrs.get('pep_ident', [])
            print(f"  Vertex {i}: label='{attrs.get('label', 'N/A')}', pep_ident={pep_ident}")
    
    print("\nEnter component numbers to visualize (0-indexed, space-separated).")
    print("Example: '0 5 10' to visualize components 0, 5, and 10")
    print("Enter 'all' to visualize all components, or 'quit'/'exit' to finish.\n")
    
    # Build clustering suffix for filenames
    clustering_suffix = ""
    if clustering_method:
        clustering_suffix = f"_{clustering_method}"
        if component_cluster_info:  # Only add weight mode if clustering was done
            clustering_suffix += f"_{clustering_weight_mode[:3]}"
    
    # Build filter method suffix for filenames
    filter_suffix = ""
    if filter_method:
        # Convert filter method to short suffix
        filter_short = "simpmod" if filter_method == "simplified-modularity" else "mod"
        filter_suffix = f"_{filter_short}"
    
    visualized_count = 0
    
    while True:
        try:
            user_input = input("Components to visualize: ").strip().lower()
            
            if user_input in ['quit', 'exit', 'q']:
                print(f"\n✓ Visualization complete. {visualized_count} component(s) visualized.")
                break
            
            if user_input == 'all':
                component_indices = list(range(len(components)))
            else:
                # Parse space-separated component numbers
                try:
                    component_indices = [int(x) for x in user_input.split()]
                except ValueError:
                    print(f"  ✗ Invalid input. Please enter numbers separated by spaces.")
                    continue
            
            # Validate component indices
            valid_indices = []
            for idx in component_indices:
                if 0 <= idx < len(components):
                    valid_indices.append(idx)
                else:
                    print(f"  ⚠ Component {idx} out of range (0-{len(components)-1})")
            
            if not valid_indices:
                print(f"  ✗ No valid component indices provided.")
                continue
            
            # Visualize each requested component
            for comp_idx in valid_indices:
                component_vertices = components[comp_idx]
                
                # Extract subgraph for this component
                subgraph = graph.induced_subgraph(component_vertices)
                
                # Create vertex attributes for subgraph
                sub_vertex_attrs = {new_id: vertex_attrs[orig_id] 
                                   for new_id, orig_id in enumerate(component_vertices)}
                
                # Create cluster map for subgraph if available
                sub_cluster_map = None
                if vertex_to_cluster_map and component_cluster_info and comp_idx in component_cluster_info:
                    sub_cluster_map = {}
                    for new_id, orig_id in enumerate(component_vertices):
                        if orig_id in vertex_to_cluster_map:
                            sub_cluster_map[new_id] = vertex_to_cluster_map[orig_id]
                
                # Get deleted vertices for this component (map original to subgraph IDs)
                deleted_in_component = set()
                comp_deletions = []
                if deleted_vertices_list:
                    original_deleted = {d['vertex_id'] for d in deleted_vertices_list if d['component_id'] == comp_idx}
                    comp_deletions = [d for d in deleted_vertices_list if d['component_id'] == comp_idx]
                    # Map original vertex IDs to subgraph vertex IDs
                    for new_id, orig_id in enumerate(component_vertices):
                        if orig_id in original_deleted:
                            deleted_in_component.add(new_id)
                
                # Add filtering status suffix (filtered/unfiltered)
                filtering_status_suffix = ""
                if filter_method:
                    # indicate if vertices were actually deleted from this component
                    filtering_status_suffix = "_filtered" if deleted_in_component else "_unfiltered"
                
                # Generate output filename with clustering and filter info
                output_image = output_image_template.format(n=comp_idx)
                output_image = Path(output_image).parent / f"{Path(output_image).stem}_component{comp_idx}{clustering_suffix}{filter_suffix}{filtering_status_suffix}{Path(output_image).suffix}"
                
                # Enhance vertex labels with pep_ident information for on-demand visualization
                pep_ident_count = 0
                empty_ident_count = 0
                shown_samples = 0
                print(f"\n  [DEBUG] Processing {len(sub_vertex_attrs)} vertices for component {comp_idx}...")
                sys.stdout.flush()
                for vertex_id in range(len(sub_vertex_attrs)):
                    if vertex_id in sub_vertex_attrs:
                        attrs = sub_vertex_attrs[vertex_id]
                        original_label = attrs.get('label', str(vertex_id))
                        
                        # Get pep_ident if available
                        pep_ident = attrs.get('pep_ident', [])
                        
                        # Debug: Show first 3 vertices with their pep_ident values
                        if shown_samples < 3:
                            print(f"    Vertex {vertex_id} ({original_label}): pep_ident={pep_ident}")
                            shown_samples += 1
                        sys.stdout.flush()
                        
                        if pep_ident and len(pep_ident) > 0:
                            # Format pep_ident for label (first identification only for brevity)
                            pep_str = str(pep_ident[0]) if isinstance(pep_ident[0], str) else str(pep_ident[0])
                            # Truncate long peptides
                            if len(pep_str) > 15:
                                pep_str = pep_str[:12] + "..."
                            # Append to label in brackets on same line
                            attrs['label'] = f"{original_label} [{pep_str}]"
                            pep_ident_count += 1
                        else:
                            # If no pep_ident, don't add anything
                            attrs['label'] = original_label
                            empty_ident_count += 1
                
                print(f"\n  Visualizing component {comp_idx} ({len(component_vertices)} vertices, {subgraph.ecount()} edges)...")
                print(f"  Vertices with peptide identifications: {pep_ident_count}/{len(sub_vertex_attrs)}")
                if empty_ident_count > 0:
                    print(f"  Vertices without peptide identifications: {empty_ident_count}")
                sys.stdout.flush()
                
                # Print summary of deleted vertices for this component if filtering was applied
                if comp_deletions:
                    print(f"\n  ⚠ DELETED VERTICES SUMMARY (Component {comp_idx}):")
                    print(f"    Total deleted: {len(comp_deletions)}")
                    # Group by cluster
                    deletions_by_cluster = {}
                    for deletion in comp_deletions:
                        cluster_id = deletion['cluster_id']
                        if cluster_id not in deletions_by_cluster:
                            deletions_by_cluster[cluster_id] = []
                        deletions_by_cluster[cluster_id].append(deletion)
                    
                    for cluster_id in sorted(deletions_by_cluster.keys()):
                        cluster_dels = deletions_by_cluster[cluster_id]
                        print(f"    Cluster {cluster_id}: {len(cluster_dels)} deleted")
                        for deletion in sorted(cluster_dels, key=lambda x: x['label']):
                            print(f"      - {deletion['label']} (from {Path(deletion['filename']).stem})")
                            print(f"        Reason: {deletion['reason']}")
                
                # Visualize this component
                success = visualize_graph(
                    subgraph,
                    sub_vertex_attrs,
                    str(output_image),
                    cluster_map=sub_cluster_map if vertex_to_cluster_map else None,
                    component_cluster_info={comp_idx: component_cluster_info[comp_idx]} if component_cluster_info and comp_idx in component_cluster_info else None,
                    layout_engine=layout_engine,
                    use_cluster_subgraphs=False,
                    layout_use_weights=layout_use_weights,
                    use_coordinate_layout=False,
                    vertex_coordinates=None,
                    deleted_vertex_ids=deleted_in_component if deleted_in_component else None
                )
                
                if success:
                    visualized_count += 1
                    print(f"  ✓ Saved: {output_image}")
        
        except KeyboardInterrupt:
            print(f"\n\n✓ Visualization cancelled. {visualized_count} component(s) visualized.")
            break
        except Exception as e:
            print(f"  ✗ Error: {e}")


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
    parser.add_argument("--test-mode", action='store_true',
                       help="Enable test mode: load only 10%% of edges for quick testing")
    parser.add_argument("--fraction", type=float, default=None,
                       help="Fraction of edges to load (0.0 to 1.0). If specified, overrides --test-mode fraction")
    parser.add_argument("--output_graphml", help="Output GraphML file path")
    parser.add_argument("--output_gml", help="Output GML file path")
    parser.add_argument("--output_edgelist", help="Output edge list file path")
    parser.add_argument("--output_analysis", help="Output analysis JSON file path")
    parser.add_argument("--output_composition", help="Output graph composition JSON file path")
    parser.add_argument("--output_image", help="Output image file path (PNG or PDF)")
    parser.add_argument("--output_histogram", help="Output degree distribution histogram file path (PNG)")
    parser.add_argument("--skip-analysis", action='store_true', 
                       help="Skip graph analysis computation and statistics (speeds up execution)")
    parser.add_argument("--skip-visualization", action='store_true', 
                       help="Skip graph visualization rendering")
    parser.add_argument("--layout_engine", type=str, default='dot',
                       choices=['dot', 'sfdp', 'fdp', 'neato', 'twopi', 'circo'],
                       help="Graphviz layout engine for visualization (default: dot)")
    parser.add_argument("--layout_use_weights", action='store_true',
                       help="Use edge distances to influence layout (neato/fdp/sfdp only)")
    parser.add_argument("--use-graphviz-clusters", action='store_true',
                       help="Group clusters into Graphviz subgraphs (default: disabled)")
    parser.add_argument("--use-coordinate-layout", action='store_true',
                       help="Use original feature coordinates as initial positions (neato/fdp only, requires --feature-data-jsons)")
    parser.add_argument("--feature-data-jsons", nargs='+',
                       help="List of feature data JSON files containing x_center, y_center coordinates for coordinate-based layout")
    
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
    parser.add_argument("--mega_complexity_threshold", type=float, default=10000,
                       help="Complexity threshold for using heuristic approach instead of exhaustive search. For components with combo_estimate > threshold, uses clusterwise modularity optimization (fast, greedy) instead of evaluating all combinations. Set to 0 or negative to disable heuristic. (default: 1e15)")
    parser.add_argument("--output-deleted-vertices", type=str, default=None,
                       help="Output path for deleted vertices JSON report (default: deleted_vertices_TIMESTAMP.json in current directory)")
    parser.add_argument("--save-deleted-vertices", action='store_true',
                       help="Attempt to recover deleted vertices by adding them to adjacent clusters if they have edges to those clusters and the cluster doesn't have another vertex from the same file")
    parser.add_argument("--disable-multiprocessing", action='store_true',
                       help="Disable multiprocessing for filtering duplicate vertices - use serial processing instead (useful for debugging or when container CPU allocation is limited)")
    
    # Visualization on demand
    parser.add_argument("--visualization-on-demand", action='store_true',
                       help="Enter interactive mode after computing: visualize specific components on demand instead of full graph")
    parser.add_argument("--component-output-template", type=str, default="component_{n}.svg",
                       help="Template for component visualization filenames (default: component_{n}.svg, extended with method/resolution info for Leiden)")
    
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
    
    # Load edges and cutoff parameters
    print(f"\nLoading edges from: {args.input_pkl}")
    timing_marks['load_edges_start'] = time.time()
    edges_full = load_edges_from_paired_features(args.input_pkl)
    timing_marks['load_edges_end'] = time.time()
    
    # Load cutoff parameters from file (can be overridden by command-line args)
    cutoff_params = load_cutoff_params_from_paired_features(args.input_pkl)
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
    
    # Pre-compute bidirectional map from FULL dataset before sampling
    print("\nPre-computing bidirectional edges from full dataset...")
    timing_marks['bidir_start'] = time.time()
    bidirectional_map = compute_bidirectional_map(edges_full)
    timing_marks['bidir_end'] = time.time()
    
    # Apply test mode or fraction sampling
    edges = edges_full
    if args.fraction is not None:
        # Custom fraction specified
        sample_fraction = args.fraction
        edges_sample = sample_edges_for_testing(edges, sample_fraction)
        print(f"  Sampling {sample_fraction*100:.1f}% of edges ({len(edges_sample)} edges)")
        edges = edges_sample
    elif args.test_mode:
        # Test mode: use 10% default
        sample_fraction = 0.1
        edges_sample = sample_edges_for_testing(edges, sample_fraction)
        print(f"  [TEST MODE] Using {sample_fraction*100:.1f}% of edges ({len(edges_sample)} edges)")
        edges = edges_sample
    
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
                                                              bidirectional_map=bidirectional_map,
                                                              ident_lookup=ident_lookup)
    timing_marks['build_graph_end'] = time.time()
    
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
    optimization_result = None  # Track optimization results
    
    # Load unpaired vertices (needed for optimization scoring if available)
    unpaired_vertices = load_unpaired_vertices_json(args.unpaired_json)
    unpaired_singletons_count = len(unpaired_vertices) if unpaired_vertices else 0
    
    if args.enable_clustering:
        print("\n" + "="*70)
        print("Clustering Analysis")
        print("="*70)
        timing_marks['clustering_start'] = time.time()
        
        # Initialize clustering data (will be populated if optimization is performed)
        best_clustering_data = None
        
        # Auto-select resolution if enabled
        if args.auto_select_resolution and args.clustering_method.lower() == 'leiden':
            # Prepare output path for iteration history
            iterations_output_path = args.output_resolution_iterations
            if not iterations_output_path:
                timestamp = time.strftime("%Y%m%d_%H%M%S")
                iterations_output_path = f"resolution_iterations_{timestamp}.json"
            
            # Run optimization - choose method based on parameter
            if args.resolution_optimization_method.lower() == 'grid-search':
                optimal_resolution, best_metrics, iteration_history, best_clustering_data = grid_search_resolution(
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
                    unpaired_singletons_count=unpaired_singletons_count,
                    initial_seed=args.random_seed if args.random_seed is not None else 42,
                    mega_complexity_threshold=args.mega_complexity_threshold
                )
            else:  # gradient-descent (default)
                optimal_resolution, best_metrics, iteration_history, best_clustering_data = select_optimal_resolution(
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
                    unpaired_singletons_count=unpaired_singletons_count,
                    initial_seed=args.random_seed if args.random_seed is not None else 42,
                    mega_complexity_threshold=args.mega_complexity_threshold
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
                'iterations': iteration_history,
                'iterations_output_path': iterations_output_path
            }
            
            # Use optimized resolution for final clustering
            final_resolution = optimal_resolution
            print(f"\n  ✓ Using optimized resolution: {optimal_resolution:.4f}")
            
            # Use the clustering from the best iteration to avoid re-clustering
            # (ensures consistency between optimization and final output)
            if best_clustering_data:
                vertex_to_cluster_map = best_clustering_data.get('vertex_to_cluster_map', {})
                component_cluster_info = best_clustering_data.get('component_cluster_info', {})
                print(f"  ✓ Using clustering from best iteration (deterministic: no re-clustering)")
                # Calculate initial cluster distribution from the saved clustering
                initial_cluster_sizes_final = []
                for comp_id, comp_info in component_cluster_info.items():
                    for cluster_id, vertex_list in comp_info.get('clusters', {}).items():
                        initial_cluster_sizes_final.append(len(vertex_list))
                initial_cluster_distribution_final = create_cluster_size_distribution(initial_cluster_sizes_final) if initial_cluster_sizes_final else {}
            else:
                # Fallback: re-cluster if clustering data not available
                print(f"  ⚠ Warning: Could not retrieve clustering from optimization, re-clustering...")
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
        else:
            final_resolution = args.resolution_parameter
            if args.auto_select_resolution and args.clustering_method.lower() != 'leiden':
                print("  ⚠ Warning: --auto-select-resolution requires --clustering-method leiden")
                print(f"  Using provided resolution: {args.resolution_parameter}\n")
            
            # Perform final clustering with selected/provided resolution
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
    output_graphml = add_cutoff_suffix_to_filename(args.output_graphml, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    output_gml = add_cutoff_suffix_to_filename(args.output_gml, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    output_edgelist = add_cutoff_suffix_to_filename(args.output_edgelist, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    output_image = add_cutoff_suffix_to_filename(args.output_image, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    output_analysis = add_cutoff_suffix_to_filename(args.output_analysis, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    output_composition = add_cutoff_suffix_to_filename(args.output_composition, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    output_histogram = add_cutoff_suffix_to_filename(args.output_histogram, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)

    # Add layout engine to output image filename for clarity
    if output_image and args.layout_engine:
        image_path = Path(output_image)
        suffix_parts = [args.layout_engine]
        if args.layout_use_weights:
            suffix_parts.append("lwe")
        if args.use_coordinate_layout:
            suffix_parts.append("crd")
        if args.enable_clustering and args.clustering_use_weights:
            suffix_parts.append(f"w{args.clustering_weight_mode[:3]}")
        if args.enable_clustering and args.use_graphviz_clusters:
            suffix_parts.append("gvcl")
        suffix = "_" + "_".join(suffix_parts)
        output_image = str(image_path.parent / f"{image_path.stem}{suffix}{image_path.suffix}")
    
    # Load feature coordinates for coordinate-based layout (if enabled)
    vertex_coordinates = None
    if args.use_coordinate_layout:
        if not args.feature_data_jsons:
            print("\\n⚠ Warning: --use-coordinate-layout specified but no --feature-data-jsons provided")
            print("  Coordinate layout will be disabled")
        else:
            coord_map = load_feature_coordinates_from_jsons(args.feature_data_jsons)
            if coord_map:
                vertex_coordinates = map_coordinates_to_vertices(coord_map, vertex_id_map, vertex_attrs)
            if not vertex_coordinates:
                print("  ⚠ Warning: No valid coordinates found, disabling coordinate layout")
                args.use_coordinate_layout = False
    
    # Save graph
    if output_graphml:
        save_graph(graph, output_graphml, format='graphml')
    
    if output_gml:
        save_graph(graph, output_gml, format='gml')
    
    if output_edgelist:
        save_graph(graph, output_edgelist, format='edgelist')
    
    # Save visualization (if not skipped)
    if not args.skip_visualization:
        # Check if visualization-on-demand mode is enabled
        if args.visualization_on_demand:
            # Apply filtering if enabled and clustering is on
            deleted_vertices_list = []
            filtered_vertex_map = vertex_to_cluster_map
            filtered_component_cluster_info = component_cluster_info
            num_clusters_with_recoveries = 0
            
            if args.enable_clustering and args.filter_duplicate_file_vertices and vertex_to_cluster_map:
                timing_marks['filter_start'] = time.time()
                filtered_vertex_map, filtered_component_cluster_info, deleted_vertices_list = filter_duplicate_file_vertices(
                    graph, vertex_to_cluster_map, component_cluster_info, vertex_attrs, args.filter_method, num_cores, not args.disable_multiprocessing,
                    mega_complexity_threshold=args.mega_complexity_threshold
                )
                timing_marks['filter_end'] = time.time()
            
            # Attempt to recover deleted vertices (add them back to adjacent clusters)
            if args.save_deleted_vertices and deleted_vertices_list:
                filtered_component_cluster_info, num_clusters_with_recoveries = save_deleted_vertices(
                    graph, deleted_vertices_list, filtered_component_cluster_info, vertex_attrs, filtered_vertex_map
                )
            
            # Get connected components for interactive visualization
            components = graph.components()
            
            # Build component output template with clustering info (only add resolution/objective for Leiden)
            component_template = args.component_output_template
            if args.enable_clustering and args.clustering_method.lower() == 'leiden':
                # For Leiden, add resolution parameter and objective function to template
                path_obj = Path(component_template)
                component_template = str(path_obj.parent / f"{path_obj.stem}_res={args.resolution_parameter}_ofunc={args.objective_function}{path_obj.suffix}")
            else:
                # For other methods, just use the template as-is
                path_obj = Path(component_template)
                if args.enable_clustering:
                    component_template = str(path_obj.parent / f"{path_obj.stem}_{args.clustering_method.lower()}{path_obj.suffix}")
            
            # Enter interactive visualization mode
            visualize_component_on_demand(
                graph,
                vertex_attrs,
                components,
                component_cluster_info=filtered_component_cluster_info if args.enable_clustering else None,
                vertex_to_cluster_map=filtered_vertex_map if args.enable_clustering else None,
                output_image_template=component_template,
                layout_engine=args.layout_engine,
                layout_use_weights=args.layout_use_weights,
                clustering_method=args.clustering_method if args.enable_clustering else None,
                clustering_weight_mode=args.clustering_weight_mode if args.enable_clustering else 'inverse',
                deleted_vertices_list=deleted_vertices_list if args.enable_clustering and args.filter_duplicate_file_vertices else None,
                filter_method=args.filter_method if args.filter_duplicate_file_vertices else None
            )
        elif output_image:
            visualize_graph(
                graph,
                vertex_attrs,
                output_image,
                cluster_map=vertex_to_cluster_map if args.enable_clustering else None,
                component_cluster_info=component_cluster_info if args.enable_clustering else None,
                layout_engine=args.layout_engine,
                use_cluster_subgraphs=args.use_graphviz_clusters,
                layout_use_weights=args.layout_use_weights,
                use_coordinate_layout=args.use_coordinate_layout,
                vertex_coordinates=vertex_coordinates
            )
        else:
            print("No --output_image specified and --visualization-on-demand not enabled. Visualization skipped.")
    else:
        print("Skipping visualization (--skip-visualization enabled)")
    
    timing_marks['export_end'] = time.time()
    
    # ===== PRINT TIMING SUMMARY =====
    print(f"\n{'='*70}")
    print("TIMING SUMMARY")
    print(f"{'='*70}")
    
    timing_phases = [
        ('Loading edges', 'load_edges_start', 'load_edges_end'),
        ('Computing bidirectional edges', 'bidir_start', 'bidir_end'),
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
    if output_histogram and analysis is not None:
        export_degree_distribution_histogram(graph, vertex_attrs, output_histogram, graph_composition=graph_composition)
    elif output_histogram and analysis is None:
        print("✗ Skipping histogram (analysis was skipped)")
    
    # Save analysis
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
        if unpaired_vertices:
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
        # Unpaired singletons are added AFTER optimization and should NOT be part of the optimization metrics
        # This ensures data consistency: iteration_history records the true optimization state, uncontaminated
        
        print(f"✓ Added {len(unpaired_vertices) if unpaired_vertices else 0} unpaired singletons to cluster structure")
        print(f"  ✓ resolution_iterations.json remains unchanged (represents pure optimization result)")
        print(f"  ✓ Output files (clusters.json) include unpaired singletons as separate components")
        
        if args.output_clusters_json:
            output_clusters = add_cutoff_suffix_to_filename(args.output_clusters_json, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
            export_clusters_to_json(vertex_to_cluster_map, component_cluster_info, vertex_attrs, output_clusters)
            
            # Also generate visualization report
            output_clusters_viz = output_clusters.replace('.json', '_report.html')
            save_cluster_visualization_report(graph, vertex_attrs, vertex_to_cluster_map, component_cluster_info, output_clusters_viz)
            
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
        
        
        # resolution_iterations.json is NOT modified after optimization
        # It represents the pure optimization result without post-processing
        
        if not args.output_clusters_json:
            print("Clustering performed but --output-clusters-json not specified. Cluster data not saved.")
    
    # Return graph for interactive use
    return graph, vertex_attrs, analysis


if __name__ == "__main__":
    main()

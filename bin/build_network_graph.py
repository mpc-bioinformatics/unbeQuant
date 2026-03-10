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
from pathlib import Path
from typing import Dict, List, Tuple
from pdb import set_trace as bp
import time

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
                       bidirectional_map: Dict = None) -> Tuple[ig.Graph, Dict]:
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
            'label': f"{Path(filename).stem}_{feature_idx}"
        }
    
    print(f"Created {len(vertices)} unique vertices from edges")
    
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


def cluster_graph_independent_components(graph: ig.Graph, vertex_attrs: Dict, 
                                          method: str = 'louvain',
                                          use_weights: bool = False,
                                          weight_mode: str = 'inverse') -> Tuple[Dict, Dict]:
    """
    Perform clustering on each connected component independently.
    
    Args:
        graph: igraph Graph object
        vertex_attrs: Dictionary of vertex attributes
        method: Clustering method ('louvain', 'walktrap', 'label_propagation', 'edge_betweenness', 'leiden', 'spinglass')
    
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
            weights = None
            if use_weights:
                weights = _compute_cluster_edge_weights(subgraph, weight_mode)
                subgraph.es['weight'] = weights

            if method.lower() == 'louvain':
                clustering = subgraph.community_multilevel(weights=weights) if use_weights else subgraph.community_multilevel()
            elif method.lower() == 'leiden':
                clustering = subgraph.community_leiden(weights=weights, objective_function='modularity') if use_weights else subgraph.community_leiden(objective_function='modularity')
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
                    'label': vertex_info['label']
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
                   vertex_coordinates: Dict = None):
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
                                  clustering_weight_mode: str = 'inverse'):
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
    """
    if not output_image_template:
        output_image_template = "component_{n}.svg"
    
    print(f"\n{'='*70}")
    print("Visualization On Demand Mode")
    print(f"{'='*70}")
    print(f"\nAvailable components: {len(components)}")
    print("\nEnter component numbers to visualize (0-indexed, space-separated).")
    print("Example: '0 5 10' to visualize components 0, 5, and 10")
    print("Enter 'all' to visualize all components, or 'quit'/'exit' to finish.\n")
    
    # Build clustering suffix for filenames
    clustering_suffix = ""
    if clustering_method:
        clustering_suffix = f"_{clustering_method}"
        if component_cluster_info:  # Only add weight mode if clustering was done
            clustering_suffix += f"_{clustering_weight_mode[:3]}"
    
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
                
                # Generate output filename with clustering info
                output_image = output_image_template.format(n=comp_idx)
                output_image = Path(output_image).parent / f"{Path(output_image).stem}_component{comp_idx}{clustering_suffix}{Path(output_image).suffix}"
                
                print(f"\n  Visualizing component {comp_idx} ({len(component_vertices)} vertices, {subgraph.ecount()} edges)...")
                
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
                    vertex_coordinates=None
                )
                
                if success:
                    visualized_count += 1
                    print(f"  ✓ Saved: {output_image}")
        
        except KeyboardInterrupt:
            print(f"\n\n✓ Visualization cancelled. {visualized_count} component(s) visualized.")
            break
        except Exception as e:
            print(f"  ✗ Error: {e}")


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
                       choices=['louvain', 'walktrap', 'label_propagation', 'edge_betweenness', 'leiden', 'spinglass'],
                       help="Clustering algorithm to use (default: louvain)")
    parser.add_argument("--clustering-use-weights", action='store_true', default=True,
                       help="Use edge distances as weights in clustering (default: enabled)")
    parser.add_argument("--clustering-no-weights", action='store_false', dest='clustering_use_weights',
                       help="Disable weighted clustering and use unweighted algorithms")
    parser.add_argument("--clustering-weight-mode", type=str, default='inverse',
                       choices=['inverse', 'distance'],
                       help="Weight mode for clustering: inverse (1/(1+distance)) or distance (default: inverse)")
    parser.add_argument("--output-clusters-json", help="Output JSON file with cluster information")
    
    # Visualization on demand
    parser.add_argument("--visualization-on-demand", action='store_true',
                       help="Enter interactive mode after computing: visualize specific components on demand instead of full graph")
    parser.add_argument("--component-output-template", type=str, default="component_{n}.svg",
                       help="Template for component visualization filenames (default: component_{n}.svg)")
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print("Building Network Graph")
    print(f"{'='*70}")
    
    # Load edges and cutoff parameters
    print(f"\nLoading edges from: {args.input_pkl}")
    edges_full = load_edges_from_paired_features(args.input_pkl)
    
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
    bidirectional_map = compute_bidirectional_map(edges_full)
    
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
    graph, vertex_id_map, vertex_attrs = build_network_graph(edges, edge_cutoff=args.edge_cutoff,
                                                              mz_cutoff=args.mz_cutoff, rt_cutoff=args.rt_cutoff,
                                                              bidirectional_map=bidirectional_map)
    
    if graph is None:
        return
    
    # Analyze graph (if not skipped)
    analysis = None
    graph_composition = None
    if not args.skip_analysis:
        print("\nAnalyzing graph...")
        analysis = analyze_graph(graph, vertex_attrs)
        print_analysis(analysis)
        
        # Also compute graph composition statistics
        print("\nAnalyzing graph composition...")
        graph_composition = analyze_graph_composition(graph, vertex_attrs)
        
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
    if args.enable_clustering:
        print("\n" + "="*70)
        print("Clustering Analysis")
        print("="*70)
        vertex_to_cluster_map, component_cluster_info = cluster_graph_independent_components(
            graph,
            vertex_attrs,
            method=args.clustering_method.lower(),
            use_weights=args.clustering_use_weights,
            weight_mode=args.clustering_weight_mode
        )
    else:
        print("\nClustering disabled (use --enable-clustering to enable)")
    
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
            # Get connected components for interactive visualization
            components = graph.components()
            
            # Enter interactive visualization mode
            visualize_component_on_demand(
                graph,
                vertex_attrs,
                components,
                component_cluster_info=component_cluster_info if args.enable_clustering else None,
                vertex_to_cluster_map=vertex_to_cluster_map if args.enable_clustering else None,
                output_image_template=args.component_output_template,
                layout_engine=args.layout_engine,
                layout_use_weights=args.layout_use_weights,
                clustering_method=args.clustering_method if args.enable_clustering else None,
                clustering_weight_mode=args.clustering_weight_mode if args.enable_clustering else 'inverse'
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
    
    # Export degree distribution histogram (if not skipped)
    if output_histogram and analysis is not None:
        export_degree_distribution_histogram(graph, vertex_attrs, output_histogram, graph_composition=graph_composition)
    elif output_histogram and analysis is None:
        print("✗ Skipping histogram (analysis was skipped)")
    
    # Save analysis
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
        if args.output_clusters_json:
            output_clusters = add_cutoff_suffix_to_filename(args.output_clusters_json, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
            export_clusters_to_json(vertex_to_cluster_map, component_cluster_info, vertex_attrs, output_clusters)
            
            # Also generate visualization report
            output_clusters_viz = output_clusters.replace('.json', '_report.html')
            save_cluster_visualization_report(graph, vertex_attrs, vertex_to_cluster_map, component_cluster_info, output_clusters_viz)
        else:
            print("Clustering performed but --output-clusters-json not specified. Cluster data not saved.")
    
    # Return graph for interactive use
    return graph, vertex_attrs, analysis


if __name__ == "__main__":
    main()

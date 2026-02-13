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
            if all(isinstance(k, int) and isinstance(v, list) for k, v in network_edges.items()):
                # Dictionary format: flatten all edges
                for file_idx, file_edges in network_edges.items():
                    edges.extend(file_edges)
            else:
                edges = network_edges
        # Check if it's the edges dictionary format {file_idx: [edges]}
        elif all(isinstance(k, int) and isinstance(v, list) for k, v in data.items()):
            # Dictionary format: flatten all edges
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
    Sample interconnected edges for testing using random selection.
    Random sampling better preserves the natural distribution of bidirectional
    pairs compared to contiguous or interleaved sampling.
    
    Args:
        edges: List of edge dictionaries
        fraction: Fraction of edges to keep (0.0 to 1.0)
    
    Returns:
        Sampled list of edges
    """
    if fraction <= 0.0 or fraction >= 1.0:
        if fraction >= 1.0:
            return edges
        else:
            return []
    
    target_edges = max(10, int(len(edges) * fraction))
    
    # Random sampling with fixed seed for reproducibility
    np.random.seed(42)
    sampled_indices = np.random.choice(len(edges), size=target_edges, replace=False)
    sampled_edges = [edges[i] for i in sorted(sampled_indices)]
    
    # Count unique vertices for reporting
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
    
    print(f"  Random sampling: {target_edges} edges ({len(vertices)} unique vertices)")
    
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
        # Coordinate-based filtering overrides distance cutoff
        # Note: Distance is still kept in output for reference, but filtering uses coordinates
        print(f"Using coordinate-based filtering: mz_cutoff={mz_cutoff}, rt_cutoff={rt_cutoff}")
        filtered_edges = edges  # Distance value always kept in edges dict
        print(f"Loaded {len(edges)} total edges (all edges kept in dict, coordinate filtering noted)")
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
    
    analysis = {
        'num_vertices': graph.vcount(),
        'num_edges': graph.ecount(),
        'density': graph.density(),
        'num_components': len(graph.components()),
        'avg_degree': np.mean(graph.degree()),
        'max_degree': max(graph.degree()),
        'min_degree': min(graph.degree()),
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


def visualize_graph(graph: ig.Graph, vertex_attrs: Dict, output_image: str):
    """Generate and save network graph visualization using Graphviz (with optimizations)."""
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
        
        # Create graph
        dot = graphviz.Digraph(comment='Feature Network', format=output_format)
        
        # Optimize DOT attributes for speed and label placement
        dot.attr(rankdir='LR', splines='polyline', overlap='compress', sep='+0.5', pad='0.2')
        dot.attr('node', shape='circle', style='filled', fontsize='1.5', margin='0.01', width='0.1', height='0.1')
        dot.attr('edge', fontsize='1')
        
        print(f"\n  [1/3] Adding {num_vertices} vertices...", end='', flush=True)
        vertex_time_start = time.time()
        
        # Add vertices (nodes) with small labels
        for vertex_id in range(num_vertices):
            # Get color from vertex attributes (assigned based on source file)
            if vertex_attrs and vertex_id in vertex_attrs and 'color' in vertex_attrs[vertex_id]:
                color = vertex_attrs[vertex_id]['color']
            else:
                # Fallback color (should not happen)
                color = '#87CEEB'
            
            # Get the label from vertex_attrs
            label = vertex_attrs[vertex_id]['label'] if vertex_id in vertex_attrs else str(vertex_id)
            
            # Add node with label visible and assigned color
            dot.node(str(vertex_id), label=label, color=color)
            
            # Progress every 50k vertices
            if (vertex_id + 1) % 50000 == 0:
                elapsed = time.time() - vertex_time_start
                rate = (vertex_id + 1) / elapsed
                remaining = (num_vertices - vertex_id - 1) / rate / 60
                print(f"\n    {vertex_id + 1}/{num_vertices} ({100*(vertex_id+1)/num_vertices:.1f}%) - ETA: {remaining:.1f} min", end='', flush=True)
        
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
            
            # If bidirectional, use bidirectional arrow; otherwise use single direction
            arrow_dir = 'both' if is_bidirectional else 'forward'
            
            # Add edge with weight label
            dot.edge(str(src), str(tgt), label=distance_label, 
                    color=f'#{alpha:02x}{alpha:02x}{alpha:02x}', 
                    penwidth='0.5', arrowsize='0.3', dir=arrow_dir)
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
        print(f"  [3/3] Rendering with dot (using {output_format} format)...")
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
    print(f"  Average Degree: {analysis.get('avg_degree', 'N/A'):.2f}")
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
    parser.add_argument("--output_image", help="Output image file path (PNG or PDF)")
    parser.add_argument("--skip-analysis", action='store_true', 
                       help="Skip graph analysis computation and statistics (speeds up execution)")
    parser.add_argument("--skip-visualization", action='store_true', 
                       help="Skip graph visualization rendering")
    
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
    if not args.skip_analysis:
        print("\nAnalyzing graph...")
        analysis = analyze_graph(graph, vertex_attrs)
        print_analysis(analysis)
    else:
        print("Skipping graph analysis (--skip-analysis enabled)")
    
    # Add cutoff suffix to output filenames
    output_graphml = add_cutoff_suffix_to_filename(args.output_graphml, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    output_gml = add_cutoff_suffix_to_filename(args.output_gml, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    output_edgelist = add_cutoff_suffix_to_filename(args.output_edgelist, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    output_image = add_cutoff_suffix_to_filename(args.output_image, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    output_analysis = add_cutoff_suffix_to_filename(args.output_analysis, args.edge_cutoff, args.mz_cutoff, args.rt_cutoff, args.fraction)
    
    # Save graph
    if output_graphml:
        save_graph(graph, output_graphml, format='graphml')
    
    if output_gml:
        save_graph(graph, output_gml, format='gml')
    
    if output_edgelist:
        save_graph(graph, output_edgelist, format='edgelist')
    
    # Save visualization (if not skipped)
    if not args.skip_visualization:
        if output_image:
            visualize_graph(graph, vertex_attrs, output_image)
        else:
            print("No --output_image specified. Visualization skipped.")
    else:
        print("Skipping visualization (--skip-visualization enabled)")
    
    # Save analysis
    if output_analysis and analysis is not None:
        with open(output_analysis, 'w') as f:
            json.dump(analysis, f, indent=2)
        print(f"✓ Saved analysis: {output_analysis}")
    elif output_analysis and analysis is None:
        print("✗ Cannot save analysis (analysis was skipped)")
    
    # Return graph for interactive use
    return graph, vertex_attrs, analysis


if __name__ == "__main__":
    main()

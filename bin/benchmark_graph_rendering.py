#!/usr/bin/env python3
"""
Benchmark graph rendering time and extrapolate to full graph size.
Tests with progressively larger subgraphs to estimate rendering time.
"""

import pickle
import time
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple
import sys

try:
    import igraph as ig
except ImportError:
    print("igraph not found. Install with: pip install python-igraph")
    exit(1)

try:
    import graphviz
except ImportError:
    print("graphviz not found. Install with: pip install graphviz")
    exit(1)


def load_edges(edges_pkl: str) -> List[Dict]:
    """Load edges from pickle file."""
    with open(edges_pkl, 'rb') as f:
        data = pickle.load(f)
    
    edges = []
    if isinstance(data, dict):
        if all(isinstance(k, int) and isinstance(v, list) for k, v in data.items()):
            for file_idx, file_edges in data.items():
                edges.extend(file_edges)
        else:
            edges = data.get('network_edges', [])
    elif isinstance(data, list):
        edges = data
    
    return edges


def build_small_graph(edges: List[Dict], num_edges: int) -> Tuple[ig.Graph, Dict]:
    """Build a graph from first N edges."""
    filtered_edges = edges[:num_edges]
    
    vertices = set()
    vertex_id_map = {}
    vertex_attrs = {}
    
    for edge in filtered_edges:
        if 'file1' in edge:
            file1_info = edge['file1']
            file2_info = edge['file2']
        else:
            file1_info = edge['current_file']
            file2_info = edge['matched_file']
        
        vertex1_key = (file1_info['filename'], file1_info['feature_idx'])
        vertex2_key = (file2_info['filename'], file2_info['feature_idx'])
        
        vertices.add(vertex1_key)
        vertices.add(vertex2_key)
    
    for idx, vertex_key in enumerate(sorted(vertices)):
        vertex_id_map[vertex_key] = idx
        filename, feature_idx = vertex_key
        vertex_attrs[idx] = {
            'filename': filename,
            'feature_idx': feature_idx,
            'label': f"{Path(filename).stem}_{feature_idx}"
        }
    
    graph = ig.Graph(len(vertices))
    
    for vertex_id, attrs in vertex_attrs.items():
        for attr_name, attr_value in attrs.items():
            graph.vs[vertex_id][attr_name] = attr_value
    
    edge_list = []
    edge_weights = []
    
    for edge in filtered_edges:
        if 'file1' in edge:
            file1_info = edge['file1']
            file2_info = edge['file2']
        else:
            file1_info = edge['current_file']
            file2_info = edge['matched_file']
        
        vertex1_key = (file1_info['filename'], file1_info['feature_idx'])
        vertex2_key = (file2_info['filename'], file2_info['feature_idx'])
        
        vertex1_id = vertex_id_map[vertex1_key]
        vertex2_id = vertex_id_map[vertex2_key]
        
        edge_list.append((vertex1_id, vertex2_id))
        edge_weights.append(edge['distance'])
    
    graph.add_edges(edge_list)
    
    for edge_id, weight in enumerate(edge_weights):
        graph.es[edge_id]['weight'] = weight
        graph.es[edge_id]['distance'] = weight
    
    return graph, vertex_attrs


def render_graph_to_graphviz(graph: ig.Graph, vertex_attrs: Dict) -> float:
    """Build and render a dot graph, return time in seconds."""
    start = time.time()
    
    dot = graphviz.Digraph(comment='Feature Network', format='png')
    dot.attr(rankdir='LR', splines='curved', overlap='false')
    dot.attr('node', shape='circle', style='filled', fontsize='8', width='0.4', height='0.4')
    
    degrees = graph.degree()
    min_degree = min(degrees) if degrees else 1
    max_degree = max(degrees) if degrees else 1
    degree_range = max_degree - min_degree if max_degree > min_degree else 1
    
    # Add vertices
    for vertex_id in range(graph.vcount()):
        degree = degrees[vertex_id]
        size = 0.3 + 1.0 * (degree - min_degree) / degree_range if degree_range > 0 else 0.3
        color = '#87CEEB'
        dot.node(str(vertex_id), label='', width=str(size), height=str(size), fillcolor=color)
    
    # Add edges
    for edge in graph.es:
        src, tgt = edge.source, edge.target
        distance = edge['distance'] if 'distance' in edge.attributes() else 1.0
        width = max(0.2, 2.0 / (1.0 + distance))
        alpha = max(15, int(255 * min(1.0, 1.0 / (1.0 + distance))))
        dot.edge(str(src), str(tgt), penwidth=str(width), color=f'#{alpha:02x}{alpha:02x}{alpha:02x}', arrowsize='0.3')
    
    dot.attr('graph', label=f'Benchmark ({graph.vcount()} vertices, {graph.ecount()} edges)')
    
    # Time the actual render
    render_start = time.time()
    temp_path = '/tmp/benchmark_graph'
    dot.render(temp_path, cleanup=True)
    render_time = time.time() - render_start
    
    return render_time


def main():
    if len(sys.argv) < 2:
        print("Usage: python benchmark_graph_rendering.py <edges.pkl>")
        sys.exit(1)
    
    edges_pkl = sys.argv[1]
    print(f"\nLoading edges from {edges_pkl}...")
    edges = load_edges(edges_pkl)
    total_edges = len(edges)
    
    print(f"Total edges: {total_edges}")
    print(f"\nRunning benchmark with progressively larger graphs...")
    print(f"(This will take several minutes)\n")
    
    # Test sizes - aim for logarithmic spacing
    test_sizes = []
    for size in [100, 500, 1000, 5000, 10000, 50000, 100000]:
        if size < total_edges:
            test_sizes.append(size)
    
    results = []
    
    for idx, num_edges in enumerate(test_sizes, 1):
        print(f"[{idx}/{len(test_sizes)}] Testing with {num_edges} edges...", end='', flush=True)
        
        graph, vertex_attrs = build_small_graph(edges, num_edges)
        num_vertices = graph.vcount()
        
        print(f" ({num_vertices} vertices)", end='', flush=True)
        
        try:
            render_time = render_graph_to_graphviz(graph, vertex_attrs)
            print(f" ✓ {render_time:.1f}s")
            results.append({
                'edges': num_edges,
                'vertices': num_vertices,
                'time': render_time
            })
        except KeyboardInterrupt:
            print(" (interrupted)")
            break
        except Exception as e:
            print(f" ✗ Error: {e}")
            continue
    
    if not results:
        print("\nNo successful benchmarks!")
        return
    
    print(f"\n{'='*70}")
    print("Benchmark Results")
    print(f"{'='*70}")
    print(f"{'Edges':>10} {'Vertices':>10} {'Time (s)':>12} {'Edges/sec':>12}")
    print("-" * 70)
    
    for r in results:
        rate = r['edges'] / r['time']
        print(f"{r['edges']:>10} {r['vertices']:>10} {r['time']:>12.1f} {rate:>12.0f}")
    
    # Extrapolate using polynomial fitting
    print(f"\n{'='*70}")
    print("Extrapolation Analysis")
    print(f"{'='*70}")
    
    edges_arr = np.array([r['edges'] for r in results])
    times_arr = np.array([r['time'] for r in results])
    
    # Try different polynomial fits
    for degree in [1, 2]:
        coeffs = np.polyfit(edges_arr, times_arr, degree)
        poly = np.poly1d(coeffs)
        
        print(f"\nPolynomial fit (degree {degree}):")
        print(f"  Equation: {poly}")
        
        # Calculate R-squared
        y_pred = poly(edges_arr)
        ss_res = np.sum((times_arr - y_pred) ** 2)
        ss_tot = np.sum((times_arr - np.mean(times_arr)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        print(f"  R²: {r_squared:.4f}")
        
        # Extrapolate to full size
        full_extrapolated = poly(total_edges)
        print(f"  Extrapolated time for {total_edges} edges: {full_extrapolated:.1f}s ({full_extrapolated/60:.1f} min)")
        
        # Also show extrapolation range
        min_edges = edges_arr.min()
        max_edges = edges_arr.max()
        print(f"  Valid range: {min_edges:.0f} - {max_edges:.0f} edges")
        print(f"  Extrapolation distance: {total_edges - max_edges:.0f} edges beyond test range")
    
    # Linear extrapolation (often more reliable at extremes)
    print(f"\nLinear extrapolation (based on last 2 points):")
    edge_diff = edges_arr[-1] - edges_arr[-2]
    time_diff = times_arr[-1] - times_arr[-2]
    slope = time_diff / edge_diff
    
    # y = y2 + slope * (x - x2)
    remaining_edges = total_edges - edges_arr[-1]
    estimated_additional_time = remaining_edges * slope
    linear_extrapolated = times_arr[-1] + estimated_additional_time
    
    print(f"  Rate: {slope:.6f} seconds per edge")
    print(f"  Remaining edges: {remaining_edges:.0f}")
    print(f"  Extrapolated time for {total_edges} edges: {linear_extrapolated:.1f}s ({linear_extrapolated/60:.1f} min)")
    
    print(f"\n{'='*70}")
    print("Summary")
    print(f"{'='*70}")
    print(f"Full graph size: {total_edges} edges")
    
    # Use the most conservative estimate
    estimates = [poly(total_edges) for degree in [1, 2]]
    estimates.append(linear_extrapolated)
    avg_estimate = np.mean(estimates)
    max_estimate = np.max(estimates)
    
    print(f"Average estimate: {avg_estimate:.1f}s ({avg_estimate/60:.1f} min)")
    print(f"Conservative (max): {max_estimate:.1f}s ({max_estimate/60:.1f} min)")
    print(f"Fastest (linear): {min(estimates):.1f}s ({min(estimates)/60:.1f} min)")


if __name__ == "__main__":
    main()

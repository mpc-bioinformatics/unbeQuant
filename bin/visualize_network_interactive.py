#!/usr/bin/env python3
"""
Fast visualization alternatives for large network graphs.
Uses igraph layouts + interactive HTML visualization instead of Graphviz.
"""

import pickle
import json
import argparse
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple
import time

try:
    import igraph as ig
except ImportError:
    print("igraph not found. Install with: pip install python-igraph")
    exit(1)

try:
    import pyvis.network as net
    PYVIS_AVAILABLE = True
except ImportError:
    PYVIS_AVAILABLE = False
    print("pyvis not found. Install with: pip install pyvis")


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


def build_network_graph(edges: List[Dict], edge_cutoff: float = float('inf')) -> Tuple[ig.Graph, Dict]:
    """Build an igraph network graph from edges."""
    
    # Filter edges by distance cutoff
    filtered_edges = [e for e in edges if e.get('distance', 0) <= edge_cutoff]
    
    print(f"Loaded {len(edges)} total edges, {len(filtered_edges)} pass distance cutoff {edge_cutoff}")
    
    if not filtered_edges:
        print("No edges pass the cutoff threshold!")
        return None, {}
    
    # Build vertex set
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
    
    print(f"Created {len(vertices)} unique vertices from edges")
    
    # Create graph
    graph = ig.Graph(len(vertices))
    
    for vertex_id, attrs in vertex_attrs.items():
        for attr_name, attr_value in attrs.items():
            graph.vs[vertex_id][attr_name] = attr_value
    
    # Add edges
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
    
    print(f"Created graph with {graph.vcount()} vertices and {graph.ecount()} edges")
    
    return graph, vertex_attrs


def export_to_graphml(graph: ig.Graph, output_path: str, vertex_attrs: Dict):
    """Export graph to GraphML format for Cytoscape/Gephi."""
    print(f"Exporting to GraphML (for Cytoscape/Gephi)...")
    
    # Add vertex attributes
    for vertex_id, attrs in vertex_attrs.items():
        for attr_name, attr_value in attrs.items():
            graph.vs[vertex_id][attr_name] = attr_value
    
    # Add edge attributes
    for edge in graph.es:
        edge['label'] = f"{edge['distance']:.3f}"
    
    graph.write_graphml(output_path)
    print(f"✓ Exported to GraphML: {output_path}")
    print(f"  Open with Cytoscape (www.cytoscape.org) or Gephi (gephi.org)")
    print(f"  These tools have better layout algorithms for large graphs")


def create_interactive_html(graph: ig.Graph, output_path: str, vertex_attrs: Dict):
    """Create interactive HTML visualization using pyvis with ALL vertices and edges."""
    if not PYVIS_AVAILABLE:
        print("✗ pyvis not available. Install with: pip install pyvis")
        return False
    
    print(f"Creating interactive HTML visualization with all {graph.vcount()} vertices and {graph.ecount()} edges...")
    start = time.time()
    
    # Create pyvis network with optimizations for large graphs
    net_vis = net.Network(height='750px', width='100%', directed=False, 
                          notebook=False, cdn_resources='remote')
    
    # Calculate degree-based sizing
    degrees = graph.degree()
    min_degree = min(degrees) if degrees else 1
    max_degree = max(degrees) if degrees else 1
    degree_range = max_degree - min_degree if max_degree > min_degree else 1
    
    # Add all vertices
    print(f"  Adding {graph.vcount()} vertices...")
    vertex_add_start = time.time()
    
    for vertex_id in range(graph.vcount()):
        if vertex_id in vertex_attrs:
            title = f"{vertex_attrs[vertex_id]['label']}\nDegree: {degrees[vertex_id]}"
        else:
            title = f"Vertex {vertex_id}\nDegree: {degrees[vertex_id]}"
        
        # Size and color by degree
        degree = degrees[vertex_id]
        size = 10 + 20 * (degree - min_degree) / degree_range if degree_range > 0 else 10
        color_val = int(100 + 155 * (degree - min_degree) / degree_range) if degree_range > 0 else 150
        color = f'#{color_val:02x}{color_val:02x}255'
        
        net_vis.add_node(vertex_id, label='', title=title, color=color, size=size)
        
        # Progress every 50k vertices
        if (vertex_id + 1) % 50000 == 0:
            elapsed = time.time() - vertex_add_start
            rate = (vertex_id + 1) / elapsed
            remaining = (graph.vcount() - vertex_id - 1) / rate
            print(f"    {vertex_id + 1}/{graph.vcount()} vertices ({100*(vertex_id+1)/graph.vcount():.1f}%) - ETA: {remaining:.1f}s", flush=True)
    
    vertex_elapsed = time.time() - vertex_add_start
    print(f"  ✓ Vertices added ({vertex_elapsed:.1f}s)")
    
    # Add all edges
    print(f"  Adding {graph.ecount()} edges...")
    edge_add_start = time.time()
    
    for edge_idx, edge in enumerate(graph.es):
        src, tgt = edge.source, edge.target
        distance = edge['distance'] if 'distance' in edge.attributes() else 1.0
        
        # Edge width inversely proportional to distance
        width = max(0.5, 2.0 / (1.0 + distance))
        
        # Edge opacity based on distance
        alpha = max(50, int(255 * min(1.0, 0.5 / (1.0 + distance))))
        color = f'rgba(128, 128, 128, {alpha/255:.2f})'
        
        net_vis.add_edge(src, tgt, weight=1/distance, title=f'Distance: {distance:.3f}', 
                        width=width, color=color)
        
        # Progress every 100k edges
        if (edge_idx + 1) % 100000 == 0:
            elapsed = time.time() - edge_add_start
            rate = (edge_idx + 1) / elapsed
            remaining = (graph.ecount() - edge_idx - 1) / rate
            print(f"    {edge_idx + 1}/{graph.ecount()} edges ({100*(edge_idx+1)/graph.ecount():.1f}%) - ETA: {remaining:.1f}s", flush=True)
    
    edge_elapsed = time.time() - edge_add_start
    print(f"  ✓ Edges added ({edge_elapsed:.1f}s)")
    
    # Configure physics - minimal for large graphs to avoid template issues
    print(f"  Configuring visualization...")
    net_vis.show_buttons(filter_=['physics'])
    
    # Save HTML directly
    try:
        print(f"  Writing HTML file...")
        net_vis.write_html(output_path)
        total_elapsed = time.time() - start
        
        print(f"✓ Created interactive HTML: {output_path} ({total_elapsed:.1f}s)")
        file_size_mb = Path(output_path).stat().st_size / (1024*1024)
        print(f"  File size: {file_size_mb:.1f} MB")
        print(f"\n  Interactive features:")
        print(f"    - Drag to pan, scroll to zoom")
        print(f"    - Click nodes/edges for details")
        print(f"    - Use physics buttons to adjust layout")
        print(f"    - Node size and color indicate degree (blue=low, red=high)")
        print(f"    ⚠ Warning: Large graphs may require significant browser memory")
        print(f"       If browser is slow, try --output_graphml for desktop tools")
        
        return True
        
    except Exception as e:
        print(f"✗ Error writing HTML: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Fast visualization for large network graphs using igraph + pyvis"
    )
    parser.add_argument("--input_pkl", required=True, 
                       help="Input pickle file (edges dictionary or paired features)")
    parser.add_argument("--edge_cutoff", type=float, default=float('inf'),
                       help="Maximum distance for edge inclusion")
    parser.add_argument("--output_graphml", help="Output GraphML file path (for Cytoscape/Gephi)")
    parser.add_argument("--output_html", help="Output interactive HTML file path (pyvis, includes all vertices/edges)")
    parser.add_argument("--test-mode", action='store_true',
                       help="Test mode: use only first 5000 edges to validate script works")
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print("Building Network Graph (Fast Interactive Visualization)")
    if args.test_mode:
        print("TEST MODE: Using only first 5000 edges")
    print(f"{'='*70}\n")
    
    # Load and build graph
    print(f"Loading edges from: {args.input_pkl}")
    edges = load_edges(args.input_pkl)
    
    if not edges:
        print("✗ No edges found!")
        return
    
    # Test mode: use only subset of edges
    if args.test_mode:
        edges = edges[:5000]
        print(f"✓ Test mode: using first 5000 edges (total available: {len(edges)})")
    
    print(f"\nBuilding graph with edge_cutoff={args.edge_cutoff}...")
    graph, vertex_attrs = build_network_graph(edges, edge_cutoff=args.edge_cutoff)
    
    if graph is None:
        return
    
    # Export options
    if args.output_graphml:
        export_to_graphml(graph, args.output_graphml, vertex_attrs)
    
    if args.output_html:
        create_interactive_html(graph, args.output_html, vertex_attrs)
    
    if not args.output_graphml and not args.output_html:
        print("\n⚠ No output format specified. Recommended:")
        print(f"  For interactive exploration: --output_html network.html")
        print(f"  For desktop tools (Cytoscape/Gephi): --output_graphml network.graphml")
    
    print(f"\n{'='*70}\n")


if __name__ == "__main__":
    main()

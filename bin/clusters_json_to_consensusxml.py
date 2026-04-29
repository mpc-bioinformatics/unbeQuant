#!/usr/bin/env python3

"""
Convert clusters from JSON format to ConsensusXML format.

This script takes the clustered features from the JSON output and converts them
to OpenMS ConsensusXML format, which can be used with other OpenMS tools.

Each cluster becomes a consensusElement with:
- ID: combination of component_id and cluster_id
- Centroid: average of all vertices in the cluster
- groupedElementList: references to all feature vertices in the cluster
"""

import json
import argparse
import sys
from datetime import datetime
from pathlib import Path
import statistics
from typing import Dict, List, Tuple
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert clusters JSON to ConsensusXML format"
    )
    parser.add_argument(
        "--input_json",
        required=True,
        help="Input JSON file with clusters"
    )
    parser.add_argument(
        "--output_xml",
        required=True,
        help="Output ConsensusXML file"
    )
    
    # Workflow parameters from map_mzml_features.nf
    # Feature pairing parameters
    # Note: Pairing method is now hardcoded to KD-tree optimization
    # Note: best_match_only is now hardcoded to True
    parser.add_argument("--mmf_mz_cutoff", type=float, default=0.2)
    parser.add_argument("--mmf_rt_cutoff", type=float, default=100)
    parser.add_argument("--mmf_edges_cutoff", type=float, default=None)
    parser.add_argument("--mmf_normalize_coordinates", type=bool, default=None)
    parser.add_argument("--mmf_postpair_normalize_coordinates", type=bool, default=False)
    parser.add_argument("--mmf_distance_calc_before_scaling", type=bool, default=False)
    parser.add_argument("--mmf_normalize_edge_distances", type=bool, default=False)
    
    # Heatmap parameters
    parser.add_argument("--mmf_round_up_to", type=int, default=2)
    parser.add_argument("--mmf_log_scale", type=bool, default=True)
    parser.add_argument("--mmf_scale_colors", type=bool, default=True)
    parser.add_argument("--mmf_feature_mode", type=str, default="CoM")
    
    # Network graph parameters
    parser.add_argument("--mmf_build_network_graph", type=bool, default=True)
    parser.add_argument("--mmf_graph_edge_cutoff", type=float, default=None)
    parser.add_argument("--mmf_graph_mz_cutoff", type=float, default=None)
    parser.add_argument("--mmf_graph_rt_cutoff", type=float, default=None)
    parser.add_argument("--mmf_graph_layout_engine", type=str, default="sfdp")
    parser.add_argument("--mmf_graph_layout_use_weights", type=bool, default=True)
    
    # Clustering parameters
    parser.add_argument("--mmf_enable_clustering", type=bool, default=True)
    parser.add_argument("--mmf_clustering_method", type=str, default="leiden")
    parser.add_argument("--mmf_clustering_use_weights", type=bool, default=False)
    parser.add_argument("--mmf_clustering_weight_mode", type=str, default="inverse")
    parser.add_argument("--mmf_clustering_resolution_parameter", type=float, default=0.25)
    parser.add_argument("--mmf_clustering_objective_function", type=str, default="CPM")
    
    # Resolution optimization parameters
    parser.add_argument("--mmf_auto_select_resolution", type=bool, default=True)
    parser.add_argument("--mmf_resolution_optimization_method", type=str, default="grid-search")
    parser.add_argument("--mmf_resolution_min", type=float, default=0.01)
    parser.add_argument("--mmf_resolution_max", type=float, default=0.5)
    parser.add_argument("--mmf_resolution_num_points", type=int, default=5)
    parser.add_argument("--mmf_resolution_metric", type=str, default="combined")
    parser.add_argument("--mmf_resolution_lambda", type=int, default=100)
    parser.add_argument("--mmf_resolution_max_iterations", type=int, default=100)
    
    # Filtering parameters
    parser.add_argument("--mmf_filter_duplicate_file_vertices", type=bool, default=True)
    parser.add_argument("--mmf_filter_method", type=str, default="modularity")
    parser.add_argument("--mmf_save_deleted_vertices", type=bool, default=True)
    
    # Reproducibility
    parser.add_argument("--mmf_random_seed", type=int, default=42)
    
    return parser.parse_args()


def bool_to_int(value) -> str:
    """Convert Python boolean to OpenMS int (1 or 0)."""
    if isinstance(value, bool):
        return "1" if value else "0"
    elif isinstance(value, str):
        return "1" if value.lower() in ['true', '1', 'yes'] else "0"
    else:
        return "1" if value else "0"

def extract_filename_base(filename: str) -> str:
    """Extract base filename without extensions."""
    return filename.split('.')[0]


def generate_consensus_element_id(component_id: int, cluster_id: int) -> str:
    """
    Generate consensus element ID by combining component and cluster IDs.
    
    Uses format: e_CCCCCCCCCCLLLLLLLLLL
    where C=component_id (10 digits), L=cluster_id (10 digits)
    """
    comp_str = str(component_id).zfill(10)
    clust_str = str(cluster_id).zfill(10)
    return f"e_{comp_str}{clust_str}"


def calculate_vertex_distances(vertices: List[Dict]) -> List[float]:
    """
    Calculate distances between all pairs of vertices in a cluster.
    Returns list of pairwise distances.
    """
    distances = []
    
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            # Distance in m/z dimension
            mz_diff = abs(vertices[i]['x_center'] - vertices[j]['x_center'])
            # Distance in RT dimension
            rt_diff = abs(vertices[i]['y_center'] - vertices[j]['y_center'])
            
            # Euclidean distance
            distance = np.sqrt(mz_diff**2 + rt_diff**2)
            distances.append(distance)
    
    return distances


def calculate_cluster_quality(vertices: List[Dict]) -> float:
    """
    Calculate quality metric for cluster.
    
    TODO: This is currently set to average distance between vertices.
    This metric should be assessed and potentially replaced with a more 
    meaningful quality score (e.g., based on modularity, intensity ratio, etc.)
    """
    if len(vertices) < 2:
        return 0.0
    
    distances = calculate_vertex_distances(vertices)
    if not distances:
        return 0.0
    
    average_distance = statistics.mean(distances)
    return average_distance


def build_map_list(data: Dict) -> Tuple[Dict[str, int], Dict[int, str]]:
    """
    Extract unique filenames and build mapping between filename and map_id.
    
    Returns:
        Tuple of (filename_to_id_map, id_to_filename_map)
    """
    filenames = set()
    
    # Collect all unique filenames from all vertices
    for component in data.values():
        for cluster in component.get('clusters', {}).values():
            for vertex in cluster.get('vertices', []):
                filename = vertex.get('filename')
                if filename:
                    filenames.add(filename)
    
    # Create sorted mapping
    filename_list = sorted(list(filenames))
    filename_to_id = {name: idx for idx, name in enumerate(filename_list)}
    id_to_filename = {idx: name for idx, name in enumerate(filename_list)}
    
    return filename_to_id, id_to_filename


def generate_unique_id(filename: str, map_id: int) -> str:
    """Generate a unique ID for each map entry."""
    # Use hash of filename and map_id to create a pseudo-unique ID
    combined = f"{filename}_{map_id}".encode()
    hash_val = hash(combined) & 0x7FFFFFFFFFFFFFFF  # Get positive 64-bit value
    return str(hash_val)


def load_json_clusters(json_file: str) -> Dict:
    """Load cluster data from JSON file."""
    with open(json_file, 'r') as f:
        return json.load(f)


def create_consensus_xml(
    data: Dict,
    filename_to_id: Dict[str, int],
    id_to_filename: Dict[int, str],
    output_file: str,
    input_json_file: str = None,
    workflow_params: Dict = None
):
    """
    Create ConsensusXML content from cluster data.
    
    Args:
        data: Cluster data dictionary
        filename_to_id: Mapping of filenames to indices
        id_to_filename: Mapping of indices to filenames
        output_file: Output XML file path
        input_json_file: Input JSON file path
        workflow_params: Dictionary of workflow parameters to include in header
    """
    
    if workflow_params is None:
        workflow_params = {}
    
    # Collect statistics for parameters
    total_components = len(data)
    total_clusters = sum(len(comp.get('clusters', {})) for comp in data.values())
    total_vertices = sum(
        len(cluster.get('vertices', []))
        for comp in data.values()
        for cluster in comp.get('clusters', {}).values()
    )
    
    # Calculate clustering statistics
    cluster_sizes = []
    for component in data.values():
        for cluster in component.get('clusters', {}).values():
            cluster_sizes.append(len(cluster.get('vertices', [])))
    
    avg_cluster_size = statistics.mean(cluster_sizes) if cluster_sizes else 0.0
    min_cluster_size = min(cluster_sizes) if cluster_sizes else 0
    max_cluster_size = max(cluster_sizes) if cluster_sizes else 0
    
    # Start building XML
    xml_lines = []
    
    # XML header
    xml_lines.append('<?xml version="1.0" encoding="ISO-8859-1"?>')
    xml_lines.append('<?xml-stylesheet type="text/xsl" href="https://www.openms.de/xml-stylesheet/ConsensusXML.xsl" ?>')
    
    # Root element with namespace
    xml_lines.append(
        '<consensusXML version="1.7" id="cm_clustered_features" '
        'experiment_type="label-free" '
        'xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
    )
    
    # Data processing section with parameters
    timestamp = datetime.now().isoformat()
    xml_lines.append('    <dataProcessing completion_time="{}">'.format(timestamp))
    xml_lines.append('        <software name="clusters_json_to_consensusxml" version="1.0" />')
    xml_lines.append('        <processingAction name="Cluster conversion from JSON" />')
    
    # Add input/output file parameters
    if input_json_file:
        xml_lines.append('        <UserParam type="string" name="parameter: in" value="{}"/>'.format(input_json_file))
    xml_lines.append('        <UserParam type="string" name="parameter: out" value="{}"/>'.format(output_file))
    
    # Add workflow parameters - Feature pairing
    xml_lines.append('        <!-- Feature pairing parameters -->')
    xml_lines.append('        <UserParam type="string" name="mmf_pairing_method" value="{}"/>'.format(workflow_params.get('mmf_pairing_method', 'kdtree')))
    xml_lines.append('        <UserParam type="int" name="mmf_best_match_only" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_best_match_only', True))))
    xml_lines.append('        <UserParam type="float" name="mmf_mz_cutoff" value="{}"/>'.format(workflow_params.get('mmf_mz_cutoff', 0.2)))
    xml_lines.append('        <UserParam type="float" name="mmf_rt_cutoff" value="{}"/>'.format(workflow_params.get('mmf_rt_cutoff', 100)))
    
    edges_cutoff = workflow_params.get('mmf_edges_cutoff')
    if edges_cutoff is not None:
        xml_lines.append('        <UserParam type="float" name="mmf_edges_cutoff" value="{}"/>'.format(edges_cutoff))
    
    norm_coords = workflow_params.get('mmf_normalize_coordinates')
    if norm_coords is not None:
        xml_lines.append('        <UserParam type="int" name="mmf_normalize_coordinates" value="{}"/>'.format(bool_to_int(norm_coords)))
    
    xml_lines.append('        <UserParam type="int" name="mmf_postpair_normalize_coordinates" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_postpair_normalize_coordinates', False))))
    xml_lines.append('        <UserParam type="int" name="mmf_distance_calc_before_scaling" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_distance_calc_before_scaling', False))))
    xml_lines.append('        <UserParam type="int" name="mmf_normalize_edge_distances" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_normalize_edge_distances', False))))
    
    # Add workflow parameters - Heatmap/Feature mode
    xml_lines.append('        <!-- Feature mode parameters -->')
    xml_lines.append('        <UserParam type="int" name="mmf_round_up_to" value="{}"/>'.format(workflow_params.get('mmf_round_up_to', 2)))
    xml_lines.append('        <UserParam type="int" name="mmf_log_scale" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_log_scale', True))))
    xml_lines.append('        <UserParam type="int" name="mmf_scale_colors" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_scale_colors', True))))
    xml_lines.append('        <UserParam type="string" name="mmf_feature_mode" value="{}"/>'.format(workflow_params.get('mmf_feature_mode', 'CoM')))
    
    # Add workflow parameters - Network graph
    xml_lines.append('        <!-- Network graph parameters -->')
    xml_lines.append('        <UserParam type="int" name="mmf_build_network_graph" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_build_network_graph', True))))
    
    graph_edge_cutoff = workflow_params.get('mmf_graph_edge_cutoff')
    if graph_edge_cutoff is not None:
        xml_lines.append('        <UserParam type="float" name="mmf_graph_edge_cutoff" value="{}"/>'.format(graph_edge_cutoff))
    
    graph_mz_cutoff = workflow_params.get('mmf_graph_mz_cutoff')
    if graph_mz_cutoff is not None:
        xml_lines.append('        <UserParam type="float" name="mmf_graph_mz_cutoff" value="{}"/>'.format(graph_mz_cutoff))
    
    graph_rt_cutoff = workflow_params.get('mmf_graph_rt_cutoff')
    if graph_rt_cutoff is not None:
        xml_lines.append('        <UserParam type="float" name="mmf_graph_rt_cutoff" value="{}"/>'.format(graph_rt_cutoff))
    
    xml_lines.append('        <UserParam type="string" name="mmf_graph_layout_engine" value="{}"/>'.format(workflow_params.get('mmf_graph_layout_engine', 'sfdp')))
    xml_lines.append('        <UserParam type="int" name="mmf_graph_layout_use_weights" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_graph_layout_use_weights', True))))
    
    # Add workflow parameters - Clustering
    xml_lines.append('        <!-- Clustering parameters -->')
    xml_lines.append('        <UserParam type="int" name="mmf_enable_clustering" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_enable_clustering', True))))
    xml_lines.append('        <UserParam type="string" name="mmf_clustering_method" value="{}"/>'.format(workflow_params.get('mmf_clustering_method', 'leiden')))
    xml_lines.append('        <UserParam type="int" name="mmf_clustering_use_weights" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_clustering_use_weights', False))))
    xml_lines.append('        <UserParam type="string" name="mmf_clustering_weight_mode" value="{}"/>'.format(workflow_params.get('mmf_clustering_weight_mode', 'inverse')))
    xml_lines.append('        <UserParam type="float" name="mmf_clustering_resolution_parameter" value="{}"/>'.format(workflow_params.get('mmf_clustering_resolution_parameter', 0.25)))
    xml_lines.append('        <UserParam type="string" name="mmf_clustering_objective_function" value="{}"/>'.format(workflow_params.get('mmf_clustering_objective_function', 'CPM')))
    
    # Add workflow parameters - Resolution optimization
    xml_lines.append('        <!-- Resolution optimization parameters -->')
    xml_lines.append('        <UserParam type="int" name="mmf_auto_select_resolution" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_auto_select_resolution', True))))
    xml_lines.append('        <UserParam type="string" name="mmf_resolution_optimization_method" value="{}"/>'.format(workflow_params.get('mmf_resolution_optimization_method', 'grid-search')))
    xml_lines.append('        <UserParam type="float" name="mmf_resolution_min" value="{}"/>'.format(workflow_params.get('mmf_resolution_min', 0.01)))
    xml_lines.append('        <UserParam type="float" name="mmf_resolution_max" value="{}"/>'.format(workflow_params.get('mmf_resolution_max', 0.5)))
    xml_lines.append('        <UserParam type="int" name="mmf_resolution_num_points" value="{}"/>'.format(workflow_params.get('mmf_resolution_num_points', 5)))
    xml_lines.append('        <UserParam type="string" name="mmf_resolution_metric" value="{}"/>'.format(workflow_params.get('mmf_resolution_metric', 'combined')))
    xml_lines.append('        <UserParam type="int" name="mmf_resolution_lambda" value="{}"/>'.format(workflow_params.get('mmf_resolution_lambda', 100)))
    xml_lines.append('        <UserParam type="int" name="mmf_resolution_max_iterations" value="{}"/>'.format(workflow_params.get('mmf_resolution_max_iterations', 100)))
    
    # Add workflow parameters - Filtering
    xml_lines.append('        <!-- Filtering parameters -->')
    xml_lines.append('        <UserParam type="int" name="mmf_filter_duplicate_file_vertices" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_filter_duplicate_file_vertices', True))))
    xml_lines.append('        <UserParam type="string" name="mmf_filter_method" value="{}"/>'.format(workflow_params.get('mmf_filter_method', 'modularity')))
    xml_lines.append('        <UserParam type="int" name="mmf_save_deleted_vertices" value="{}"/>'.format(bool_to_int(workflow_params.get('mmf_save_deleted_vertices', True))))
    
    # Add workflow parameters - Reproducibility
    xml_lines.append('        <!-- Reproducibility -->')
    xml_lines.append('        <UserParam type="int" name="mmf_random_seed" value="{}"/>'.format(workflow_params.get('mmf_random_seed', 42)))
    
    # Add clustering statistics parameters
    xml_lines.append('        <!-- Clustering statistics -->')
    xml_lines.append('        <UserParam type="int" name="parameter: num_files" value="{}"/>'.format(len(id_to_filename)))
    xml_lines.append('        <UserParam type="int" name="parameter: num_components" value="{}"/>'.format(total_components))
    xml_lines.append('        <UserParam type="int" name="parameter: num_clusters" value="{}"/>'.format(total_clusters))
    xml_lines.append('        <UserParam type="int" name="parameter: num_vertices" value="{}"/>'.format(total_vertices))
    xml_lines.append('        <UserParam type="float" name="parameter: avg_cluster_size" value="{}"/>'.format(avg_cluster_size))
    xml_lines.append('        <UserParam type="int" name="parameter: min_cluster_size" value="{}"/>'.format(min_cluster_size))
    xml_lines.append('        <UserParam type="int" name="parameter: max_cluster_size" value="{}"/>'.format(max_cluster_size))
    xml_lines.append('        <UserParam type="string" name="parameter: quality_metric" value="average_pairwise_euclidean_distance"/>')
    
    # Add file list
    file_list = [id_to_filename[i] for i in sorted(id_to_filename.keys())]
    xml_lines.append('        <UserParam type="stringList" name="parameter: files" value="[{}]"/>'.format(','.join(file_list)))
    
    xml_lines.append('    </dataProcessing>')
    
    # Map list
    xml_lines.append('    <mapList count="{}">'.format(len(id_to_filename)))
    for map_id in sorted(id_to_filename.keys()):
        filename = id_to_filename[map_id]
        unique_id = generate_unique_id(filename, map_id)
        # For mzML filename, use the base name
        mzml_filename = f"{filename}.mzML"
        xml_lines.append(
            '        <map id="{}" name="{}" unique_id="{}" label="" size="0"/>'.format(
                map_id, mzml_filename, unique_id
            )
        )
    xml_lines.append('    </mapList>')
    
    # Consensus element list
    total_elements = 0
    for component in data.values():
        total_elements += len(component.get('clusters', {}))
    
    xml_lines.append('    <consensusElementList>')
    
    # Process each component and its clusters
    for component_key, component in data.items():
        component_id = component.get('component_id', 0)
        
        for cluster_key, cluster in component.get('clusters', {}).items():
            cluster_id = cluster.get('cluster_id', 0)
            vertices = cluster.get('vertices', [])
            
            if not vertices:
                continue
            
            # Generate unique consensus element ID
            consensus_id = generate_consensus_element_id(component_id, cluster_id)
            
            # Calculate centroid and quality
            avg_mz = statistics.mean([v['x_center'] for v in vertices])
            avg_rt = statistics.mean([v['y_center'] for v in vertices])
            avg_intensity = statistics.mean([v['intensity'] for v in vertices])
            quality = calculate_cluster_quality(vertices)
            
            # Get charge from first vertex (or average if needed)
            charge = vertices[0].get('charge', 0)
            
            # Start consensus element
            xml_lines.append(
                '        <consensusElement id="{}" quality="{}" charge="{}">'.format(
                    consensus_id, quality, charge
                )
            )
            
            # Centroid
            xml_lines.append(
                '            <centroid rt="{}" mz="{}" it="{}"/>'.format(
                    avg_rt, avg_mz, avg_intensity
                )
            )
            
            # Grouped element list (references to individual features)
            xml_lines.append('            <groupedElementList>')
            
            for vertex in vertices:
                filename = vertex.get('filename')
                openms_fid = vertex.get('openms_fid', '')
                
                # Extract the numeric ID from openms_fid (remove "f_" prefix)
                feature_id = openms_fid.replace('f_', '') if openms_fid.startswith('f_') else openms_fid
                
                map_id = filename_to_id.get(filename, -1)
                
                if map_id == -1:
                    print(f"Warning: Could not find map_id for filename {filename}")
                    continue
                
                rt = vertex.get('y_center')
                mz = vertex.get('x_center')
                intensity = vertex.get('intensity')
                vertex_charge = vertex.get('charge', 0)
                
                xml_lines.append(
                    '                <element map="{}" id="{}" rt="{}" mz="{}" it="{}" charge="{}"/>'.format(
                        map_id, feature_id, rt, mz, intensity, vertex_charge
                    )
                )
            
            xml_lines.append('            </groupedElementList>')
            xml_lines.append('        </consensusElement>')
    
    xml_lines.append('    </consensusElementList>')
    xml_lines.append('</consensusXML>')
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(xml_lines))
    
    return total_elements


def main():
    args = parse_args()
    
    # Load JSON data
    print(f"Loading cluster data from {args.input_json}...")
    data = load_json_clusters(args.input_json)
    
    # Build filename mapping
    print("Building filename mapping...")
    filename_to_id, id_to_filename = build_map_list(data)
    
    print(f"  Found {len(id_to_filename)} unique files:")
    for map_id in sorted(id_to_filename.keys()):
        print(f"    Map {map_id}: {id_to_filename[map_id]}")
    
    # Collect workflow parameters
    workflow_params = {
        'mmf_pairing_method': 'kdtree',  # Hardcoded to KD-tree optimization
        'mmf_best_match_only': True,
        'mmf_mz_cutoff': args.mmf_mz_cutoff,
        'mmf_rt_cutoff': args.mmf_rt_cutoff,
        'mmf_edges_cutoff': args.mmf_edges_cutoff,
        'mmf_normalize_coordinates': args.mmf_normalize_coordinates,
        'mmf_postpair_normalize_coordinates': args.mmf_postpair_normalize_coordinates,
        'mmf_distance_calc_before_scaling': args.mmf_distance_calc_before_scaling,
        'mmf_normalize_edge_distances': args.mmf_normalize_edge_distances,
        'mmf_round_up_to': args.mmf_round_up_to,
        'mmf_log_scale': args.mmf_log_scale,
        'mmf_scale_colors': args.mmf_scale_colors,
        'mmf_feature_mode': args.mmf_feature_mode,
        'mmf_build_network_graph': args.mmf_build_network_graph,
        'mmf_graph_edge_cutoff': args.mmf_graph_edge_cutoff,
        'mmf_graph_mz_cutoff': args.mmf_graph_mz_cutoff,
        'mmf_graph_rt_cutoff': args.mmf_graph_rt_cutoff,
        'mmf_graph_layout_engine': args.mmf_graph_layout_engine,
        'mmf_graph_layout_use_weights': args.mmf_graph_layout_use_weights,
        'mmf_enable_clustering': args.mmf_enable_clustering,
        'mmf_clustering_method': args.mmf_clustering_method,
        'mmf_clustering_use_weights': args.mmf_clustering_use_weights,
        'mmf_clustering_weight_mode': args.mmf_clustering_weight_mode,
        'mmf_clustering_resolution_parameter': args.mmf_clustering_resolution_parameter,
        'mmf_clustering_objective_function': args.mmf_clustering_objective_function,
        'mmf_auto_select_resolution': args.mmf_auto_select_resolution,
        'mmf_resolution_optimization_method': args.mmf_resolution_optimization_method,
        'mmf_resolution_min': args.mmf_resolution_min,
        'mmf_resolution_max': args.mmf_resolution_max,
        'mmf_resolution_num_points': args.mmf_resolution_num_points,
        'mmf_resolution_metric': args.mmf_resolution_metric,
        'mmf_resolution_lambda': args.mmf_resolution_lambda,
        'mmf_resolution_max_iterations': args.mmf_resolution_max_iterations,
        'mmf_filter_duplicate_file_vertices': args.mmf_filter_duplicate_file_vertices,
        'mmf_filter_method': args.mmf_filter_method,
        'mmf_save_deleted_vertices': args.mmf_save_deleted_vertices,
        'mmf_random_seed': args.mmf_random_seed,
    }
    
    # Create ConsensusXML
    print(f"Converting clusters to ConsensusXML format...")
    total_elements = create_consensus_xml(data, filename_to_id, id_to_filename, args.output_xml, args.input_json, workflow_params)
    
    print(f"✓ Successfully created ConsensusXML with {total_elements} consensus elements")
    print(f"  Output: {args.output_xml}")


if __name__ == "__main__":
    main()

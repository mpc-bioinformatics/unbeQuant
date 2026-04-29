#!/usr/bin/env python3
"""
Generate summary report of feature pairing results.
"""

import argparse
import json
import os
from typing import Dict, List
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import base64
from io import BytesIO
import pandas as pd


def detect_conflicting_pep_idents(vertices: List[Dict]) -> bool:
    """
    Detect if a cluster has conflicting pep_idents.
    
    A conflict exists if vertices with pep_idents don't share at least one common pep_ident.
    Single vertices or vertices without pep_idents don't have conflicts.
    
    Logic:
    - Collect pep_idents from all vertices (handle both list and string types)
    - Find the intersection of all pep_ident sets
    - If intersection is non-empty: NOT conflicting (all share at least one)
    - If intersection is empty: CONFLICTING (no common pep_ident)
    
    Args:
        vertices: List of vertex dictionaries with 'pep_ident' field
        
    Returns:
        True if conflicts exist, False otherwise
    """
    # Collect pep_ident sets for each vertex that has pep_idents
    pep_ident_sets = []
    
    for vertex in vertices:
        pep_ident = vertex.get('pep_ident')
        if pep_ident:
            if isinstance(pep_ident, list):
                # Convert list to set
                pep_set = set(pep_ident) if pep_ident else set()
            else:
                # Single string pep_ident
                pep_set = {pep_ident}
            
            if pep_set:  # Only add non-empty sets
                pep_ident_sets.append(pep_set)
    
    # If 0 or 1 vertices have pep_idents, no conflict possible
    if len(pep_ident_sets) <= 1:
        return False
    
    # Find intersection of all pep_ident sets (what they all share)
    common_pep_idents = pep_ident_sets[0]
    for pep_set in pep_ident_sets[1:]:
        common_pep_idents = common_pep_idents.intersection(pep_set)
        if not common_pep_idents:
            # Early exit if no common pep_idents found
            break
    
    # If there's any common pep_ident, no conflict
    # If intersection is empty, there's a conflict
    return not common_pep_idents


def detect_conflicting_pep_idents_between_clusters(clusters: Dict) -> bool:
    """
    Detect if clusters within a component share identical pep_idents.
    
    Checks if two or more clusters contain the same pep_ident, which indicates
    duplicate pep_idents across clusters within the component.
    
    Args:
        clusters: Dictionary of cluster_id -> cluster_data (each with 'vertices' list)
        
    Returns:
        True if identical pep_idents found in different clusters, False otherwise
    """
    # Collect pep_ident sets for each cluster
    cluster_pep_idents = {}  # Maps cluster_id -> set of pep_idents in that cluster
    
    for cluster_key, cluster_data in clusters.items():
        vertices = cluster_data.get('vertices', [])
        
        # Collect all unique pep_idents from vertices in this cluster
        pep_idents = set()
        for vertex in vertices:
            pep_ident = vertex.get('pep_ident')
            if pep_ident:
                if isinstance(pep_ident, list):
                    pep_idents.update(pep_ident)
                else:
                    pep_idents.add(pep_ident)
        
        if pep_idents:
            cluster_pep_idents[cluster_key] = pep_idents
    
    # If fewer than 2 clusters have pep_idents, no duplicates possible
    if len(cluster_pep_idents) <= 1:
        return False
    
    # Check if any pep_idents are shared between different clusters
    cluster_ids = list(cluster_pep_idents.keys())
    for i in range(len(cluster_ids)):
        for j in range(i + 1, len(cluster_ids)):
            cluster_i_peps = cluster_pep_idents[cluster_ids[i]]
            cluster_j_peps = cluster_pep_idents[cluster_ids[j]]
            
            # Find intersection - if non-empty, we have shared pep_idents
            shared = cluster_i_peps.intersection(cluster_j_peps)
            if shared:
                return True  # Found identical pep_idents in different clusters
    
    return False  # No shared pep_idents found


def generate_pairing_report(paired_json: str, output_summary: str, output_stats: str, edges_json: str = None, output_html: str = None, graph_composition_json: str = None, clusters_json: str = None, filtered_clusters_json: str = None, deleted_vertices_json: str = None, resolution_iterations_json: str = None, parameters: Dict = None):
    """Generate summary report and statistics from paired features."""
    
    # Initialize early for use throughout function
    conflicting_clusters_list = []
    filtered_cluster_summary = {}
    recovered_cluster_summary = {}
    
    with open(paired_json, 'r') as f:
        data = json.load(f)
    
    # Load graph composition if provided
    graph_composition = None
    if graph_composition_json and os.path.exists(graph_composition_json):
        try:
            with open(graph_composition_json, 'r') as f:
                graph_composition = json.load(f)
            print(f"  Loaded graph composition from {os.path.basename(graph_composition_json)}")
        except Exception as e:
            print(f"  Warning: Could not load graph composition: {e}")
    
    # Load cluster data if provided
    clusters_data = None
    if clusters_json and os.path.exists(clusters_json):
        try:
            with open(clusters_json, 'r') as f:
                clusters_data = json.load(f)
            print(f"  Loaded cluster data from {os.path.basename(clusters_json)}")
        except Exception as e:
            print(f"  Warning: Could not load cluster data: {e}")
    
    # Load filtered cluster data if provided
    filtered_clusters_data = None
    if filtered_clusters_json and os.path.exists(filtered_clusters_json):
        try:
            with open(filtered_clusters_json, 'r') as f:
                filtered_clusters_data = json.load(f)
            print(f"  Loaded filtered cluster data from {os.path.basename(filtered_clusters_json)}")
        except Exception as e:
            print(f"  Warning: Could not load filtered cluster data: {e}")
    
    # Initialize deletion metrics (will be populated from resolution_iterations.json)
    total_deleted = 0
    total_recovered = 0
    vertices_deleted_percentage = 0.0
    filter_method = None
    
    # Load unpaired vertices count for report
    # Load total_unpaired from paired_features_unpaired.json if available
    unpaired_singletons_count = 0
    
    # Get the real path (follow symlinks) of paired_json to find unpaired file
    paired_json_real = os.path.realpath(paired_json)
    paired_dir = os.path.dirname(paired_json_real)
    
    # Try multiple strategies to find unpaired JSON
    unpaired_candidates = [
        os.path.join(paired_dir, 'paired_features_unpaired.json'),
        paired_json.replace('paired_features.json', 'paired_features_unpaired.json'),
        os.path.join(paired_dir, 'paired_features.unpaired.json'),
    ]
    
    unpaired_json_path = None
    for candidate in unpaired_candidates:
        if os.path.exists(candidate):
            unpaired_json_path = candidate
            break
    
    print(f"  [DEBUG] paired_json: {paired_json}")
    print(f"  [DEBUG] paired_json_real: {paired_json_real}")
    print(f"  [DEBUG] paired_dir: {paired_dir}")
    print(f"  [DEBUG] unpaired_json_path: {unpaired_json_path}")
    print(f"  [DEBUG] exists: {os.path.exists(unpaired_json_path) if unpaired_json_path else False}")
    if unpaired_json_path and os.path.exists(unpaired_json_path):
        try:
            with open(unpaired_json_path, 'r') as f:
                unpaired_data = json.load(f)
            # Get total_unpaired from the top-level key
            unpaired_singletons_count = unpaired_data.get('total_unpaired', 0)
            print(f"  Loaded unpaired vertices: {unpaired_singletons_count}")
        except Exception as e:
            print(f"  Note: Could not load unpaired vertices file: {e}")
    
    # Load resolution optimization iterations if provided
    resolution_iterations_data = None
    optimization_result = None
    optimization_was_performed = False
    
    if resolution_iterations_json and os.path.exists(resolution_iterations_json):
        try:
            with open(resolution_iterations_json, 'r') as f:
                resolution_iterations_data = json.load(f)
            
            iterations = resolution_iterations_data.get('iterations', [])
            optimization_result = {
                'iterations': iterations,
                'config': resolution_iterations_data.get('optimization_config', {})
            }
            # Check if optimization actually ran (has iterations)
            optimization_was_performed = len(iterations) > 0
            if optimization_was_performed:
                # Find best iteration by score
                best_iter = max(iterations, key=lambda x: x.get('score', float('-inf')))
                best_res = best_iter.get('resolution', 'Unknown')
                print(f"  Loaded resolution optimization results (optimal resolution: {best_res})")
                
                # Extract deletion metrics from best iteration
                total_deleted = best_iter.get('vertices_deleted', 0)
                total_recovered = best_iter.get('vertices_recovered', 0)
                # Get filter_method from optimization_config
                config = resolution_iterations_data.get('optimization_config', {})
                filter_method = config.get('filter_method')
                if filter_method:
                    print(f"  Deletion metrics (best iteration): {total_deleted} deleted, {total_recovered} recovered")
                    print(f"  Filter method used: {filter_method}")
            else:
                best_res = 'N/A'
                print(f"  Resolution iterations file found but no optimization performed (placeholder file)")
        except Exception as e:
            print(f"  Warning: Could not load resolution iterations: {e}")
    
    # Pre-compute resolution string for use in all HTML sections
    # If optimization was performed, use the optimal resolution
    # Otherwise, use the fixed clustering resolution parameter
    optimal_resolution_str = "N/A"
    if optimization_result and optimization_was_performed:
        iterations = optimization_result.get('iterations', [])
        best_iter = max(iterations, key=lambda x: x.get('score', float('-inf')))
        res_val = best_iter.get('resolution', 'N/A')
        optimal_resolution_str = f"{res_val:.4f}" if isinstance(res_val, (int, float)) else str(res_val)
    elif parameters and parameters.get('mmf_clustering_resolution_parameter'):
        # Use fixed resolution parameter if optimization wasn't performed
        res_val = parameters.get('mmf_clustering_resolution_parameter')
        optimal_resolution_str = f"{res_val:.4f}" if isinstance(res_val, (int, float)) else str(res_val)
    
    # Check if this is an edges file format (has numeric keys with lists of edges)
    is_edges_format = isinstance(data, dict) and all(isinstance(v, list) for v in data.values())
    
    # Extract pairing parameters
    pairing_parameters = {}
    cutoff_params = {}
    
    if is_edges_format:
        # Handle edges format directly
        paired_features = []
        for charge_key, edges_list in data.items():
            paired_features.append({'num_files_matched': 0, 'pep_ident': None})
        pep_ident_stats = {}
        total_features_before_filter = len(paired_features)
        total_features_after_filter = len(paired_features)
        # Use the paired_json as edges_json if edges_json not provided
        if not edges_json:
            edges_json = paired_json
    # Handle both old format (list) and new format (dict with statistics)
    elif isinstance(data, dict) and 'paired_features' in data:
        paired_features = data['paired_features']
        pep_ident_stats = data.get('pep_ident_statistics', {})
        total_features_before_filter = data.get('total_features_before_filtering', data.get('total_features_before_filter', len(paired_features)))
        total_features_after_filter = data.get('total_features_after_filter', sum(len(f.get('individual_features', [])) for f in paired_features))
        pairing_parameters = data.get('pairing_parameters', {})
        cutoff_params = data.get('cutoff_params', {})
    else:
        paired_features = data
        pep_ident_stats = {}
        total_features_before_filter = sum(len(f.get('individual_features', [])) if isinstance(f, dict) else 1 for f in paired_features)
        total_features_after_filter = total_features_before_filter
    
    # Calculate percentage of vertices deleted (based on features after filtering)
    if total_features_after_filter > 0:
        vertices_deleted_percentage = (total_deleted / total_features_after_filter) * 100
    else:
        vertices_deleted_percentage = 0.0
    
    print(f"  Processing {len(paired_features)} paired feature groups...")
    
    # Extract distances from edges JSON if provided
    distances = []
    if edges_json and os.path.exists(edges_json):
        try:
            with open(edges_json, 'r') as f:
                edges_data = json.load(f)
            # Flatten the nested structure to collect all distances
            for charge_key, edges_list in edges_data.items():
                for edge in edges_list:
                    if 'distance' in edge:
                        distances.append(edge['distance'])
            print(f"  Extracted {len(distances)} distances from edges file")
        except Exception as e:
            print(f"  Warning: Could not extract distances from edges file: {e}")
    
    # Calculate statistics
    total_features = len(paired_features)
    files_matched = [f.get('num_files_matched', 0) for f in paired_features if isinstance(f, dict)]
    
    # Extract intra-group distances from new match structure (distances_to_group_features)
    intra_group_distances = []
    for match in paired_features:
        if isinstance(match, dict) and 'individual_features' in match:
            for feature in match['individual_features']:
                distances_dict = feature.get('distances_to_group_features', {})
                intra_group_distances.extend(distances_dict.values())
    
    avg_intra_distance = sum(intra_group_distances) / len(intra_group_distances) if intra_group_distances else 0
    max_intra_distance = max(intra_group_distances) if intra_group_distances else 0
    min_intra_distance = min(intra_group_distances) if intra_group_distances else 0
    
    # Count features matched across different numbers of files
    # and analyze file uniqueness and multiplicity within each group
    file_match_distribution = {}
    file_match_uniqueness = {}  # Track file uniqueness and multiplicity
    
    for f in paired_features:
        if isinstance(f, dict):
            num_files = f.get('num_files_matched', 0)
            file_match_distribution[num_files] = file_match_distribution.get(num_files, 0) + 1
            
            # Analyze file uniqueness and multiplicity for this group
            if num_files not in file_match_uniqueness:
                file_match_uniqueness[num_files] = {
                    'all_unique': 0,  # All files unique (one feature per file)
                    'with_duplicates': 0,  # At least one file has multiple features
                    'total': 0,
                    'file_composition': {},  # file_list -> count
                    'multiplicity_distribution': {},  # e.g., "1x2,1x1" -> count (1 file with 2, 1 with 1)
                    'file_appearance': {}  # filename -> count (how many times each file appears in matches)
                }
            
            file_match_uniqueness[num_files]['total'] += 1
            
            # Count file occurrences in this group
            individual_features = f.get('individual_features', [])
            if individual_features:
                file_counts = {}  # filename -> count
                for feat in individual_features:
                    fname = feat.get('filename')
                    file_counts[fname] = file_counts.get(fname, 0) + 1
                
                # Track file appearance frequency
                for fname in file_counts.keys():
                    if fname not in file_match_uniqueness[num_files]['file_appearance']:
                        file_match_uniqueness[num_files]['file_appearance'][fname] = 0
                    file_match_uniqueness[num_files]['file_appearance'][fname] += 1
                
                # Check if all files are unique (each appears exactly once)
                multiplicity_list = sorted(file_counts.values(), reverse=True)
                is_all_unique = all(count == 1 for count in multiplicity_list)
                
                if is_all_unique:
                    file_match_uniqueness[num_files]['all_unique'] += 1
                else:
                    file_match_uniqueness[num_files]['with_duplicates'] += 1
                
                # Track file composition
                files_in_group = sorted(set(fname for fname in file_counts.keys()))
                file_comp_str = ', '.join(files_in_group)
                if file_comp_str not in file_match_uniqueness[num_files]['file_composition']:
                    file_match_uniqueness[num_files]['file_composition'][file_comp_str] = 0
                file_match_uniqueness[num_files]['file_composition'][file_comp_str] += 1
                
                # Track multiplicity pattern (e.g., "2x2,1x1" for 2 files with 2 features, 1 with 1)
                multiplicity_pattern = ','.join([f"{count}x{multiplicity_list.count(count)}" for count in sorted(set(multiplicity_list), reverse=True)])
                if multiplicity_pattern not in file_match_uniqueness[num_files]['multiplicity_distribution']:
                    file_match_uniqueness[num_files]['multiplicity_distribution'][multiplicity_pattern] = 0
                file_match_uniqueness[num_files]['multiplicity_distribution'][multiplicity_pattern] += 1
    
    # Count features with peptide identifications
    features_with_ident = sum(1 for f in paired_features if isinstance(f, dict) and f.get('pep_ident'))
    
    # Extract file linkage scores from cluster data if available
    file_linkage_scores = None
    if clusters_data and isinstance(clusters_data, dict):
        file_linkage_scores = clusters_data.get('file_linkage_scores')
    
    # Generate HTML report if requested
    if output_html:
        generate_html_report(
            output_html,
            total_features,
            total_features_before_filter,
            total_features_after_filter,
            features_with_ident,
            avg_intra_distance,
            max_intra_distance,
            min_intra_distance,
            file_match_distribution,
            file_match_uniqueness,
            pep_ident_stats,
            distances,
            graph_composition,
            clusters_data,
            filtered_clusters_data,
            total_deleted,
            total_recovered,
            vertices_deleted_percentage,
            filter_method,
            intra_group_distances,
            pairing_parameters,
            cutoff_params,
            optimization_result,
            parameters,
            filtered_clusters_json,
            conflicting_clusters_list,
            filtered_cluster_summary,
            recovered_cluster_summary,
            optimal_resolution_str,
            file_linkage_scores,
            unpaired_singletons_count
        )
    
    # Generate summary text file
    with open(output_summary, 'w') as f:
        f.write("="*70 + "\n")
        f.write("FEATURE PAIRING SUMMARY REPORT\n")
        f.write("="*70 + "\n\n")
        
        # Write pairing parameters at the beginning
        if pairing_parameters:
            f.write("PAIRING PARAMETERS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Pairing Method: {pairing_parameters.get('pairing_method', 'unknown')}\n")
            f.write(f"Best Match Only: True\n")
            f.write(f"Normalize Coordinates: {pairing_parameters.get('normalize_coordinates', True)}\n")
            f.write(f"Distance Calc Before Scaling: {pairing_parameters.get('distance_calc_before_scaling', False)}\n")
            f.write(f"Match Cutoff: {pairing_parameters.get('match_cutoff', 0.0)}\n")
            f.write(f"Input Files: {pairing_parameters.get('input_files', 'unknown')}\n")
            f.write(f"Skip Matchfinder: {pairing_parameters.get('skip_matchfinder', False)}\n\n")
        
        if cutoff_params:
            f.write("FILTERING PARAMETERS\n")
            f.write("-" * 70 + "\n")
            if cutoff_params.get('mz_cutoff') is not None:
                f.write(f"M/Z Cutoff: {cutoff_params.get('mz_cutoff')}\n")
            if cutoff_params.get('rt_cutoff') is not None:
                f.write(f"RT Cutoff: {cutoff_params.get('rt_cutoff')}\n")
            if cutoff_params.get('edges_cutoff') is not None:
                f.write(f"Edges Distance Cutoff: {cutoff_params.get('edges_cutoff')}\n")
            if not any(v is not None for v in cutoff_params.values()):
                f.write("No filtering applied\n")
            f.write("\n")
        
        # Write clustering parameters if available
        if parameters.get('enable_clustering') or parameters.get('clustering_method'):
            f.write("CLUSTERING PARAMETERS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Clustering Enabled: {parameters.get('enable_clustering', False)}\n")
            f.write(f"Clustering Method: {parameters.get('clustering_method', 'N/A')}\n")
            f.write(f"Use Weights: {parameters.get('clustering_use_weights', False)}\n")
            f.write(f"Weight Mode: {parameters.get('clustering_weight_mode', 'N/A')}\n")
            if parameters.get('clustering_resolution_parameter') is not None:
                f.write(f"Resolution Parameter: {parameters.get('clustering_resolution_parameter')}\n")
            if parameters.get('clustering_objective_function'):
                f.write(f"Objective Function: {parameters.get('clustering_objective_function')}\n")
            f.write("\n")
        
        f.write("OVERVIEW\n")
        f.write("-" * 70 + "\n")
        f.write(f"Total paired feature groups: {total_features}\n")
        f.write(f"Total features before filtering: {total_features_before_filter}\n")
        f.write(f"Total features after filtering: {total_features_after_filter}\n")
        f.write(f"Features with peptide identifications: {features_with_ident} ({100*features_with_ident/total_features:.1f}%)\n\n")
        
        f.write("INTRA-GROUP DISTANCE STATISTICS\n")
        f.write("-" * 70 + "\n")
        f.write(f"Average intra-group distance: {avg_intra_distance:.4f}\n")
        f.write(f"Maximum intra-group distance: {max_intra_distance:.4f}\n")
        f.write(f"Minimum intra-group distance: {min_intra_distance:.4f}\n")
        f.write(f"Total intra-group pairwise distances: {len(intra_group_distances)}\n\n")
        
        if distances:
            f.write("DISTANCE STATISTICS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Total paired edges: {len(distances)}\n")
            f.write(f"Average distance: {sum(distances)/len(distances):.4f}\n")
            f.write(f"Median distance: {sorted(distances)[len(distances)//2]:.4f}\n")
            f.write(f"Min distance: {min(distances):.4f}\n")
            f.write(f"Max distance: {max(distances):.4f}\n\n")
        
        # Add graph composition statistics if available
        if graph_composition:
            comp_summary = graph_composition.get('composition_summary', {})
            comp_dist = graph_composition.get('component_composition', {})
            if comp_summary:
                f.write("GRAPH COMPOSITION STATISTICS\n")
                f.write("-" * 70 + "\n")
                f.write(f"Total connected components: {comp_summary.get('total_components')}\n")
                f.write(f"Total features in components: {comp_summary.get('total_features')}\n")
                f.write(f"Average features per component: {comp_summary.get('avg_features_per_component', 0):.2f}\n")
                f.write(f"Average files per component: {comp_summary.get('avg_files_per_component', 0):.2f}\n")
                f.write(f"Max features in any component: {comp_summary.get('max_features_in_component')}\n")
                f.write(f"Max files in any component: {comp_summary.get('max_files_in_component')}\n")
                
                # Show distribution by composition
                if comp_dist:
                    f.write("\nComponent distribution by (features, files):\n")
                    for comp_key in sorted(comp_dist.keys(), key=lambda x: tuple(map(int, x.split('_')))):
                        count = comp_dist[comp_key]
                        parts = comp_key.split('_')
                        f.write(f"  {parts[0]} features from {parts[1]} file(s): {count} component{'s' if count != 1 else ''}\n")
                f.write("\n")
        
        f.write("FILE MATCHING DISTRIBUTION\n")
        f.write("-" * 70 + "\n")
        for num_files in sorted(file_match_distribution.keys()):
            count = file_match_distribution[num_files]
            percentage = 100 * count / total_features
            f.write(f"Matched across {num_files} file(s): {count} ({percentage:.1f}%)\n")
        
        # Add peptide identification statistics if available
        if pep_ident_stats:
            f.write("\nPEPTIDE IDENTIFICATION STATISTICS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Groups with all features sharing at least one common pep_ident: {pep_ident_stats.get('groups_with_common_pep_idents', 0)}\n")
            f.write(f"Groups with all features having mismatching pep_idents (no common ones): {pep_ident_stats.get('groups_with_mismatching_pep_idents', 0)}\n")
            f.write(f"Groups with no pep_idents: {pep_ident_stats.get('groups_with_no_pep_idents', 0)}\n")
            f.write(f"Groups with partial pep_idents (some features have none): {pep_ident_stats.get('groups_with_partial_matching', 0)}\n")
            if pep_ident_stats.get('pep_ident_match_distribution'):
                f.write("\nMatching distribution:\n")
                for key, count in pep_ident_stats['pep_ident_match_distribution'].items():
                    f.write(f"  {key}: {count}\n")
        
        # Add conflicting clusters if any exist
        if conflicting_clusters_list:
            f.write("\nCLUSTERS WITH CONFLICTING PEP-IDENTS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Total clusters with conflicts: {len(conflicting_clusters_list)}\n\n")
            f.write("Component ID | Cluster ID | Total Vertices | Vertices with Pep-Idents\n")
            f.write("-" * 70 + "\n")
            for cluster_info in conflicting_clusters_list:
                component_id = str(cluster_info.get('component_id', 'N/A'))
                cluster_id = str(cluster_info.get('cluster_id', 'N/A'))
                num_vertices = cluster_info.get('num_vertices', 0)
                num_with_peps = cluster_info.get('num_vertices_with_pep_idents', 0)
                f.write(f"{component_id:<14} {cluster_id:<12} {num_vertices:<15} {num_with_peps}\n")
        
        # Add filtered clustering summary if available
        if filtered_cluster_summary:
            f.write("\nFILTERED CLUSTERING SUMMARY\n")
            f.write("-" * 70 + "\n")
            f.write(f"Total clusters: {filtered_cluster_summary.get('total_clusters', 0)}\n")
            f.write(f"Average cluster size: {filtered_cluster_summary.get('avg_cluster_size', 0):.2f}\n")
            f.write(f"Min cluster size: {filtered_cluster_summary.get('min_cluster_size', 0)}\n")
            f.write(f"Max cluster size: {filtered_cluster_summary.get('max_cluster_size', 0)}\n")
            f.write(f"Median cluster size: {filtered_cluster_summary.get('median_cluster_size', 0)}\n")
            f.write(f"Components with clustering: {filtered_cluster_summary.get('components_with_clustering', 0)}\n")
            f.write(f"Single-cluster components: {filtered_cluster_summary.get('single_cluster_components', 0)}\n")
            f.write(f"Total vertices deleted: {filtered_cluster_summary.get('total_vertices_deleted', 0)}\n")
            f.write(f"Clusters with conflicting pep-idents: {filtered_cluster_summary.get('clusters_with_conflicting_pep_idents', 0)}\n")
            f.write(f"Components with cluster pep-ident duplicates: {filtered_cluster_summary.get('components_with_cluster_pep_ident_duplicates', 0)}\n")
            if filtered_cluster_summary.get('avg_modularity') is not None:
                f.write(f"Average modularity: {filtered_cluster_summary.get('avg_modularity'):.3f}\n")
        
        f.write("\n" + "="*70 + "\n")
    
    # Generate statistics CSV file
    with open(output_stats, 'w') as f:
        f.write("metric,value\n")
        f.write(f"total_features_before_filter,{total_features_before_filter}\n")
        f.write(f"total_features_after_filter,{total_features_after_filter}\n")
        f.write(f"total_features,{total_features}\n")
        f.write(f"features_with_identifications,{features_with_ident}\n")
        f.write(f"average_intra_group_distance,{avg_intra_distance:.6f}\n")
        f.write(f"max_intra_group_distance,{max_intra_distance:.6f}\n")
        f.write(f"min_intra_group_distance,{min_intra_distance:.6f}\n")
        f.write(f"total_intra_group_distances,{len(intra_group_distances)}\n")
        
        # Write distance statistics
        if distances:
            f.write(f"total_paired_edges,{len(distances)}\n")
            f.write(f"average_distance,{sum(distances)/len(distances):.6f}\n")
            f.write(f"median_distance,{sorted(distances)[len(distances)//2]:.6f}\n")
            f.write(f"min_distance,{min(distances):.6f}\n")
            f.write(f"max_distance,{max(distances):.6f}\n")
            f.write(f"q1_distance,{sorted(distances)[len(distances)//4]:.6f}\n")
            f.write(f"q3_distance,{sorted(distances)[3*len(distances)//4]:.6f}\n")
        
        for num_files in sorted(file_match_distribution.keys()):
            count = file_match_distribution[num_files]
            f.write(f"matched_across_{num_files}_files,{count}\n")
        
        # Write clustering parameters
        f.write(f"clustering_enabled,{parameters.get('enable_clustering', False)}\n")
        f.write(f"clustering_method,{parameters.get('clustering_method', 'N/A')}\n")
        f.write(f"clustering_use_weights,{parameters.get('clustering_use_weights', False)}\n")
        f.write(f"clustering_weight_mode,{parameters.get('clustering_weight_mode', 'N/A')}\n")
        if parameters.get('clustering_resolution_parameter') is not None:
            f.write(f"clustering_resolution_parameter,{parameters.get('clustering_resolution_parameter')}\n")
        if parameters.get('clustering_objective_function'):
            f.write(f"clustering_objective_function,{parameters.get('clustering_objective_function')}\n")
        
        # Write graph composition statistics
        if graph_composition:
            comp_summary = graph_composition.get('composition_summary', {})
            if comp_summary:
                f.write(f"total_components,{comp_summary.get('total_components')}\n")
                f.write(f"total_features_in_components,{comp_summary.get('total_features')}\n")
                f.write(f"avg_features_per_component,{comp_summary.get('avg_features_per_component', 0):.6f}\n")
                f.write(f"avg_files_per_component,{comp_summary.get('avg_files_per_component', 0):.6f}\n")
                f.write(f"max_features_in_component,{comp_summary.get('max_features_in_component')}\n")
                f.write(f"max_files_in_component,{comp_summary.get('max_files_in_component')}\n")
                
                # Write individual composition distribution
                comp_dist = graph_composition.get('component_composition', {})
                for comp_key, count in comp_dist.items():
                    f.write(f"component_composition_{comp_key},{count}\n")
        
        # Write peptide identification statistics
        if pep_ident_stats:
            f.write(f"groups_with_common_pep_idents,{pep_ident_stats.get('groups_with_common_pep_idents', 0)}\n")
            f.write(f"groups_with_mismatching_pep_idents,{pep_ident_stats.get('groups_with_mismatching_pep_idents', 0)}\n")
            f.write(f"groups_with_no_pep_idents,{pep_ident_stats.get('groups_with_no_pep_idents', 0)}\n")
            f.write(f"groups_with_partial_matching,{pep_ident_stats.get('groups_with_partial_matching', 0)}\n")
            
            for key, count in pep_ident_stats.get('pep_ident_match_distribution', {}).items():
                f.write(f"pep_ident_distribution_{key},{count}\n")
        
        # Write filtered clustering statistics
        if filtered_cluster_summary:
            f.write(f"total_clusters,{filtered_cluster_summary.get('total_clusters', 0)}\n")
            f.write(f"avg_cluster_size,{filtered_cluster_summary.get('avg_cluster_size', 0):.6f}\n")
            f.write(f"min_cluster_size,{filtered_cluster_summary.get('min_cluster_size', 0)}\n")
            f.write(f"max_cluster_size,{filtered_cluster_summary.get('max_cluster_size', 0)}\n")
            f.write(f"median_cluster_size,{filtered_cluster_summary.get('median_cluster_size', 0)}\n")
            f.write(f"components_with_clustering,{filtered_cluster_summary.get('components_with_clustering', 0)}\n")
            f.write(f"single_cluster_components,{filtered_cluster_summary.get('single_cluster_components', 0)}\n")
            f.write(f"total_vertices_deleted,{filtered_cluster_summary.get('total_vertices_deleted', 0)}\n")
            f.write(f"clusters_with_conflicting_pep_idents,{filtered_cluster_summary.get('clusters_with_conflicting_pep_idents', 0)}\n")
            f.write(f"components_with_cluster_pep_ident_duplicates,{filtered_cluster_summary.get('components_with_cluster_pep_ident_duplicates', 0)}\n")
            if filtered_cluster_summary.get('avg_modularity') is not None:
                f.write(f"avg_modularity,{filtered_cluster_summary.get('avg_modularity'):.6f}\n")
    
    print(f"✓ Summary report generated")
    if output_html:
        print(f"✓ HTML report generated")
    print(f"✓ Statistics CSV generated")


def create_distance_histogram_base64(distances: List[float]) -> str:
    """Create a distance distribution histogram and return as base64 encoded image."""
    if not distances or len(distances) == 0:
        return None
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create histogram
    ax.hist(distances, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax.set_xlabel('Distance', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title('Distance Distribution Among Paired Features', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    # Add statistics as text
    mean_dist = sum(distances) / len(distances)
    median_dist = sorted(distances)[len(distances)//2]
    ax.axvline(mean_dist, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_dist:.3f}')
    ax.axvline(median_dist, color='green', linestyle='--', linewidth=2, label=f'Median: {median_dist:.3f}')
    ax.legend(fontsize=10)
    
    # Convert to base64
    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=100, bbox_inches='tight')
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode()
    plt.close(fig)
    
    return image_base64


def create_gridsearch_plots_base64(iterations: List[Dict], total_features_after_filter: int, show_absolute: bool = False, optimal_resolution: float = None) -> str:
    """Create gridsearch optimization plots and return as base64 encoded image.
    
    Generates four subplots:
    1. Average cluster size (multicluster only)
    2. Amount of remaining features (vertices not deleted)
    3. Count of single-vertex clusters
    4. Both pep-ident metrics (dual y-axes):
       - Clusters with conflicting pep-idents (left axis)
       - Components with cluster pep-ident duplicates (right axis)
    
    Args:
        iterations: List of iteration dictionaries with metrics
        total_features_after_filter: Total features for percentage calculations
        show_absolute: If True, show absolute values; if False, show relative percentages
        optimal_resolution: Optional optimal resolution value to mark with a vertical line
    """
    if not iterations or len(iterations) == 0:
        return None
    
    # Sort iterations by resolution in ascending order for consistent visualization
    sorted_iterations = sorted(iterations, key=lambda x: x.get('resolution', 0))
    
    # Extract data from iterations
    resolutions = []
    avg_sizes = []
    vertices_deleted_list = []
    single_vertex_clusters = []
    conflicting_pep_idents = []
    component_pep_duplicates = []
    
    for iteration in sorted_iterations:
        resolutions.append(iteration.get('resolution', 0))
        avg_sizes.append(iteration.get('avg_cluster_size', 0))
        vertices_deleted_list.append(iteration.get('vertices_deleted', 0))
        single_vertex_clusters.append(iteration.get('num_single_vertex_clusters', 0))
        conflicting_pep_idents.append(iteration.get('clusters_with_conflicting_pep_idents', 0))
        component_pep_duplicates.append(iteration.get('components_with_cluster_pep_ident_duplicates', 0))
    
    if not resolutions or not avg_sizes or not vertices_deleted_list:
        return None
    
    # Generate metrics for plotting
    if show_absolute:
        # Use absolute values
        plot_avg_sizes = avg_sizes
        remaining_features = [total_features_after_filter - deleted for deleted in vertices_deleted_list]
        plot_single_clusters = single_vertex_clusters
        plot_conflicting_peps = conflicting_pep_idents
        plot_component_dups = component_pep_duplicates
        
        avg_size_ylabel = 'Average Cluster Size'
        avg_size_title = 'Average Cluster Size (All Clusters, After Filtering & Recovery)'
        remaining_ylabel = 'Remaining Features'
        remaining_title = 'Amount of Remaining Features'
        single_ylabel = 'Single-Vertex Clusters'
        single_title = 'Count of Single-Vertex Clusters'
        pep_title = 'Pep-Ident Metrics (Absolute Values)'
        conflict_ylabel = 'Clusters with Conflicts'
        dup_ylabel = 'Components with Duplicates'
    else:
        # Use relative metrics (percentages)
        max_avg_size = max(avg_sizes) if max(avg_sizes) > 0 else 1
        plot_avg_sizes = [(size / max_avg_size) * 100 for size in avg_sizes]
        remaining_features = [
            100 - ((deleted / total_features_after_filter) * 100) 
            if total_features_after_filter > 0 else 0
            for deleted in vertices_deleted_list
        ]
        
        # Relative single-vertex clusters (percentage of total clusters)
        total_clusters_list = [iteration.get('num_clusters', 1) for iteration in iterations]
        plot_single_clusters = [
            (single / total * 100) if total > 0 else 0 
            for single, total in zip(single_vertex_clusters, total_clusters_list)
        ]
        
        # Relative pep-ident metrics (percentage of max)
        max_conflicting = max(conflicting_pep_idents) if conflicting_pep_idents else 1
        max_duplicates = max(component_pep_duplicates) if component_pep_duplicates else 1
        
        plot_conflicting_peps = [
            (conf / max_conflicting * 100) if max_conflicting > 0 else 0
            for conf in conflicting_pep_idents
        ]
        plot_component_dups = [
            (dup / max_duplicates * 100) if max_duplicates > 0 else 0
            for dup in component_pep_duplicates
        ]
        
        avg_size_ylabel = 'Relative Avg Size (%)'
        avg_size_title = 'Relative Average Cluster Size (All Clusters, After Filtering & Recovery)'
        remaining_ylabel = 'Remaining Features (%)'
        remaining_title = 'Relative Amount of Remaining Features'
        single_ylabel = 'Single-Vertex Clusters (%)'
        single_title = 'Relative Count of Single-Vertex Clusters'
        pep_title = 'Pep-Ident Metrics (Relative %)'
        conflict_ylabel = 'Conflicting Clusters (%)'
        dup_ylabel = 'Duplicate Components (%)'
    
    # Create figure with four subplots
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 16))
    title_suffix = '(Absolute Values)' if show_absolute else '(Relative Values)'
    fig.suptitle(f'Resolution Optimization Grid Search Results {title_suffix}', fontsize=16, fontweight='bold', y=0.995)
    
    # Color scheme
    color1 = '#667eea'
    color2 = '#f093fb'
    color3 = '#50c878'
    color4 = '#ff6b6b'
    color5 = '#ffa94d'
    
    # Plot 1: Average cluster size
    ax1.plot(resolutions, plot_avg_sizes, marker='o', linewidth=2.5, markersize=7, color=color1, label='Avg Size')
    ax1.fill_between(resolutions, plot_avg_sizes, alpha=0.3, color=color1)
    ax1.set_ylabel(avg_size_ylabel, fontsize=12, fontweight='bold')
    ax1.set_title(avg_size_title, fontsize=13, fontweight='bold', pad=15)
    ax1.grid(True, alpha=0.3, linestyle='--')
    if not show_absolute:
        ax1.set_ylim([0, 105])
    
    # Add optimal resolution vertical line
    if optimal_resolution is not None:
        ax1.axvline(x=optimal_resolution, color='red', linestyle=':', linewidth=1.5, alpha=0.7, label=f'Optimal: {optimal_resolution:.4f}')
        ax1.legend(loc='best', fontsize=10)
    
    # Add value labels on points
    for x, y in zip(resolutions, plot_avg_sizes):
        label_fmt = f'{y:.1f}%' if not show_absolute else f'{y:.2f}'
        ax1.annotate(label_fmt, xy=(x, y), xytext=(0, 5), textcoords='offset points', 
                     ha='center', fontsize=9, alpha=0.7)
    
    # Plot 2: Remaining features
    ax2.plot(resolutions, remaining_features, marker='s', linewidth=2.5, markersize=7, color=color2, label='Remaining')
    ax2.fill_between(resolutions, remaining_features, alpha=0.3, color=color2)
    ax2.set_ylabel(remaining_ylabel, fontsize=12, fontweight='bold')
    ax2.set_title(remaining_title, fontsize=13, fontweight='bold', pad=15)
    ax2.grid(True, alpha=0.3, linestyle='--')
    if not show_absolute:
        ax2.set_ylim([0, 105])
    
    # Add optimal resolution vertical line
    if optimal_resolution is not None:
        ax2.axvline(x=optimal_resolution, color='red', linestyle=':', linewidth=1.5, alpha=0.7, label=f'Optimal: {optimal_resolution:.4f}')
        ax2.legend(loc='best', fontsize=10)
    
    # Add value labels on points
    for x, y in zip(resolutions, remaining_features):
        label_fmt = f'{y:.1f}%' if not show_absolute else f'{int(y)}'
        ax2.annotate(label_fmt, xy=(x, y), xytext=(0, 5), textcoords='offset points', 
                     ha='center', fontsize=9, alpha=0.7)
    
    # Plot 3: Single-vertex clusters
    ax3.plot(resolutions, plot_single_clusters, marker='^', linewidth=2.5, markersize=7, color=color3, label='Single-Vertex')
    ax3.fill_between(resolutions, plot_single_clusters, alpha=0.3, color=color3)
    ax3.set_ylabel(single_ylabel, fontsize=12, fontweight='bold')
    ax3.set_title(single_title, fontsize=13, fontweight='bold', pad=15)
    ax3.grid(True, alpha=0.3, linestyle='--')
    if not show_absolute:
        ax3.set_ylim([0, 105])
    
    # Add optimal resolution vertical line
    if optimal_resolution is not None:
        ax3.axvline(x=optimal_resolution, color='red', linestyle=':', linewidth=1.5, alpha=0.7, label=f'Optimal: {optimal_resolution:.4f}')
        ax3.legend(loc='best', fontsize=10)
    
    # Add value labels on points
    for x, y in zip(resolutions, plot_single_clusters):
        label_fmt = f'{y:.1f}%' if not show_absolute else f'{int(y)}'
        ax3.annotate(label_fmt, xy=(x, y), xytext=(0, 5), textcoords='offset points', 
                     ha='center', fontsize=9, alpha=0.7)
    
    # Plot 4: Both pep-ident metrics with dual y-axes
    ax4.set_title(pep_title, fontsize=13, fontweight='bold', pad=15)
    ax4.set_xlabel('Resolution Parameter', fontsize=12, fontweight='bold')
    
    # Left y-axis: Conflicting clusters
    ax4.plot(resolutions, plot_conflicting_peps, marker='D', linewidth=2.5, markersize=7, 
             color=color4, label='Clusters with Conflicts')
    ax4.fill_between(resolutions, plot_conflicting_peps, alpha=0.2, color=color4)
    ax4.set_ylabel(conflict_ylabel, fontsize=11, fontweight='bold', color=color4)
    ax4.tick_params(axis='y', labelcolor=color4)
    if not show_absolute:
        ax4.set_ylim([0, 105])
    
    # Add value labels for left plot
    for x, y in zip(resolutions, plot_conflicting_peps):
        label_fmt = f'{y:.1f}%' if not show_absolute else f'{int(y)}'
        ax4.annotate(label_fmt, xy=(x, y), xytext=(0, -15), textcoords='offset points', 
                     ha='center', fontsize=8, alpha=0.7, color=color4)
    
    # Right y-axis: Component duplicates
    ax4_right = ax4.twinx()
    ax4_right.plot(resolutions, plot_component_dups, marker='v', linewidth=2.5, markersize=7, 
                   color=color5, label='Components with Duplicates')
    ax4_right.fill_between(resolutions, plot_component_dups, alpha=0.2, color=color5)
    ax4_right.set_ylabel(dup_ylabel, fontsize=11, fontweight='bold', color=color5)
    ax4_right.tick_params(axis='y', labelcolor=color5)
    if not show_absolute:
        ax4_right.set_ylim([0, 105])
    
    # Add value labels for right plot
    for x, y in zip(resolutions, plot_component_dups):
        label_fmt = f'{y:.1f}%' if not show_absolute else f'{int(y)}'
        ax4_right.annotate(label_fmt, xy=(x, y), xytext=(0, 5), textcoords='offset points', 
                           ha='center', fontsize=8, alpha=0.7, color=color5)
    
    ax4.grid(True, alpha=0.3, linestyle='--')
    
    # Add optimal resolution vertical line
    if optimal_resolution is not None:
        ax4.axvline(x=optimal_resolution, color='red', linestyle=':', linewidth=1.5, alpha=0.7, label=f'Optimal: {optimal_resolution:.4f}')
    
    # Add combined legend
    lines1, labels1 = ax4.get_legend_handles_labels()
    lines2, labels2 = ax4_right.get_legend_handles_labels()
    ax4.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=10)
    
    # Adjust layout
    plt.tight_layout()
    
    # Convert to base64
    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=100, bbox_inches='tight')
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode()
    plt.close(fig)
    
    return image_base64


def generate_html_report(
    output_html: str,
    total_features: int,
    total_features_before_filter: int,
    total_features_after_filter: int,
    features_with_ident: int,
    avg_intra_distance: float,
    max_intra_distance: float,
    min_intra_distance: float,
    file_match_distribution: Dict,
    file_match_uniqueness: Dict,
    pep_ident_stats: Dict,
    distances: List[float],
    graph_composition: Dict = None,
    clusters_data: Dict = None,
    filtered_clusters_data: Dict = None,
    total_deleted: int = 0,
    total_recovered: int = 0,
    vertices_deleted_percentage: float = 0.0,
    filter_method: str = None,
    intra_group_distances: List[float] = None,
    pairing_parameters: Dict = None,
    cutoff_params: Dict = None,
    optimization_result: Dict = None,
    parameters: Dict = None,
    filtered_clusters_json: str = None,
    conflicting_clusters_list: List[Dict] = None,
    filtered_cluster_summary: Dict = None,
    recovered_cluster_summary: Dict = None,
    optimal_resolution_str: str = "N/A",
    file_linkage_scores: Dict = None,
    unpaired_singletons_count: int = 0
):
    """Generate a beautiful HTML report with visualizations."""
    
    # Ensure summary dictionaries have default values
    if filtered_cluster_summary is None:
        filtered_cluster_summary = {}
    if recovered_cluster_summary is None:
        recovered_cluster_summary = {}
    
    # Create histogram
    histogram_base64 = create_distance_histogram_base64(distances) if distances else None
    
    # Calculate distance statistics
    dist_stats = {}
    if distances:
        dist_stats = {
            'count': len(distances),
            'mean': sum(distances) / len(distances),
            'median': sorted(distances)[len(distances)//2],
            'min': min(distances),
            'max': max(distances),
            'q1': sorted(distances)[len(distances)//4],
            'q3': sorted(distances)[3*len(distances)//4],
        }
    
    # Generate pairing parameters HTML
    params_html = ""
    if pairing_parameters:
        params_html = f"""
        <div class="stats-section">
            <h3>Pairing Parameters</h3>
            <div class="stats-grid">
                <div class="stat-box">
                    <div class="stat-value">{pairing_parameters.get('pairing_method', 'unknown').upper()}</div>
                    <div class="stat-label">Pairing Method</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">YES</div>
                    <div class="stat-label">Best Match Only</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{'YES' if pairing_parameters.get('normalize_coordinates', True) else 'NO'}</div>
                    <div class="stat-label">Normalize Coordinates</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{'YES' if pairing_parameters.get('distance_calc_before_scaling', False) else 'NO'}</div>
                    <div class="stat-label">Distance Before Scaling</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{pairing_parameters.get('match_cutoff', 0.0)}</div>
                    <div class="stat-label">Match Cutoff</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{pairing_parameters.get('input_files', 'unknown')}</div>
                    <div class="stat-label">Input Files</div>
                </div>
            </div>
        </div>
        """
    
    # Generate filtering parameters HTML
    filter_html = ""
    if cutoff_params and any(v is not None for v in cutoff_params.values()):
        filter_items_html = ""
        if cutoff_params.get('mz_cutoff') is not None:
            filter_items_html += f"<div>M/Z Cutoff: <strong>{cutoff_params.get('mz_cutoff')}</strong></div>"
        if cutoff_params.get('rt_cutoff') is not None:
            filter_items_html += f"<div>RT Cutoff: <strong>{cutoff_params.get('rt_cutoff')}</strong></div>"
        if cutoff_params.get('edges_cutoff') is not None:
            filter_items_html += f"<div>Edges Distance Cutoff: <strong>{cutoff_params.get('edges_cutoff')}</strong></div>"
        
        if filter_items_html:
            filter_html = f"""
            <div class="stats-section">
                <h3>Filtering Parameters</h3>
                <div style="background: #f8f9fa; padding: 15px; border-radius: 8px; font-size: 1.05em;">
                    {filter_items_html}
                </div>
            </div>
            """
    
    # Generate file distribution HTML with clickable bars
    file_dist_html = """
    <div id="fileDistChart" class="distribution-section">
    """
    
    # Generate data for JavaScript
    file_dist_data = {}
    for num_files in sorted(file_match_distribution.keys()):
        count = file_match_distribution[num_files]
        percentage = 100 * count / total_features
        
        uniqueness_info = file_match_uniqueness.get(num_files, {})
        all_unique_count = uniqueness_info.get('all_unique', 0)
        with_duplicates_count = uniqueness_info.get('with_duplicates', 0)
        file_compositions = uniqueness_info.get('file_composition', {})
        multiplicity_dist = uniqueness_info.get('multiplicity_distribution', {})
        
        # Create composition breakdown
        composition_breakdown = []
        for file_comp, comp_count in sorted(file_compositions.items(), key=lambda x: -x[1]):
            composition_breakdown.append({
                'composition': file_comp,
                'count': comp_count,
                'percentage': 100 * comp_count / count if count > 0 else 0
            })
        
        # Create multiplicity breakdown
        multiplicity_breakdown = []
        for mult_pattern, mult_count in sorted(multiplicity_dist.items(), key=lambda x: -x[1]):
            multiplicity_breakdown.append({
                'pattern': mult_pattern,
                'count': mult_count,
                'percentage': 100 * mult_count / count if count > 0 else 0
            })
        
        # Create file appearance breakdown
        file_appearance_dist = uniqueness_info.get('file_appearance', {})
        file_appearance_breakdown = []
        for fname, appearance_count in sorted(file_appearance_dist.items(), key=lambda x: -x[1]):
            file_appearance_breakdown.append({
                'filename': fname,
                'count': appearance_count,
                'percentage': 100 * appearance_count / count if count > 0 else 0
            })
        
        file_dist_data[num_files] = {
            'total': count,
            'percentage': percentage,
            'all_unique': all_unique_count,
            'with_duplicates': with_duplicates_count,
            'composition': composition_breakdown,
            'multiplicity': multiplicity_breakdown,
            'file_appearance': file_appearance_breakdown
        }
        
        bar_width = min(95, percentage)
        file_dist_html += f"""
        <div class="distribution-row clickable-bar" onclick="showFileDistributionDetails({num_files})" data-num-files="{num_files}" style="cursor: pointer;">
            <div class="dist-label">{num_files} file(s)</div>
            <div class="dist-bar-container">
                <div class="dist-bar" style="width: {bar_width}%"></div>
            </div>
            <div class="dist-value">{count} ({percentage:.1f}%)</div>
            <div class="dist-click-hint">📊 Click for details</div>
        </div>
        """
    
    file_dist_html += """
    </div>
    <script id="fileDistData" type="application/json">
    """ + json.dumps(file_dist_data) + """
    </script>
    """
    
    # Generate cluster size distribution HTML (if clustering data available)
    cluster_dist_html = ""
    cluster_dist_data = {}
    cluster_summary = {}
    
    # Determine if optimization was performed (has iterations data)
    optimization_was_performed = optimization_result and len(optimization_result.get('iterations', [])) > 0
    
    # If optimization was performed, use best iteration cluster distribution
    # This ensures the report shows metrics for the best resolution
    if optimization_result and optimization_was_performed:
        iterations = optimization_result.get('iterations', [])
        best_iteration = max(iterations, key=lambda x: x.get('score', float('-inf')))
        # Use initial distribution (clustering completed, before filtering)
        size_dist = best_iteration.get('initial_cluster_distribution', {})

        
        # Build cluster_dist_html from distribution
        if size_dist:
            cluster_dist_html = """
    <div id="clusterDistChart" class="distribution-section">
    """
            
            total_clusters = sum(int(v) for v in size_dist.values()) if size_dist else 0
            total_with_unpaired = total_clusters + unpaired_singletons_count
            
            for size in sorted([int(s) for s in size_dist.keys()]):
                size_str = str(size)
                count = int(size_dist[size_str]) if size_str in size_dist else 0
                percentage = 100 * count / total_clusters if total_clusters > 0 else 0
                
                # Special handling for size-1: show just percentage bar
                if size == 1 and unpaired_singletons_count > 0:
                    total_size1 = count + unpaired_singletons_count
                    bar_width = min(95, 100 * total_size1 / total_with_unpaired)
                    
                    cluster_dist_data[size] = {
                        'count': count,
                        'percentage': 100 * count / total_clusters if total_clusters > 0 else 0,
                        'unpaired': unpaired_singletons_count,
                    }
                    
                    cluster_dist_html += f"""
        <div class="distribution-row" style="display: flex; align-items: center; margin: 8px 0; font-size: 0.9em;">
            <div style="min-width: 120px; font-weight: 500; color: #333;">{size} vertex/vertices</div>
            <div style="flex: 1; margin: 0 10px;">
                <div style="background: #e0e0e0; border-radius: 4px; height: 20px; overflow: hidden;">
                    <div style="background: linear-gradient(90deg, #667eea, #764ba2); width: {bar_width}%; height: 100%;"></div>
                </div>
            </div>
            <div style="min-width: 150px; text-align: right; color: #666;">{total_size1} cluster{'s' if total_size1 != 1 else ''} ({100 * total_size1 / total_with_unpaired:.1f}%)</div>
        </div>
        """
                else:
                    bar_width = min(95, 100 * count / total_with_unpaired) if total_with_unpaired > 0 else 0
                    cluster_dist_data[size] = {
                        'count': count,
                        'percentage': percentage,
                    }
                    
                    cluster_dist_html += f"""
        <div class="distribution-row" style="display: flex; align-items: center; margin: 8px 0; font-size: 0.9em;">
            <div style="min-width: 120px; font-weight: 500; color: #333;">{size} vertex/vertices</div>
            <div style="flex: 1; margin: 0 10px;">
                <div style="background: #e0e0e0; border-radius: 4px; height: 20px; overflow: hidden;">
                    <div style="background: linear-gradient(90deg, #667eea, #764ba2); width: {bar_width}%; height: 100%;"></div>
                </div>
            </div>
            <div style="min-width: 150px; text-align: right; color: #666;">{count} cluster{'s' if count != 1 else ''} ({percentage:.1f}%)</div>
        </div>
        """
            
            cluster_dist_html += """
    </div>
    <script id="clusterDistData" type="application/json">
    """ + json.dumps(cluster_dist_data) + """
    </script>
    """
        
        # Build summary from best iteration
        if size_dist and total_clusters > 0:
            cluster_sizes = []
            for size_str, count in size_dist.items():
                cluster_sizes.extend([int(size_str)] * int(count))
            
            # Include unpaired singletons in average cluster size calculation
            total_size_sum = sum(cluster_sizes) + unpaired_singletons_count
            total_cluster_count = len(cluster_sizes) + unpaired_singletons_count
            
            cluster_summary = {
                'total_clusters': total_clusters,
                'avg_cluster_size': total_size_sum / total_cluster_count if total_cluster_count > 0 else 0,
                'min_cluster_size': min(cluster_sizes) if cluster_sizes else 0,
                'max_cluster_size': max(cluster_sizes) if cluster_sizes else 0,
                'median_cluster_size': sorted(cluster_sizes)[len(cluster_sizes)//2] if cluster_sizes else 0,
                'components_with_clustering': best_iteration.get('num_clusters', 0),  # Approximate
                'single_cluster_components': sum(1 for size_str, count in size_dist.items() if int(size_str) == 1),
                'avg_modularity': best_iteration.get('avg_modularity'),
            }
    
    elif clusters_data:
        # Process clusters to build size distribution
        cluster_sizes = []
        component_clusters = {}
        modularity_values = []
        single_cluster_components = 0
        
        for comp_key, comp_data in clusters_data.items():
            if isinstance(comp_data, dict) and 'clusters' in comp_data:
                component_id = comp_data.get('component_id', comp_key)
                method = comp_data.get('method', 'unknown')
                modularity = comp_data.get('modularity')
                clusters = comp_data.get('clusters', {})
                
                # Track single cluster components separately (modularity = 0.0)
                if len(clusters) == 1:
                    single_cluster_components += 1
                elif modularity is not None and modularity != 0.0:
                    modularity_values.append(modularity)
                
                component_clusters[component_id] = {
                    'method': method,
                    'clusters': {},
                    'total_clusters': len(clusters)
                }
                
                for cluster_key, cluster_data in clusters.items():
                    num_vertices = cluster_data.get('num_vertices', 0)
                    vertices = cluster_data.get('vertices', [])
                    
                    # Check for duplicates (same file appearing multiple times in cluster)
                    file_counts = {}
                    for vertex in vertices:
                        fname = vertex.get('filename', 'unknown')
                        file_counts[fname] = file_counts.get(fname, 0) + 1
                    
                    has_duplicates = any(count > 1 for count in file_counts.values())
                    
                    cluster_sizes.append(num_vertices)
                    component_clusters[component_id]['clusters'][cluster_key] = {
                        'size': num_vertices,
                        'has_duplicates': has_duplicates,
                        'file_counts': file_counts,
                        'vertices': vertices
                    }
        
        # Build size distribution
        size_dist = {}
        for size in cluster_sizes:
            size_dist[size] = size_dist.get(size, 0) + 1
        
        # Calculate summary statistics
        if cluster_sizes:
            # Include unpaired singletons in average cluster size calculation
            total_size_sum = sum(cluster_sizes) + unpaired_singletons_count
            total_cluster_count = len(cluster_sizes) + unpaired_singletons_count
            
            cluster_summary = {
                'total_clusters': len(cluster_sizes),
                'avg_cluster_size': total_size_sum / total_cluster_count if total_cluster_count > 0 else 0,
                'min_cluster_size': min(cluster_sizes),
                'max_cluster_size': max(cluster_sizes),
                'median_cluster_size': sorted(cluster_sizes)[len(cluster_sizes)//2],
                'components_with_clustering': len(component_clusters),
                'single_cluster_components': single_cluster_components,
                'avg_modularity': sum(modularity_values) / len(modularity_values) if modularity_values else None,
                'modularity_values_count': len(modularity_values)
            }
        
        # Generate cluster distribution HTML
        if size_dist:
            cluster_dist_html = """
    <div id="clusterDistChart" class="distribution-section">
    """
            
            total_clusters = sum(size_dist.values())+ unpaired_singletons_count if unpaired_singletons_count > 0 else sum(size_dist.values())
            
            for size in sorted(size_dist.keys()):
                count = size_dist[size] 
                percentage = 100 * count / total_clusters
                bar_width = min(95, percentage)
                
                # Calculate statistics for clusters of this size
                clusters_with_duplicates = 0
                clusters_without_duplicates = 0
                total_files_in_clusters = 0
                components_with_this_size = set()
                
                for comp_id, comp_info in component_clusters.items():
                    for cluster_key, cluster_info in comp_info['clusters'].items():
                        if cluster_info['size'] == size:
                            components_with_this_size.add(comp_id)
                            if cluster_info['has_duplicates']:
                                clusters_with_duplicates += 1
                            else:
                                clusters_without_duplicates += 1
                            total_files_in_clusters += len(cluster_info['file_counts'])
                
                avg_files_per_cluster = total_files_in_clusters / count if count > 0 else 0
                
                cluster_dist_data[size] = {
                    'count': count,
                    'percentage': percentage,
                    'clusters_with_duplicates': clusters_with_duplicates,
                    'clusters_without_duplicates': clusters_without_duplicates,
                    'avg_files_per_cluster': avg_files_per_cluster,
                    'components_represented': len(components_with_this_size)
                }
                
                cluster_dist_html += f"""
        <div class="distribution-row clickable-bar" onclick="showClusterDistributionDetails({size})" data-cluster-size="{size}" style="cursor: pointer;">
            <div class="dist-label">{size} vertex/vertices</div>
            <div class="dist-bar-container">
                <div class="dist-bar" style="width: {bar_width}%"></div>
            </div>
            <div class="dist-value">{count} cluster{'s' if count != 1 else ''} ({percentage:.1f}%)</div>
            <div class="dist-click-hint">📊 Click for details</div>
        </div>
        """
            
            cluster_dist_html += """
    </div>
    <script id="clusterDistData" type="application/json">
    """ + json.dumps(cluster_dist_data) + """
    </script>
    """
    
    # Generate filtered cluster size distribution HTML (if filtered clustering data available)
    filtered_cluster_dist_html = ""
    filtered_cluster_dist_data = {}
    deletion_summary_html = ""
    filtered_size_dist = {}  # Initialize to empty dict
    total_deleted = 0  # Initialize to 0
    
    # If optimization was performed, use best iteration filtered distribution
    if optimization_result and optimization_was_performed:
        iterations = optimization_result.get('iterations', [])
        best_iteration = max(iterations, key=lambda x: x.get('score', float('-inf')))
        # Use filtered distribution (after filtering, before recovery)
        filtered_size_dist = best_iteration.get('filtered_cluster_distribution', {})
        total_deleted = best_iteration.get('vertices_deleted', 0)
        
        # Build filtered_cluster_dist_html from distribution (this is after filtering but before recovery)
        if filtered_size_dist:
            # Add a note showing deleted vertices
            deleted_note = ""
            if total_deleted > 0:
                deleted_note = f"""
    <div style="background: #fff3cd; padding: 12px; margin-bottom: 15px; border-radius: 4px; color: #856404; font-size: 0.95em; border-left: 4px solid #ffc107;">
        <strong>⚠ Note:</strong> This distribution shows clusters <strong>after filtering</strong> (duplicate file vertices removed).
        <strong style="color: #d9534f;">{total_deleted} vertices were deleted</strong> during filtering due to modularity optimization.
        These vertices will be shown below as <strong>recovered singletons</strong> after being re-added.
    </div>"""
            
            filtered_cluster_dist_html = f"""
    <div id="filteredClusterDistChart" class="distribution-section">
    {deleted_note}
    """
            
            total_filtered_clusters = sum(int(v) for v in filtered_size_dist.values()) if filtered_size_dist else 0
            total_filtered_with_unpaired = total_filtered_clusters + unpaired_singletons_count
            
            for size in sorted([int(s) for s in filtered_size_dist.keys()]):
                size_str = str(size)
                count = int(filtered_size_dist[size_str]) if size_str in filtered_size_dist else 0
                percentage = 100 * count / total_filtered_clusters if total_filtered_clusters > 0 else 0
                
                # Special handling for size-1: show just percentage bar
                if size == 1 and (total_deleted > 0 or unpaired_singletons_count > 0):
                    total_size1 = count + total_deleted + unpaired_singletons_count
                    bar_width = min(95, 100 * total_size1 / total_filtered_with_unpaired)
                    
                    filtered_cluster_dist_data[size] = {
                        'count': count,
                        'percentage': 100 * count / total_filtered_clusters if total_filtered_clusters > 0 else 0,
                        'deleted': total_deleted,
                        'unpaired': unpaired_singletons_count,
                    }
                    
                    filtered_cluster_dist_html += f"""
        <div class="distribution-row" style="display: flex; align-items: center; margin: 8px 0; font-size: 0.9em;">
            <div style="min-width: 120px; font-weight: 500; color: #333;">{size} vertex/vertices</div>
            <div style="flex: 1; margin: 0 10px;">
                <div style="background: #e0e0e0; border-radius: 4px; height: 20px; overflow: hidden;">
                    <div style="background: linear-gradient(90deg, #667eea, #764ba2); width: {bar_width}%; height: 100%;"></div>
                </div>
            </div>
            <div style="min-width: 150px; text-align: right; color: #666;">{total_size1} cluster{'s' if total_size1 != 1 else ''} ({100 * total_size1 / total_filtered_with_unpaired:.1f}%)</div>
        </div>
        """
                else:
                    bar_width = min(95, 100 * count / total_filtered_with_unpaired) if total_filtered_with_unpaired > 0 else 0
                    
                    filtered_cluster_dist_data[size] = {
                        'count': count,
                        'percentage': percentage,
                    }
                    
                    filtered_cluster_dist_html += f"""
        <div class="distribution-row" style="display: flex; align-items: center; margin: 8px 0; font-size: 0.9em;">
            <div style="min-width: 120px; font-weight: 500; color: #333;">{size} vertex/vertices</div>
            <div style="flex: 1; margin: 0 10px;">
                <div style="background: #e0e0e0; border-radius: 4px; height: 20px; overflow: hidden;">
                    <div style="background: linear-gradient(90deg, #667eea, #764ba2); width: {bar_width}%; height: 100%;"></div>
                </div>
            </div>
            <div style="min-width: 150px; text-align: right; color: #666;">{count} cluster{'s' if count != 1 else ''} ({percentage:.1f}%)</div>
        </div>
        """
            
            filtered_cluster_dist_html += """
    </div>
    <script id="filteredClusterDistData" type="application/json">
    """ + json.dumps(filtered_cluster_dist_data) + """
    </script>
    """
        
        # Build filtered summary from best iteration
        if filtered_size_dist and total_filtered_clusters > 0:
            filtered_cluster_sizes = []
            for size_str, count in filtered_size_dist.items():
                filtered_cluster_sizes.extend([int(size_str)] * int(count))
            
            # Include unpaired singletons in average cluster size calculation
            total_size_sum = sum(filtered_cluster_sizes) + unpaired_singletons_count
            total_cluster_count = len(filtered_cluster_sizes) + unpaired_singletons_count
            
            filtered_cluster_summary = {
                'total_clusters': total_filtered_clusters,
                'total_vertices_deleted': total_deleted,
                'avg_cluster_size': total_size_sum / total_cluster_count if total_cluster_count > 0 else 0,
                'min_cluster_size': min(filtered_cluster_sizes) if filtered_cluster_sizes else 0,
                'max_cluster_size': max(filtered_cluster_sizes) if filtered_cluster_sizes else 0,
                'median_cluster_size': sorted(filtered_cluster_sizes)[len(filtered_cluster_sizes)//2] if filtered_cluster_sizes else 0,
                'components_with_clustering': best_iteration.get('num_clusters', 0),  # Approximate
                'single_cluster_components': sum(1 for size_str, count in filtered_size_dist.items() if int(size_str) == 1),
                'clusters_with_conflicting_pep_idents': best_iteration.get('clusters_with_conflicting_pep_idents', 0),
                'components_with_cluster_pep_ident_duplicates': best_iteration.get('components_with_cluster_pep_ident_duplicates', 0),
                'avg_modularity': best_iteration.get('avg_modularity'),
            }
    
    elif filtered_clusters_data:
        # Process filtered clusters to build size distribution
        filtered_cluster_sizes = []
        filtered_component_clusters = {}
        filtered_modularity_values = []
        filtered_single_cluster_components = 0
        filtered_clusters_with_conflicts = 0
        filtered_components_with_cluster_conflicts = 0
        
        for comp_key, comp_data in filtered_clusters_data.items():
            if isinstance(comp_data, dict) and 'clusters' in comp_data:
                component_id = comp_data.get('component_id', comp_key)
                method = comp_data.get('method', 'unknown')
                modularity = comp_data.get('modularity')
                clusters = comp_data.get('clusters', {})
                
                # Track single cluster components separately (modularity = 0.0)
                if len(clusters) == 1:
                    filtered_single_cluster_components += 1
                elif modularity is not None and modularity != 0.0:
                    filtered_modularity_values.append(modularity)
                
                filtered_component_clusters[component_id] = {
                    'method': method,
                    'clusters': {},
                    'total_clusters': len(clusters)
                }
                
                for cluster_key, cluster_data in clusters.items():
                    num_vertices = cluster_data.get('num_vertices', 0)
                    vertices = cluster_data.get('vertices', [])
                    
                    # Check for conflicting pep_idents
                    has_conflicts = detect_conflicting_pep_idents(vertices)
                    if has_conflicts:
                        filtered_clusters_with_conflicts += 1
                        # Track this conflicting cluster
                        conflicting_clusters_list.append({
                            'component_id': component_id,
                            'cluster_id': cluster_key,
                            'num_vertices': num_vertices,
                            'num_vertices_with_pep_idents': sum(1 for v in vertices if v.get('pep_ident'))
                        })
                    
                    # Check for duplicates (same file appearing multiple times in cluster)
                    file_counts = {}
                    for vertex in vertices:
                        fname = vertex.get('filename', 'unknown')
                        file_counts[fname] = file_counts.get(fname, 0) + 1
                    
                    has_duplicates = any(count > 1 for count in file_counts.values())
                    
                    filtered_cluster_sizes.append(num_vertices)
                    filtered_component_clusters[component_id]['clusters'][cluster_key] = {
                        'size': num_vertices,
                        'has_duplicates': has_duplicates,
                        'file_counts': file_counts,
                        'vertices': vertices
                    }
                
                # Check for conflicts between clusters in this component
                if len(clusters) > 1:  # Only check if component has multiple clusters
                    has_cluster_conflicts = detect_conflicting_pep_idents_between_clusters(clusters)
                    if has_cluster_conflicts:
                        filtered_components_with_cluster_conflicts += 1
        
        # Note: filtered_size_dist is from best iteration (pre-calculated, pre-recovery)
        # Do NOT recalculate from JSON as it would give post-recovery data
        
        # Calculate summary statistics from filtered_cluster_sizes (from actual JSON)
        if filtered_cluster_sizes:
            # Include unpaired singletons in average cluster size calculation
            total_size_sum = sum(filtered_cluster_sizes) + unpaired_singletons_count
            total_cluster_count = len(filtered_cluster_sizes) + unpaired_singletons_count
            
            filtered_cluster_summary = {
                'total_clusters': len(filtered_cluster_sizes),
                'avg_cluster_size': total_size_sum / total_cluster_count if total_cluster_count > 0 else 0,
                'min_cluster_size': min(filtered_cluster_sizes),
                'max_cluster_size': max(filtered_cluster_sizes),
                'median_cluster_size': sorted(filtered_cluster_sizes)[len(filtered_cluster_sizes)//2],
                'components_with_clustering': len(filtered_component_clusters),
                'single_cluster_components': filtered_single_cluster_components,
                'total_vertices_deleted': total_deleted,
                'clusters_with_conflicting_pep_idents': filtered_clusters_with_conflicts,
                'components_with_cluster_pep_ident_duplicates': filtered_components_with_cluster_conflicts,
                'avg_modularity': sum(filtered_modularity_values) / len(filtered_modularity_values) if filtered_modularity_values else None,
                'modularity_values_count': len(filtered_modularity_values)
            }
        
        # Generate filtered cluster distribution HTML
        if filtered_size_dist:
            deleted_note = ""
            if total_deleted > 0:
                deleted_note = f"""
    <div style="background: #fff3cd; padding: 12px; margin-bottom: 15px; border-radius: 4px; color: #856404; font-size: 0.95em; border-left: 4px solid #ffc107;">
        <strong>⚠ Note:</strong> This distribution shows clusters <strong>after filtering</strong> (duplicate file vertices removed).
        <strong style="color: #d9534f;">{total_deleted} vertices were deleted</strong> during filtering due to doubling occurrence of features from the same file within the same cluster.
        These vertices will be shown below as <strong>recovered singletons</strong> after being re-added.
    </div>"""
            
            filtered_cluster_dist_html = f"""
    <div id="filteredClusterDistChart" class="distribution-section">
    {deleted_note}
    """
            
            total_filtered_clusters = sum(filtered_size_dist.values())
            total_filtered_with_unpaired = total_filtered_clusters + unpaired_singletons_count
            
            for size in sorted(filtered_size_dist.keys()):
                count = filtered_size_dist[size]
                percentage = 100 * count / total_filtered_with_unpaired
                
                # Special handling for size-1: show just percentage bar
                if size == 1 and (total_deleted > 0 or unpaired_singletons_count > 0):
                    total_size1 = count + total_deleted + unpaired_singletons_count
                    bar_width = min(95, 100 * total_size1 / total_filtered_with_unpaired)
                    
                    filtered_cluster_dist_data[size] = {
                        'count': count,
                        'percentage': 100 * count / total_filtered_with_unpaired if total_filtered_with_unpaired > 0 else 0,
                        'deleted': total_deleted,
                        'unpaired': unpaired_singletons_count,
                    }
                    
                    filtered_cluster_dist_html += f"""
        <div class="distribution-row" style="display: flex; align-items: center; margin: 8px 0; font-size: 0.9em;">
            <div style="min-width: 120px; font-weight: 500; color: #333;">{size} vertex/vertices</div>
            <div style="flex: 1; margin: 0 10px;">
                <div style="background: #e0e0e0; border-radius: 4px; height: 20px; overflow: hidden;">
                    <div style="background: linear-gradient(90deg, #667eea, #764ba2); width: {bar_width}%; height: 100%;"></div>
                </div>
            </div>
            <div style="min-width: 150px; text-align: right; color: #666;">{total_size1} cluster{'s' if total_size1 != 1 else ''} ({100 * total_size1 / total_filtered_with_unpaired:.1f}%)</div>
        </div>
        """
                else:
                    bar_width = min(95, 100 * count / total_filtered_with_unpaired) if total_filtered_with_unpaired > 0 else 0
                    
                    filtered_cluster_dist_data[size] = {
                        'count': count,
                        'percentage': percentage,
                    }
                    
                    filtered_cluster_dist_html += f"""
        <div class="distribution-row clickable-bar" onclick="showFilteredClusterDistributionDetails({size})" data-cluster-size="{size}" style="cursor: pointer;">
            <div class="dist-label">{size} vertex/vertices</div>
            <div class="dist-bar-container">
                <div class="dist-bar" style="width: {bar_width}%"></div>
            </div>
            <div class="dist-value" style="color: #666;">{count} cluster{'s' if count != 1 else ''} ({percentage:.1f}%)</div>
            <div class="dist-click-hint">📊 Click for details</div>
        </div>
        """
            
            filtered_cluster_dist_html += """
    </div>
    <script id="filteredClusterDistData" type="application/json">
    """ + json.dumps(filtered_cluster_dist_data) + """
    </script>
    """
        
        # Generate deletion summary
        if total_deleted > 0:
            filter_method_display = filter_method if filter_method else "N/A"
            vertices_recovered_display = total_recovered if total_recovered > 0 else "N/A"
            deletion_summary_html = f"""
        <div class="stats-section">
            <h3>Duplicate Vertex Filtering Summary</h3>
            {f'<p style="color: #666; font-size: 0.95em; margin: 10px 0;">Based on optimal resolution: {optimal_resolution_str}</p>' if optimization_result else ''}
            <div class="stats-grid">
                <div class="stat-box">
                    <div class="stat-value">{total_deleted}</div>
                    <div class="stat-label">Vertices Deleted</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{vertices_recovered_display}</div>
                    <div class="stat-label">Vertices Saved/Recovered</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{len(filtered_component_clusters)}</div>
                    <div class="stat-label">Affected Components</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{filtered_cluster_summary.get('total_clusters', 'N/A')}</div>
                    <div class="stat-label">Clusters (Filtered)</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{filter_method_display}</div>
                    <div class="stat-label">Filter Method</div>
                </div>
            </div>
        </div>
    """
    
    # Generate recovered cluster distribution (if vertices were recovered)
    recovered_cluster_dist_html = ""
    recovered_cluster_summary = {}
    
    # Get baseline values from optimization result if available
    if optimization_result and optimization_was_performed:
        iterations = optimization_result.get('iterations', [])
        best_iteration = max(iterations, key=lambda x: x.get('score', float('-inf')))
    else:
        best_iteration = {}
    
    # Use the three pre-calculated distributions from best iteration
    # IMPORTANT: These represent the TRUE STATE after each step of optimization
    # When optimization was performed, these are the ONLY authoritative source:
    #
    # 1. initial_cluster_distribution:  Distribution BEFORE filtering
    #    (clustering completed, no filtering yet)
    #
    # 2. filtered_cluster_distribution: Distribution AFTER filtering, BEFORE recovery
    #    (duplicate vertices removed, but deleted vertices not yet re-added)
    #
    # 3. cluster_recovered_final:       Distribution AFTER recovery (FINAL STATE)
    #    (deleted vertices added back as new single-vertex clusters)
    #    This is what "Recovered Cluster Size Distribution" should display.
    #
    initial_size_dist = best_iteration.get('initial_cluster_distribution', {})
    filtered_size_dist_for_recovery = best_iteration.get('filtered_cluster_distribution', {})
    recovered_size_dist = best_iteration.get('cluster_recovered_final', {})
    
    total_deleted = best_iteration.get('vertices_deleted', 0)
    total_recovered_vertices = best_iteration.get('recovered_vertices', 0) if 'recovered_vertices' in best_iteration else total_recovered
    
    # Load recovered clusters JSON if available (for reference only)
    recovered_clusters_data = {}
    recovered_from_actual_json = False
    
    # NOTE: When optimization was performed, always use the optimization result
    # as the source of truth, since that's what was actually used for scoring.
    # Only try to load external JSON files if optimization did NOT provide the data.
    if total_recovered > 0 and not optimization_was_performed:
        # Only load from external files if we don't have optimization data
        if filtered_clusters_json:
            recovered_clusters_path = filtered_clusters_json.replace('_filtered.json', '_recovered.json')
            try:
                with open(recovered_clusters_path, 'r') as f:
                    recovered_clusters_data = json.load(f)
                    recovered_from_actual_json = True
            except FileNotFoundError:
                pass  # Will use optimization result if available
    
        
    # Build recovered_cluster_dist_html from distribution (this is after recovery)
    if recovered_size_dist:
            total_filtered_clusters_before = sum(int(v) for v in filtered_size_dist_for_recovery.values()) if filtered_size_dist_for_recovery else 0
            total_recovered_clusters = sum(int(v) for v in recovered_size_dist.values()) if recovered_size_dist else 0
            
            recovery_note = ""
            if total_deleted > 0:
                unrecovered = total_deleted - total_recovered_vertices
                recovery_note = f"""
    <div style="background: #d4edda; padding: 12px; margin-bottom: 15px; border-radius: 4px; color: #155724; font-size: 0.95em; border-left: 4px solid #5cb85c;">
        <strong>✓ Recovery Complete:</strong> <strong style="color: #5cb85c;">{total_recovered_vertices} of {total_deleted} deleted vertices</strong> were successfully re-added to existing clusters.
        Remaining <strong>{unrecovered} vertices</strong> became <strong>new singleton clusters</strong> (shown as "1 vertex/vertices" entries below).
        Total clusters increased from <strong>{total_filtered_clusters_before}</strong> (after filtering) to <strong>{total_recovered_clusters}</strong> (after recovery).
    </div>"""
            
            recovered_cluster_dist_html = f"""
    <div id="recoveredClusterDistChart" class="distribution-section">
    {recovery_note}
    """
            
            recovered_cluster_dist_data = {}
            total_recovered_with_unpaired = total_recovered_clusters + unpaired_singletons_count
            
            for size in sorted([int(s) for s in recovered_size_dist.keys()]):
                size_str = str(size)
                count = int(recovered_size_dist[size_str]) if size_str in recovered_size_dist else 0
                percentage = 100 * count / total_recovered_clusters if total_recovered_clusters > 0 else 0
                
                # Special handling for size-1: show net loss (unpaired + unrecovered deleted vertices)
                if size == 1:
                    net_loss_vertices = unpaired_singletons_count + (total_deleted - total_recovered_vertices)
                    if net_loss_vertices > 0:
                        # Calculate proportions for bar display
                        total_size1_with_loss = count + net_loss_vertices
                        recovered_part_width = 100 * count / total_size1_with_loss if total_size1_with_loss > 0 else 0
                        loss_width = 100 - recovered_part_width
                        
                        # Bar width proportion of all clusters (including unpaired + net loss)
                        bar_width = min(95, 100 * total_size1_with_loss / total_recovered_with_unpaired)
                        
                        recovered_cluster_dist_data[size] = {
                            'count': count,
                            'percentage': 100 * count / total_recovered_clusters,
                            'net_loss': net_loss_vertices,
                            'net_loss_breakdown': f'{unpaired_singletons_count} unpaired + {total_deleted - total_recovered_vertices} unrecovered',
                        }
                        
                        recovered_cluster_dist_html += f"""
        <div class="distribution-row" style="display: flex; align-items: center; margin: 8px 0; font-size: 0.9em;">
            <div style="min-width: 120px; font-weight: 500; color: #333;">{size} vertex/vertices</div>
            <div style="flex: 1; margin: 0 10px;">
                <div style="background: #e0e0e0; border-radius: 4px; height: 20px; overflow: hidden; display: flex;">
                    <div style="background: linear-gradient(90deg, #667eea, #764ba2); width: {bar_width}%; height: 100%;"></div>
                </div>
            </div>
            <div style="min-width: 150px; text-align: right; color: #666;">{count} cluster{'s' if count != 1 else ''} ({percentage:.1f}%) </div>
        </div>
        """
                    else:
                        bar_width = min(95, 100 * count / total_recovered_with_unpaired) if total_recovered_with_unpaired > 0 else 0
                        recovered_cluster_dist_data[size] = {
                            'count': count,
                            'percentage': percentage,
                        }
                        
                        recovered_cluster_dist_html += f"""
        <div class="distribution-row" data-cluster-size="{size}" style="display: flex; align-items: center; margin: 8px 0; font-size: 0.9em;">
            <div style="min-width: 120px; font-weight: 500; color: #333;">{size} vertex/vertices</div>
            <div style="flex: 1; margin: 0 10px;">
                <div style="background: #e0e0e0; border-radius: 4px; height: 20px; overflow: hidden;">
                    <div style="background: linear-gradient(90deg, #667eea, #764ba2); width: {bar_width}%; height: 100%;"></div>
                </div>
            </div>
            <div style="min-width: 150px; text-align: right; color: #666;">{count} cluster{'s' if count != 1 else ''} ({percentage:.1f}%)</div>
        </div>
        """
                else:
                    bar_width = min(95, 100 * count / total_recovered_with_unpaired) if total_recovered_with_unpaired > 0 else 0
                    
                    recovered_cluster_dist_data[size] = {
                        'count': count,
                        'percentage': percentage,
                    }
                    
                    recovered_cluster_dist_html += f"""
        <div class="distribution-row" data-cluster-size="{size}" style="display: flex; align-items: center; margin: 8px 0; font-size: 0.9em;">
            <div style="min-width: 120px; font-weight: 500; color: #333;">{size} vertex/vertices</div>
            <div style="flex: 1; margin: 0 10px;">
                <div style="background: #e0e0e0; border-radius: 4px; height: 20px; overflow: hidden;">
                    <div style="background: linear-gradient(90deg, #667eea, #764ba2); width: {bar_width}%; height: 100%;"></div>
                </div>
            </div>
            <div style="min-width: 150px; text-align: right; color: #666;">{count} cluster{'s' if count != 1 else ''} ({percentage:.1f}%)</div>
        </div>
        """
            
            recovered_cluster_dist_html += """
    </div>
    <script id="recoveredClusterDistData" type="application/json">
    """ + json.dumps(recovered_cluster_dist_data) + """
    </script>
    """
            # Build recovered summary from best iteration (if not already built from actual JSON)
            if not recovered_cluster_summary and recovered_size_dist and total_recovered_clusters > 0:
                recovered_cluster_sizes = []
                for size_str, count in recovered_size_dist.items():
                    recovered_cluster_sizes.extend([int(size_str)] * int(count))
                
                # Include unpaired singletons in average cluster size calculation
                total_size_sum = sum(recovered_cluster_sizes) + unpaired_singletons_count
                total_cluster_count = len(recovered_cluster_sizes) + unpaired_singletons_count
                
                recovered_cluster_summary = {
                    'total_clusters': total_recovered_clusters,
                    'avg_cluster_size': total_size_sum / total_cluster_count if total_cluster_count > 0 else 0,
                    'min_cluster_size': min(recovered_cluster_sizes) if recovered_cluster_sizes else 0,
                    'max_cluster_size': max(recovered_cluster_sizes) if recovered_cluster_sizes else 0,
                    'median_cluster_size': sorted(recovered_cluster_sizes)[len(recovered_cluster_sizes)//2] if recovered_cluster_sizes else 0,
                    'components_with_clustering': best_iteration.get('num_clusters', 0),  # Approximate
                    'single_cluster_components': sum(1 for size_str, count in recovered_size_dist.items() if int(size_str) == 1),
                    'total_vertices_recovered': total_recovered,
                    'avg_modularity': best_iteration.get('avg_modularity'),
                }
    
    # Generate identification statistics HTML
    ident_html = ""
    if pep_ident_stats:
        ident_html = f"""
        <div class="stats-section">
            <h3>Peptide Identification Statistics</h3>
            <div class="stats-grid">
                <div class="stat-box">
                    <div class="stat-value">{pep_ident_stats.get('groups_with_common_pep_idents', 0)}</div>
                    <div class="stat-label">Common Identifications</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{pep_ident_stats.get('groups_with_mismatching_pep_idents', 0)}</div>
                    <div class="stat-label">Mismatching Identifications</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{pep_ident_stats.get('groups_with_no_pep_idents', 0)}</div>
                    <div class="stat-label">No Identifications</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{pep_ident_stats.get('groups_with_partial_matching', 0)}</div>
                    <div class="stat-label">Partial Matching</div>
                </div>
            </div>
        </div>
        """
    
    # Generate parameters section HTML
    params_html = ""
    if parameters:
        # Parameter tooltips
        tooltips = {
            'mz_cutoff': 'Maximum m/z value difference allowed for pairing two features',
            'rt_cutoff': 'Maximum retention time difference (seconds) allowed for pairing',
            'edges_cutoff': 'Maximum Euclidean distance in scaled space for feature pairing',
            'best_match_only': 'Keep only the best scoring match for each feature pair',
            'pairing_method': 'Pairing method is hardcoded to KD-tree based nearest-neighbor matching',
            'rt_start_trim': 'Exclude features with RT values below this threshold (seconds)',
            'rt_end_trim': 'Exclude features with RT values above this threshold (seconds)',
            'feature_mode': 'Feature center calculation method: CoM (all isotopes), first_mean (first isotope only), or single_peak (highest intensity peak)',
            'postpair_normalize_coordinates': 'Rescale coordinate axes to balance ranges after pairing for better distance calculation',
            'compare_rt_alignment': 'Generate and compare RT correction formulas across file pairs',
            'rt_correction_mode': 'How to apply RT corrections: none (disabled), pre-kdtree (correct RT before matching), or per-distance (correct during distance calculation)',
            'rt_alignment_method': 'Fitting method for RT correction curves: polynomial, spline, rbf, piecewise, or loess',
            'rt_alignment_polynomial_degree': 'Polynomial degree used for fitting RT correction curves (higher = more complex fitting)',
            'rt_alignment_spline_degree': 'Polynomial degree for each spline segment (1=linear, 2=quadratic, 3=cubic)',
            'rt_alignment_spline_smoothing': 'Smoothing factor for splines (lower = tighter fit, higher = smoother)',
            'rt_alignment_loess_fraction': 'LOESS smoothing fraction - proportion of data for local regression (range 0.01-1.0)',
            'rt_alignment_loess_iterations': 'LOESS iterations for robustness (higher = more robust to outliers)',
            'rt_alignment_decimal_points': 'Decimal places preserved for polynomial coefficients in equations',
        }
        
        params_html = """
        <div class="stats-section">
            <h3>Analysis Configuration</h3>
        """
        
        # Pairing parameters (feature matching configuration)
        pairing_params = [
            ('mz_cutoff', 'M/Z Cutoff'),
            ('rt_cutoff', 'RT Cutoff'),
            ('edges_cutoff', 'Euclidean Cutoff'),
            ('best_match_only', 'Best Match Only'),
            ('pairing_method', 'Pairing Method'),
            ('rt_start_trim', 'RT Start Trim'),
            ('rt_end_trim', 'RT End Trim'),
            ('feature_mode', 'Feature Center Mode'),
            ('postpair_normalize_coordinates', 'Post-Pair Normalization'),
        ]
        
        if any(parameters.get(p[0]) for p in pairing_params):
            params_html += """
            <div class="params-subsection">
                <h4 class="params-subtitle">Pairing Configuration</h4>
                <div class="stats-grid">
            """
            for param_key, param_label in pairing_params:
                value = parameters.get(param_key)
                if value is not None:
                    display_value = str(value) if not isinstance(value, bool) else ('Yes' if value else 'No')
                    tooltip = tooltips.get(param_key, '')
                    tooltip_attr = f' title="{tooltip}"' if tooltip else ''
                    params_html += f"""
                    <div class="stat-box"{tooltip_attr}>
                        <div class="stat-label">{param_label}</div>
                        <div class="stat-value" style="font-size: 1.5em; margin-top: 8px;">{display_value}</div>
                    </div>
                    """
            params_html += """
                </div>
            </div>
            """
        
        # Alignment Configuration parameters - build dynamically based on fitting method
        alignment_params = [
            ('compare_rt_alignment', 'RT Alignment Comparison'),
            ('rt_correction_mode', 'RT Correction Mode'),
            ('rt_alignment_method', 'Fitting Method'),
        ]
        
        # Add method-specific parameters
        fitting_method = parameters.get('rt_alignment_method', 'polynomial')
        if fitting_method == 'polynomial':
            alignment_params.append(('rt_alignment_polynomial_degree', 'Polynomial Degree'))
        elif fitting_method == 'spline':
            alignment_params.extend([
                ('rt_alignment_spline_degree', 'Spline Degree'),
                ('rt_alignment_spline_smoothing', 'Spline Smoothing'),
            ])
        elif fitting_method == 'loess':
            alignment_params.extend([
                ('rt_alignment_loess_fraction', 'LOESS Fraction'),
                ('rt_alignment_loess_iterations', 'LOESS Iterations'),
            ])
        
        alignment_params.append(('rt_alignment_decimal_points', 'Decimal Points'))
        
        if any(parameters.get(p[0]) for p in alignment_params):
            params_html += """
            <div class="params-subsection">
                <h4 class="params-subtitle">Alignment Configuration</h4>
                <div class="stats-grid">
            """
            for param_key, param_label in alignment_params:
                value = parameters.get(param_key)
                if value is not None:
                    display_value = str(value) if not isinstance(value, bool) else ('Yes' if value else 'No')
                    tooltip = tooltips.get(param_key, '')
                    tooltip_attr = f' title="{tooltip}"' if tooltip else ''
                    params_html += f"""
                    <div class="stat-box"{tooltip_attr}>
                        <div class="stat-label">{param_label}</div>
                        <div class="stat-value" style="font-size: 1.5em; margin-top: 8px;">{display_value}</div>
                    </div>
                    """
            params_html += """
                </div>
            </div>
            """
        
        # Clustering parameter tooltips
        clustering_tooltips = {
            'enable_clustering': 'Enable community detection to identify clusters in the feature network',
            'clustering_method': 'Algorithm for community detection: edge_betweenness (slow, accurate), louvain, walktrap, label_propagation',
            'clustering_use_weights': 'Use edge distance weights in clustering algorithm for better community structure',
            'clustering_weight_mode': 'Weight calculation: inverse (1/(1+distance)) gives higher weight to closer features, distance uses raw distances',
            'clustering_resolution_parameter': 'Resolution parameter for Leiden clustering method (0-1 range; lower=fewer larger clusters, higher=more smaller clusters)',
            'clustering_objective_function': 'Objective function used in Leiden clustering: CPM (default) optimizes modularity with multi-level assessment, modularity uses traditional modularity optimization',
        }
        
        # Clustering parameters
        clustering_params = [
            ('enable_clustering', 'Clustering Enabled'),
            ('clustering_method', 'Clustering Method'),
            ('clustering_use_weights', 'Use Weights'),
            ('clustering_weight_mode', 'Weight Mode'),
            ('clustering_resolution_parameter', 'Resolution Parameter'),
            ('clustering_objective_function', 'Objective Function'),
        ]
        
        if parameters.get('enable_clustering') or parameters.get('clustering_method'):
            params_html += """
            <div class="params-subsection">
                <h4 class="params-subtitle">Clustering Configuration</h4>
                <div class="stats-grid">
            """
            for param_key, param_label in clustering_params:
                value = parameters.get(param_key)
                if value is not None:
                    display_value = str(value) if not isinstance(value, bool) else ('Yes' if value else 'No')
                    tooltip = clustering_tooltips.get(param_key, '')
                    tooltip_attr = f' title="{tooltip}"' if tooltip else ''
                    params_html += f"""
                    <div class="stat-box"{tooltip_attr}>
                        <div class="stat-label">{param_label}</div>
                        <div class="stat-value" style="font-size: 1.5em; margin-top: 8px;">{display_value}</div>
                    </div>
                    """
            params_html += """
                </div>
            </div>
            """
        
        params_html += """
        </div>
        """
    
    # Generate graph composition HTML
    comp_html = ""
    if graph_composition:
        comp_summary = graph_composition.get('composition_summary', {})
        comp_dist = graph_composition.get('component_composition', {})
        
        if comp_summary:
            # Create composition distribution table
            comp_rows_html = ""
            for comp_key in sorted(comp_dist.keys(), key=lambda x: tuple(map(int, x.split('_')))):
                count = comp_dist[comp_key]
                parts = comp_key.split('_')
                num_features = parts[0]
                num_files = parts[1]
                comp_rows_html += f"""
                <div class="comp-row">
                    <span class="comp-label">{num_features} features from {num_files} file(s)</span>
                    <span class="comp-count">{count} component{'s' if count != 1 else ''}</span>
                </div>
                """
            
            comp_html = f"""
            <div class="stats-section">
                <h3>Graph Composition</h3>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{comp_summary.get('total_components', 0)}</div>
                        <div class="stat-label">Total Components</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{comp_summary.get('avg_features_per_component', 0):.2f}</div>
                        <div class="stat-label">Avg Features/Component</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{comp_summary.get('avg_files_per_component', 0):.2f}</div>
                        <div class="stat-label">Avg Files/Component</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{comp_summary.get('max_features_in_component', 0)}</div>
                        <div class="stat-label">Max Features/Component</div>
                    </div>
                </div>
                <h4 style="margin-top: 20px; color: #333;">Component Distribution</h4>
                <div class="composition-list">
                    {comp_rows_html}
                </div>
            </div>
            """
    
    # Generate distance statistics section
    dist_section_html = ""
    if dist_stats:
        dist_section_html = f"""
        <div class="stats-section">
            <h3>Distance Statistics</h3>
            <div class="stats-grid">
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['count']}</div>
                    <div class="stat-label">Total Paired Edges</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['mean']:.4f}</div>
                    <div class="stat-label">Mean Distance</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['median']:.4f}</div>
                    <div class="stat-label">Median Distance</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['min']:.4f}</div>
                    <div class="stat-label">Min Distance</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['max']:.4f}</div>
                    <div class="stat-label">Max Distance</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['q3']-dist_stats['q1']:.4f}</div>
                    <div class="stat-label">IQR (Q3-Q1)</div>
                </div>
            </div>
        </div>
        """
    
    ident_percentage = 100 * features_with_ident / total_features if total_features > 0 else 0
    
    # Generate conflicting clusters list HTML
    conflicting_clusters_html = ""
    if conflicting_clusters_list:
        conflicting_rows = ""
        for cluster_info in conflicting_clusters_list:
            component_id = cluster_info.get('component_id', 'N/A')
            cluster_id = cluster_info.get('cluster_id', 'N/A')
            num_vertices = cluster_info.get('num_vertices', 0)
            num_with_peps = cluster_info.get('num_vertices_with_pep_idents', 0)
            
            conflicting_rows += f"""
            <tr>
                <td style="text-align: center;">{component_id}</td>
                <td style="text-align: center;">{cluster_id}</td>
                <td style="text-align: center;">{num_vertices}</td>
                <td style="text-align: center;">{num_with_peps}</td>
            </tr>
            """
        
        conflicting_clusters_html = f"""
        <div class="stats-section">
            <h3>Clusters with Conflicting Pep-Idents</h3>
            <p style="color: #666; font-size: 0.95em; margin: 10px 0;">List of all {len(conflicting_clusters_list)} cluster(s) containing conflicting peptide identifications{f' (optimal resolution: {optimal_resolution_str})' if optimization_result else ''}</p>
            <div style="overflow-x: auto; margin-top: 15px;">
                <table style="width: 100%; border-collapse: collapse;">
                    <thead>
                        <tr style="background: #667eea; color: white;">
                            <th style="padding: 12px; text-align: center; font-weight: 600;">Component ID</th>
                            <th style="padding: 12px; text-align: center; font-weight: 600;">Cluster ID</th>
                            <th style="padding: 12px; text-align: center; font-weight: 600;">Total Vertices</th>
                            <th style="padding: 12px; text-align: center; font-weight: 600;">Vertices with Pep-Idents</th>
                        </tr>
                    </thead>
                    <tbody>
                        {conflicting_rows}
                    </tbody>
                </table>
            </div>
        </div>
        """
    
    # Generate optimization results HTML
    optimization_html = ""
    gridsearch_plots_relative = None
    gridsearch_plots_absolute = None
    if optimization_result:
        iterations = optimization_result.get('iterations', [])
        config = optimization_result.get('config', {})
        
        # Find best iteration by score
        if iterations:
            best_iteration = max(iterations, key=lambda x: x.get('score', float('-inf')))
            optimal_resolution = best_iteration.get('resolution')
        else:
            best_iteration = {}
            optimal_resolution = None
        
        if iterations:
            # Generate gridsearch plots for both relative and absolute values
            gridsearch_plots_relative = create_gridsearch_plots_base64(iterations, total_features_after_filter, show_absolute=False, optimal_resolution=optimal_resolution)
            gridsearch_plots_absolute = create_gridsearch_plots_base64(iterations, total_features_after_filter, show_absolute=True, optimal_resolution=optimal_resolution)
            
            # Helper function to create cluster size distribution HTML for optimization iterations
            def create_resolution_distribution_html(cluster_size_distribution, resolution_id, unpaired_count=0, total_deleted=0, recovered_vertices=0):
                """
                Generate HTML for cluster size distribution chart for a specific resolution.
                
                Args:
                    cluster_size_distribution: Dictionary mapping size -> count
                    resolution_id: Unique identifier for this resolution (for JavaScript)
                    unpaired_count: Number of unpaired singleton vertices
                    total_deleted: Total vertices deleted
                    recovered_vertices: Number of vertices recovered
                
                Returns:
                    HTML string with the distribution chart
                """
                if not cluster_size_distribution:
                    return '<div style="padding: 10px; color: #999;">No cluster distribution data available</div>'
                
                # Create working copy of distribution
                working_dist = dict(cluster_size_distribution)
                total_clusters = sum(working_dist.values())
                if total_clusters == 0:
                    return '<div style="padding: 10px; color: #999;">No clusters to display</div>'
                
                dist_html = '<div style="background: #f8f9fa; border-radius: 6px; padding: 15px;">'
                dist_data = {}
                
                for size in sorted(working_dist.keys(), key=int):
                    count = working_dist[size]
                    print("DDBUG: ", working_dist.keys())
                    percentage = 100 * count / total_clusters
                    bar_width = min(95, percentage)
                    
                    dist_data[size] = {
                        'count': count,
                        'percentage': percentage
                    }
                    
                    dist_html += f"""
                    <div class="distribution-row" style="display: flex; align-items: center; margin-bottom: 8px; font-size: 0.9em;">
                        <div style="min-width: 80px; font-weight: 500; color: #333;">{size} vertex/vertices</div>
                        <div style="flex: 1; margin: 0 10px;">
                            <div style="background: #e0e0e0; border-radius: 4px; height: 20px; overflow: hidden;">
                                <div style="background: linear-gradient(90deg, #667eea, #764ba2); width: {bar_width}%; height: 100%; transition: width 0.3s ease;"></div>
                            </div>
                        </div>
                        <div style="min-width: 150px; text-align: right; color: #666;">{count} cluster{'s' if count != 1 else ''} ({percentage:.1f}%)</div>
                    </div>
                    """
                
                dist_html += '</div>'
                return dist_html
            
            # Create iterations table (sorted by resolution in ascending order)
            iterations_table_rows = ""
            for iter_data in sorted(iterations, key=lambda x: x.get('resolution', 0)):
                iter_num = iter_data.get('iteration', 'N/A')
                resolution = iter_data.get('resolution', 0)
                deleted = iter_data.get('vertices_deleted', 0)
                # Use per-iteration average cluster size (calculated during optimization)
                avg_size = iter_data.get('avg_cluster_size', 0)
                modularity = iter_data.get('avg_modularity')
                score = iter_data.get('score', 0)
                # Use the recovered cluster distribution (after recovery) for display
                cluster_size_dist = iter_data.get('cluster_recovered_final', {})
                
                # Highlight best iteration
                is_best = abs(resolution - best_iteration.get('resolution', 0)) < 0.0001
                row_class = "optimization-row-best" if is_best else "optimization-row"
                
                modularity_str = f"{modularity:.4f}" if modularity is not None else "N/A"
                
                # Create unique ID for this iteration's distribution
                dist_id = f"res_dist_{resolution:.4f}_{iter_num}".replace('.', '_')
                has_distribution = bool(cluster_size_dist)
                
                # Build button and best status HTML outside of f-string to avoid backslash issues
                button_html = f'<button onclick="toggleDistribution(\'{dist_id}\')" style="padding: 5px 10px; background: #667eea; color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 0.85em; font-weight: 500;">📊 Show</button>' if has_distribution else 'N/A'
                best_html = '<span style="font-weight: bold; color: #27ae60; margin-left: 10px;">✓ Best</span>' if is_best else ''
                
                # Main data row with toggle button
                iterations_table_rows += f"""
                <tr class="{row_class}">
                    <td style="text-align: center;">{iter_num}</td>
                    <td style="text-align: center;">{resolution:.4f}</td>
                    <td style="text-align: center;">{deleted}</td>
                    <td style="text-align: center;">{avg_size:.3f}</td>
                    <td style="text-align: center;">{modularity_str}</td>
                    <td style="text-align: center;">{score:.4f}</td>
                    <td style="text-align: center;">
                        {button_html}
                        {best_html}
                    </td>
                </tr>
                """
                
                # Hidden expandable row with cluster distribution
                if has_distribution:
                    recovered_verts = iter_data.get('recovered_vertices', 0)
                    dist_html = create_resolution_distribution_html(cluster_size_dist, dist_id, unpaired_singletons_count, deleted, recovered_verts)
                    iterations_table_rows += f"""
                <tr id="{dist_id}" style="display: none;">
                    <td colspan="7" style="padding: 15px; background: #f5f5f5; border-top: 2px solid #667eea;">
                        <div style="margin-bottom: 10px; font-weight: 600; color: #333;">Cluster Size Distribution @ Resolution {resolution:.4f}</div>
                        {dist_html}
                    </td>
                </tr>
                """
            
            # Pre-compute formatting for best iteration values
            avg_modularity = best_iteration.get('avg_modularity')
            avg_modularity_str = f"{avg_modularity:.4f}" if isinstance(avg_modularity, float) else "N/A"
            
            resolution = best_iteration.get('resolution')
            resolution_str = f"{resolution:.4f}" if isinstance(resolution, (int, float)) else "N/A"
            
            initial_resolution = config.get('initial_resolution')
            initial_resolution_str = f"{initial_resolution:.4f}" if isinstance(initial_resolution, (int, float)) else "N/A"
            
            # Use final recovered avg_cluster_size (after recovery and single-vertex cluster addition)
            avg_cluster_size = recovered_cluster_summary.get('avg_cluster_size', best_iteration.get('avg_cluster_size')) if recovered_cluster_summary else best_iteration.get('avg_cluster_size')
            avg_cluster_size_str = f"{avg_cluster_size:.3f}" if isinstance(avg_cluster_size, (int, float)) else "N/A"
            
            score = best_iteration.get('score')
            score_str = f"{score:.4f}" if isinstance(score, (int, float)) else "N/A"
            
            lambda_val = config.get('lambda')
            lambda_str = f"{lambda_val:.2f}" if isinstance(lambda_val, (int, float)) else "N/A"
            
            # Build plots HTML if available
            plots_html = ""
            if gridsearch_plots_relative or gridsearch_plots_absolute:
                plots_html = """
                <div style="margin-top: 30px;">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;">
                        <h4 style="margin: 0; color: #333;">Grid Search Visualization</h4>
                        <div style="display: flex; gap: 10px;">
                            <button id="toggleRelativeBtn" onclick="toggleOptimizationPlots('relative')" style="padding: 8px 16px; background: #667eea; color: white; border: none; border-radius: 6px; cursor: pointer; font-weight: 500; transition: all 0.3s ease;" class="active-plot-btn">Relative Values</button>
                            <button id="toggleAbsoluteBtn" onclick="toggleOptimizationPlots('absolute')" style="padding: 8px 16px; background: #ddd; color: #333; border: none; border-radius: 6px; cursor: pointer; font-weight: 500; transition: all 0.3s ease;">Absolute Values</button>
                        </div>
                    </div>
                    <div id="optimizationPlotsContainer" style="margin-top: 20px; text-align: center;">
                        """ + (f'<img id="optimizationPlotsRelative" src="data:image/png;base64,{gridsearch_plots_relative}" style="max-width: 100%; height: auto; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1); display: block;">' if gridsearch_plots_relative else '') + """
                        """ + (f'<img id="optimizationPlotsAbsolute" src="data:image/png;base64,{gridsearch_plots_absolute}" style="max-width: 100%; height: auto; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1); display: none;">' if gridsearch_plots_absolute else '') + """
                    </div>
                </div>
                """
            
            optimization_html = f"""
            <div class="stats-section">
                <h3>Resolution Optimization Results</h3>
                <p style="color: #666; font-size: 0.95em; margin: 10px 0;">Gradient descent optimization to find optimal resolution parameter</p>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{config.get('optimization_metric', 'combined').upper()}</div>
                        <div class="stat-label">Optimization Metric</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{resolution_str}</div>
                        <div class="stat-label">Optimal Resolution</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{initial_resolution_str}</div>
                        <div class="stat-label">Initial Resolution</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{avg_cluster_size_str}</div>
                        <div class="stat-label">Avg Size (All Clusters, After Filtering & Recovery)</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{avg_modularity_str}</div>
                        <div class="stat-label">Avg Modularity</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{best_iteration.get('vertices_deleted', 'N/A')}</div>
                        <div class="stat-label">Vertices Deleted</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{best_iteration.get('num_clusters', 'N/A')}</div>
                        <div class="stat-label">Clusters (Optimal)</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{score_str}</div>
                        <div class="stat-label">Optimization Score</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{len(iterations)}</div>
                        <div class="stat-label">Iterations Evaluated</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{lambda_str}</div>
                        <div class="stat-label">Lambda Parameter</div>
                    </div>
                </div>
                {plots_html}
                <h4 style="margin-top: 30px; color: #333;">Optimization Iterations</h4>
                <p style="color: #666; font-size: 0.9em; margin: 10px 0 15px 0;"><em>Each iteration's cluster distribution (shown via "📊 Show" button) includes recovery - i.e., deleted vertices have been re-added as singleton clusters. These distributions should match the final "Recovered Cluster Size Distribution" for the best iteration.</em></p>
                <div style="overflow-x: auto; margin-top: 15px;">
                    <table style="width: 100%; border-collapse: collapse;">
                        <thead>
                            <tr style="background: #667eea; color: white;">
                                <th style="padding: 10px; text-align: center;">Iteration</th>
                                <th style="padding: 10px; text-align: center;">Resolution</th>
                                <th style="padding: 10px; text-align: center;">Deleted</th>
                                <th style="padding: 10px; text-align: center;">Avg Size</th>
                                <th style="padding: 10px; text-align: center;">Modularity</th>
                                <th style="padding: 10px; text-align: center;">Score</th>
                                <th style="padding: 10px; text-align: center;">Status</th>
                            </tr>
                        </thead>
                        <tbody>
                            {iterations_table_rows}
                        </tbody>
                    </table>
                </div>
            </div>
            """
    
    
    # Generate file linkage scores HTML
    file_linkage_scores_html = ""
    if file_linkage_scores:
        table_rows = []
        for filename, score in sorted(file_linkage_scores.items(), key=lambda x: x[1], reverse=True):
            bar_width = score * 100
            row = f"""
                            <tr style="border-bottom: 1px solid #e0e0e0;">
                                <td style="padding: 12px; color: #333;">{filename}</td>
                                <td style="padding: 12px; text-align: right; color: #667eea; font-weight: 600;">{score:.4f}</td>
                                <td style="padding: 12px; text-align: right;">
                                    <div style="width: 200px; height: 24px; background: #e0e0e0; border-radius: 4px; overflow: hidden; display: inline-block;">
                                        <div style="width: {bar_width:.1f}%; height: 100%; background: linear-gradient(90deg, #667eea, #764ba2); transition: width 0.3s;">
                                        </div>
                                    </div>
                                </td>
                            </tr>"""
            table_rows.append(row)
        
        rows_html = "".join(table_rows)
        file_linkage_scores_html = f"""<div class="stats-section">
                <h3>File Linkage Scores</h3>
                <p style="color: #666; font-size: 0.95em; margin: 10px 0;">Score from 0.0 (isolated) to 1.0 (matched across all files). Calculated based on final cluster composition after filtering and recovery.</p>
                <div style="background: #f8f9fa; border-radius: 8px; overflow: hidden;">
                    <table style="width: 100%; border-collapse: collapse;">
                        <thead>
                            <tr style="background: #667eea; color: white;">
                                <th style="padding: 12px; text-align: left; font-weight: 600;">File Name</th>
                                <th style="padding: 12px; text-align: right; font-weight: 600;">Linkage Score</th>
                                <th style="padding: 12px; text-align: right; font-weight: 600;">Visualization</th>
                            </tr>
                        </thead>
                        <tbody>
                            {rows_html}
                        </tbody>
                    </table>
                </div>
            </div>"""
            
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Feature Pairing Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 12px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}
        
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        
        .header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        
        .content {{
            padding: 40px;
        }}
        
        .stats-section {{
            margin-bottom: 40px;
        }}
        
        .stats-section h3 {{
            color: #333;
            margin-bottom: 20px;
            font-size: 1.5em;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .stat-box {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            transition: transform 0.3s ease, box-shadow 0.3s ease;
        }}
        
        .stat-box:hover {{
            transform: translateY(-5px);
            box-shadow: 0 8px 25px rgba(0,0,0,0.15);
        }}
        
        .stat-value {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
            margin-bottom: 10px;
        }}
        
        .stat-label {{
            font-size: 0.9em;
            color: #555;
            font-weight: 500;
        }}
        
        .overview-box {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 40px;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
        }}
        
        .overview-item {{
            text-align: center;
        }}
        
        .overview-value {{
            font-size: 2.5em;
            font-weight: bold;
            margin-bottom: 10px;
        }}
        
        .overview-label {{
            font-size: 0.95em;
            opacity: 0.9;
        }}
        
        .distribution-section {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        
        .distribution-row {{
            display: flex;
            align-items: center;
            margin-bottom: 15px;
            gap: 15px;
        }}
        
        .dist-label {{
            min-width: 100px;
            font-weight: 500;
            color: #333;
        }}
        
        .dist-bar-container {{
            flex: 1;
            height: 30px;
            background: #ddd;
            border-radius: 5px;
            overflow: hidden;
        }}
        
        .dist-bar {{
            height: 100%;
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            transition: width 0.3s ease;
        }}
        
        .dist-value {{
            min-width: 120px;
            text-align: right;
            font-weight: 500;
            color: #666;
        }}
        
        .clickable-bar {{
            transition: all 0.3s ease;
            border-radius: 8px;
            padding: 10px !important;
            margin-bottom: 10px !important;
        }}
        
        .clickable-bar:hover {{
            background: rgba(102, 126, 234, 0.1);
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.2);
            transform: translateX(5px);
        }}
        
        .dist-click-hint {{
            font-size: 0.8em;
            color: #999;
            margin-left: auto;
            margin-right: 10px;
            opacity: 0;
            transition: opacity 0.3s ease;
        }}
        
        .clickable-bar:hover .dist-click-hint {{
            opacity: 1;
        }}
        
        /* Optimization Table Styles */
        .optimization-row {{
            background: #f8f9fa;
            transition: background 0.3s ease;
        }}
        
        .optimization-row:hover {{
            background: #e8eaf6;
        }}
        
        .optimization-row-best {{
            background: #c8e6c9;
            font-weight: bold;
        }}
        
        .optimization-row-best:hover {{
            background: #a5d6a7;
        }}
        
        /* Modal Styles */
        .modal {{
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0, 0, 0, 0.4);
            animation: fadeIn 0.3s ease;
        }}
        
        @keyframes fadeIn {{
            from {{ opacity: 0; }}
            to {{ opacity: 1; }}
        }}
        
        .modal-content {{
            background: white;
            margin: 5% auto;
            padding: 30px;
            border-radius: 12px;
            width: 90%;
            max-width: 600px;
            box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
            animation: slideDown 0.3s ease;
            max-height: 80vh;
            overflow-y: auto;
        }}
        
        @keyframes slideDown {{
            from {{
                transform: translateY(-50px);
                opacity: 0;
            }}
            to {{
                transform: translateY(0);
                opacity: 1;
            }}
        }}
        
        .modal-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
            padding-bottom: 15px;
            border-bottom: 2px solid #667eea;
        }}
        
        .modal-title {{
            font-size: 1.5em;
            font-weight: bold;
            color: #333;
        }}
        
        .close-btn {{
            font-size: 1.5em;
            cursor: pointer;
            color: #999;
            transition: color 0.3s ease;
        }}
        
        .close-btn:hover {{
            color: #333;
        }}
        
        .modal-detail-section {{
            margin-bottom: 25px;
        }}
        
        .modal-detail-title {{
            font-size: 1.1em;
            font-weight: 600;
            color: #667eea;
            margin-bottom: 12px;
            border-left: 4px solid #667eea;
            padding-left: 10px;
        }}
        
        .detail-stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
            margin-bottom: 15px;
        }}
        
        .detail-stat-box {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        
        .detail-stat-num {{
            font-size: 1.8em;
            font-weight: bold;
            color: #667eea;
        }}
        
        .detail-stat-label {{
            font-size: 0.85em;
            color: #555;
            margin-top: 5px;
        }}
        
        .composition-table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border: 1px solid #ddd;
            border-radius: 8px;
            overflow: hidden;
        }}
        
        .composition-table thead {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
        }}
        
        .composition-table th {{
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}
        
        .composition-table td {{
            padding: 12px;
            border-bottom: 1px solid #eee;
        }}
        
        .composition-table tbody tr:hover {{
            background: #f5f5f5;
        }}
        
        .composition-table tbody tr:last-child td {{
            border-bottom: none;
        }}
        
        .has-duplicates {{
            color: #e74c3c !important;
            background-color: #fadbd8 !important;
            font-weight: bold;
        }}
        
        .no-duplicates {{
            color: #27ae60 !important;
            background-color: #d5f4e6 !important;
            font-weight: bold;
        }}
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin-top: 15px;
        }}
        
        .comp-row {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 10px 0;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        .comp-row:last-child {{
            border-bottom: none;
        }}
        
        .comp-label {{
            font-weight: 500;
            color: #333;
        }}
        
        .comp-count {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.9em;
            font-weight: 500;
        }}
        
        .histogram-container {{
            text-align: center;
            margin: 40px 0;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 10px;
        }}
        
        .histogram-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }}
        
        .histogram-container h3 {{
            margin-bottom: 20px;
            color: #333;
        }}
        
        .footer {{
            background: #f8f9fa;
            padding: 20px 40px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
            border-top: 1px solid #ddd;
        }}
        
        .file-appearance-row {{
            display: flex;
            align-items: center;
            margin-bottom: 15px;
            gap: 15px;
            padding: 10px;
            border-radius: 6px;
            background: #f9f9f9;
        }}
        
        .file-app-label {{
            min-width: 150px;
            font-weight: 500;
            color: #333;
            font-size: 0.95em;
            word-break: break-word;
        }}
        
        .file-app-bar-container {{
            flex: 1;
            height: 28px;
            background: #ddd;
            border-radius: 5px;
            overflow: hidden;
        }}
        
        .file-app-bar {{
            height: 100%;
            background: linear-gradient(90deg, #45B7D1 0%, #667eea 100%);
            transition: width 0.3s ease;
        }}
        
        .file-app-value {{
            min-width: 120px;
            text-align: right;
            font-weight: 500;
            color: #666;
            font-size: 0.9em;
        }}
        
        @media (max-width: 768px) {{
            .stats-grid {{
                grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            }}
            
            .header h1 {{
                font-size: 1.8em;
            }}
            
            .distribution-row {{
                flex-direction: column;
                align-items: flex-start;
            }}
            
            .dist-bar-container {{
                width: 100%;
            }}
            
            .file-appearance-row {{
                flex-direction: column;
                align-items: flex-start;
            }}
            
            .file-app-bar-container {{
                width: 100%;
            }}
            
            .file-app-label {{
                min-width: auto;
            }}
        }}
        
        /* Parameters section styling */
        .params-subsection {{
            margin-bottom: 30px;
        }}
        
        .params-subtitle {{
            color: #667eea;
            margin-bottom: 20px;
            font-size: 1.2em;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }}
        
        /* Toggle Button Styles */
        .toggle-button {{
            padding: 10px 18px;
            border: none;
            border-radius: 6px;
            cursor: pointer;
            font-weight: 500;
            font-size: 0.95em;
            transition: all 0.3s ease;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        
        .toggle-button.active {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.3);
        }}
        
        .toggle-button:not(.active) {{
            background: #e8eef5;
            color: #555;
        }}
        
        .toggle-button:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        }}
        
        .toggle-button:active {{
            transform: translateY(0);
        }}
        
        .toggle-buttons-container {{
            display: flex;
            gap: 12px;
            flex-wrap: wrap;
        }}
        
        .active-plot-btn {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
            color: white !important;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>📊 Feature Pairing Report</h1>
            <p>Comprehensive Analysis of Paired Features</p>
        </div>
        
        <div class="content">
            <!-- Overview Section -->
            <div class="overview-box">
                <div class="overview-item">
                    <div class="overview-value">{total_features:,}</div>
                    <div class="overview-label">Paired Groups</div>
                </div>
                <div class="overview-item">
                    <div class="overview-value">{features_with_ident:,}</div>
                    <div class="overview-label">With Identifications</div>
                </div>
                <div class="overview-item">
                    <div class="overview-value">{ident_percentage:.1f}%</div>
                    <div class="overview-label">Identification Rate</div>
                </div>
                <div class="overview-item">
                    <div class="overview-value">{avg_intra_distance:.4f}</div>
                    <div class="overview-label">Avg Intra-group Distance</div>
                </div>
            </div>
            
            <!-- Pairing Parameters Section -->
            {params_html}
            
            <!-- Filtering Parameters Section -->
            {filter_html}
            
            <!-- File Linkage Scores Section -->
            {file_linkage_scores_html}
            
            <!-- Feature Summary Section -->
            <div class="stats-section">
                <h3>Feature Summary</h3>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{total_features_before_filter:,}</div>
                        <div class="stat-label">Features Before Filtering</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{total_features_after_filter:,}</div>
                        <div class="stat-label">Features After Filtering</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{max_intra_distance:.4f}</div>
                        <div class="stat-label">Max Intra-group Distance</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{min_intra_distance:.4f}</div>
                        <div class="stat-label">Min Intra-group Distance</div>
                    </div>
                </div>
            </div>
            
            <!-- Distance Statistics -->
            {dist_section_html}
            
            <!-- Distance Histogram -->
            {f'<div class="histogram-container"><h3>Distance Distribution</h3><img src="data:image/png;base64,{histogram_base64}" alt="Distance Distribution Histogram"></div>' if histogram_base64 else ''}
            
            <!-- Resolution Optimization Results (if available) -->
            {optimization_html}
            
            <!-- File Matching Distribution -->
            <div class="stats-section">
                <h3>File Matching Distribution</h3>
                <div class="distribution-section">
                    {file_dist_html}
                </div>
            </div>
            
            <!-- Cluster Size Distribution -->
            {f'''<div class="stats-section">
                <h3>Cluster Size Distribution</h3>
                {f'<p style="color: #666; font-size: 0.95em; margin: 10px 0;">Based on optimal resolution: {optimal_resolution_str}</p>' if optimization_result else ''}
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("total_clusters", 0) + unpaired_singletons_count:,}</div>
                        <div class="stat-label">Total Clusters (incl. unpaired)</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("avg_cluster_size", 0):.2f}</div>
                        <div class="stat-label">Avg Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("min_cluster_size", 0)}</div>
                        <div class="stat-label">Min Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("max_cluster_size", 0)}</div>
                        <div class="stat-label">Max Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("median_cluster_size", 0)}</div>
                        <div class="stat-label">Median Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("components_with_clustering", 0)}</div>
                        <div class="stat-label">Components with Clustering</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("single_cluster_components", 0)}</div>
                        <div class="stat-label">Components with one Cluster</div>
                    </div>
                    {'<div class="stat-box">' + f'<div class="stat-value">{cluster_summary.get("avg_modularity", 0):.3f}</div>' + '<div class="stat-label">Average Modularity</div>' + '</div>' if cluster_summary.get("avg_modularity") is not None else ''}
                </div>
                <div class="distribution-section">
                    {cluster_dist_html}
                </div>
            </div>''' if cluster_summary else ''}
            
            <!-- Filtered Cluster Distribution (if filtering was applied) -->
            {f'''
            <div id="filteredClustersSection" class="stats-section">
                <h3>Filtered Cluster Size Distribution</h3>
                <p style="color: #666; font-size: 0.95em; margin: 10px 0;">After removing duplicate file vertices{f' (optimal resolution: {optimal_resolution_str})' if optimization_result else ''}</p>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{filtered_cluster_summary.get("total_clusters", 0) + unpaired_singletons_count:,}</div>
                        <div class="stat-label">Total Clusters (incl. unpaired)</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{filtered_cluster_summary.get("avg_cluster_size", 0):.2f}</div>
                        <div class="stat-label">Average Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{filtered_cluster_summary.get("min_cluster_size", 0)}</div>
                        <div class="stat-label">Min Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{filtered_cluster_summary.get("max_cluster_size", 0)}</div>
                        <div class="stat-label">Max Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{filtered_cluster_summary.get("median_cluster_size", 0)}</div>
                        <div class="stat-label">Median Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{filtered_cluster_summary.get("components_with_clustering", 0)}</div>
                        <div class="stat-label">Components with Clustering</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{filtered_cluster_summary.get("single_cluster_components", 0)}</div>
                        <div class="stat-label">Components with one Cluster</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{filtered_cluster_summary.get("total_vertices_deleted", total_deleted)}</div>
                        <div class="stat-label">Vertices Deleted (Count)</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{vertices_deleted_percentage:.2f}%</div>
                        <div class="stat-label">Vertices Deleted (Percentage)</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{filtered_cluster_summary.get("clusters_with_conflicting_pep_idents", 0)}</div>
                        <div class="stat-label">Clusters with Conflicting Pep-Idents</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{filtered_cluster_summary.get("components_with_cluster_pep_ident_duplicates", 0)}</div>
                        <div class="stat-label">Components with Cluster Pep-Ident Duplicates</div>
                    </div>
                    {'<div class="stat-box">' + f'<div class="stat-value">{filtered_cluster_summary.get("avg_modularity", 0):.3f}</div>' + '<div class="stat-label">Average Modularity</div>' + '</div>' if filtered_cluster_summary.get("avg_modularity") is not None else ''}
                </div>
                <div class="distribution-section">
                    {filtered_cluster_dist_html}
                </div>
            </div>''' if filtered_cluster_summary else ''}
            
            <!-- Deletion Summary -->
            {deletion_summary_html}
            
            <!-- Recovered Cluster Distribution (if vertices were recovered) -->
            {f'''
            <div id="recoveredClustersSection" class="stats-section">
                <h3>Recovered Cluster Size Distribution</h3>
                <p style="color: #666; font-size: 0.95em; margin: 10px 0;">After re-adding saved duplicate file vertices{f' at optimal resolution {optimal_resolution_str}' if optimization_result else ''}. <em style="font-weight: 600; color: #555;">Source: {"✓ Optimization result (authoritative)" if optimization_was_performed else "External analysis"}</em></p>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{recovered_cluster_summary.get("total_clusters", 0) :,}</div>
                        <div class="stat-label">Total Clusters (incl. unpaired)</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{recovered_cluster_summary.get("avg_cluster_size", 0):.2f}</div>
                        <div class="stat-label">Average Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{recovered_cluster_summary.get("min_cluster_size", 0)}</div>
                        <div class="stat-label">Min Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{recovered_cluster_summary.get("max_cluster_size", 0)}</div>
                        <div class="stat-label">Max Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{recovered_cluster_summary.get("median_cluster_size", 0)}</div>
                        <div class="stat-label">Median Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{recovered_cluster_summary.get("components_with_clustering", 0)}</div>
                        <div class="stat-label">Components with Clustering</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{recovered_cluster_summary.get("single_cluster_components", 0)}</div>
                        <div class="stat-label">Components with one Cluster</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{recovered_cluster_summary.get("total_vertices_recovered", 0)}</div>
                        <div class="stat-label">Vertices Re-added</div>
                    </div>
                    {'<div class="stat-box">' + f'<div class="stat-value">{recovered_cluster_summary.get("avg_modularity", 0):.3f}</div>' + '<div class="stat-label">Average Modularity</div>' + '</div>' if recovered_cluster_summary.get("avg_modularity") is not None else ''}
                </div>
                <div class="distribution-section">
                    {recovered_cluster_dist_html}
                </div>
            </div>''' if recovered_cluster_summary else ''}
            
            <!-- Clusters with Conflicting Pep-Idents -->
            {conflicting_clusters_html}
            
            <!-- Graph Composition Statistics -->
            {comp_html}
            
            <!-- Peptide Identification Statistics -->
            {ident_html}
            
            <!-- Analysis Configuration -->
            {params_html}
        </div>
        
        <div class="footer">
            Generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')} | Feature Pairing Analysis Report
        </div>
    </div>
    
    <!-- Modal for file distribution details -->
    <div id="fileDistModal" class="modal">
        <div class="modal-content">
            <div class="modal-header">
                <div class="modal-title" id="modalTitle">File Matching Details</div>
                <span class="close-btn" onclick="closeFileDistModal()">&times;</span>
            </div>
            <div id="modalBody"></div>
        </div>
    </div>
    
    <script>
        // Load file distribution data
        const fileDistData = JSON.parse(document.getElementById('fileDistData').textContent);
        
        function showFileDistributionDetails(numFiles) {{
            const data = fileDistData[numFiles];
            if (!data) return;
            
            // Build modal content
            let modalBody = document.getElementById('modalBody');
            
            let html = `
                <div class="modal-detail-section">
                    <div class="modal-detail-title">Overview</div>
                    <div class="detail-stats">
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.total}}</div>
                            <div class="detail-stat-label">Total Groups</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.percentage.toFixed(1)}}%</div>
                            <div class="detail-stat-label">Of All Groups</div>
                        </div>
                    </div>
                </div>
                
                <div class="modal-detail-section">
                    <div class="modal-detail-title">File Uniqueness & Multiplicity</div>
                    <div class="detail-stats">
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.all_unique}}</div>
                            <div class="detail-stat-label">All Files Unique</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.with_duplicates}}</div>
                            <div class="detail-stat-label">With Duplicates</div>
                        </div>
                    </div>
                    <p style="color: #666; font-size: 0.9em; margin-top: 10px;">
                        <strong>All Files Unique:</strong> Each of the ${{numFiles}} files appears exactly once<br>
                        <strong>With Duplicates:</strong> At least one file appears multiple times in the group
                    </p>
                </div>
            `;
            
            if (data.composition && data.composition.length > 0) {{
                html += `
                    <div class="modal-detail-section">
                        <div class="modal-detail-title">File Composition Breakdown</div>
                        <table class="composition-table">
                            <thead>
                                <tr>
                                    <th>Files Involved</th>
                                    <th style="text-align: right;">Count</th>
                                    <th style="text-align: right;">Percentage</th>
                                </tr>
                            </thead>
                            <tbody>
                `;
                
                data.composition.forEach(comp => {{
                    html += `
                        <tr>
                            <td>${{comp.composition}}</td>
                            <td style="text-align: right;"><strong>${{comp.count}}</strong></td>
                            <td style="text-align: right;">${{comp.percentage.toFixed(1)}}%</td>
                        </tr>
                    `;
                }});
                
                html += `
                            </tbody>
                        </table>
                    </div>
                `;
            }}
            
            if (data.multiplicity && data.multiplicity.length > 0) {{
                html += `
                    <div class="modal-detail-section">
                        <div class="modal-detail-title">Feature Multiplicity per File</div>
                        <p style="color: #666; font-size: 0.85em; margin-bottom: 12px;">
                            Shows how many features appear from each file: e.g., "2x2,1x1" means 2 files with 2 features each, 1 file with 1 feature
                        </p>
                        <table class="composition-table">
                            <thead>
                                <tr>
                                    <th>Multiplicity Pattern</th>
                                    <th style="text-align: right;">Count</th>
                                    <th style="text-align: right;">Percentage</th>
                                </tr>
                            </thead>
                            <tbody>
                `;
                
                data.multiplicity.forEach(mult => {{
                    html += `
                        <tr>
                            <td><strong>${{mult.pattern}}</strong></td>
                            <td style="text-align: right;"><strong>${{mult.count}}</strong></td>
                            <td style="text-align: right;">${{mult.percentage.toFixed(1)}}%</td>
                        </tr>
                    `;
                }});
                
                html += `
                            </tbody>
                        </table>
                    </div>
                `;
            }}
            
            if (data.file_appearance && data.file_appearance.length > 0) {{
                html += `
                    <div class="modal-detail-section">
                        <div class="modal-detail-title">File Appearance in Matches</div>
                        <p style="color: #666; font-size: 0.85em; margin-bottom: 15px;">
                            Shows which files appear in what percentage of the ${{data.total}} matching groups
                        </p>
                `;
                
                data.file_appearance.forEach(file => {{
                    const barWidth = Math.min(95, file.percentage);
                    html += `
                        <div class="file-appearance-row">
                            <div class="file-app-label">${{file.filename}}</div>
                            <div class="file-app-bar-container">
                                <div class="file-app-bar" style="width: ${{barWidth}}%"></div>
                            </div>
                            <div class="file-app-value">${{file.count}} (${{file.percentage.toFixed(1)}}%)</div>
                        </div>
                    `;
                }});
                
                html += `
                    </div>
                `;
            }}
            
            document.getElementById('modalTitle').textContent = `File Distribution Details: ${{numFiles}} File(s)`;
            modalBody.innerHTML = html;
            
            // Show modal
            document.getElementById('fileDistModal').style.display = 'block';
        }}
        
        function closeFileDistModal() {{
            document.getElementById('fileDistModal').style.display = 'none';
        }}
        
        // Close modal when clicking outside
        window.onclick = function(event) {{
            const modal = document.getElementById('fileDistModal');
            if (event.target === modal) {{
                modal.style.display = 'none';
            }}
        }}
    </script>
    
    <!-- Modal for cluster distribution details -->
    <div id="clusterDistModal" class="modal">
        <div class="modal-content">
            <div class="modal-header">
                <div class="modal-title" id="clusterModalTitle">Cluster Size Distribution Details</div>
                <span class="close-btn" onclick="closeClusterDistModal()">&times;</span>
            </div>
            <div id="clusterModalBody"></div>
        </div>
    </div>
    
    <script>
        // Load cluster distribution data (if available)
        let clusterDistData = {{}};
        const clusterDataElement = document.getElementById('clusterDistData');
        if (clusterDataElement) {{
            clusterDistData = JSON.parse(clusterDataElement.textContent);
        }}
        
        function showClusterDistributionDetails(clusterSize) {{
            const data = clusterDistData[clusterSize];
            if (!data) return;
            
            let modalBody = document.getElementById('clusterModalBody');
            
            const duplicatePercentage = data.count > 0 ? (100 * data.clusters_with_duplicates / data.count) : 0;
            const cleanPercentage = data.count > 0 ? (100 * data.clusters_without_duplicates / data.count) : 0;
            
            let html = `
                <div class="modal-detail-section">
                    <div class="modal-detail-title">Overview</div>
                    <div class="detail-stats">
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.count}}</div>
                            <div class="detail-stat-label">Clusters of This Size</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.percentage.toFixed(1)}}%</div>
                            <div class="detail-stat-label">Of All Clusters</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{clusterSize}}</div>
                            <div class="detail-stat-label">Feature(s) per Cluster</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.avg_files_per_cluster.toFixed(2)}}</div>
                            <div class="detail-stat-label">Avg Files per Cluster</div>
                        </div>
                    </div>
                </div>
                
                <div class="modal-detail-section">
                    <div class="modal-detail-title">Duplicate Detection Statistics</div>
                    <div class="detail-stats">
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.clusters_with_duplicates}}</div>
                            <div class="detail-stat-label">Clusters with Duplicates</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{duplicatePercentage.toFixed(1)}}%</div>
                            <div class="detail-stat-label">Percentage</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.clusters_without_duplicates}}</div>
                            <div class="detail-stat-label">Clean Clusters</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{cleanPercentage.toFixed(1)}}%</div>
                            <div class="detail-stat-label">Percentage</div>
                        </div>
                    </div>
                    <p style="color: #666; font-size: 0.85em; margin-top: 10px;">
                        <strong>Clusters with Duplicates:</strong> Contain features from the same file appearing multiple times<br>
                        <strong>Clean Clusters:</strong> Each file appears at most once
                    </p>
                </div>
                
                <div class="modal-detail-section">
                    <div class="modal-detail-title">Component Distribution</div>
                    <div class="detail-stats">
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.components_represented}}</div>
                            <div class="detail-stat-label">Components with This Cluster Size</div>
                        </div>
                    </div>
                    <p style="color: #666; font-size: 0.85em; margin-top: 10px;">
                        This cluster size is represented across ${{data.components_represented}} different components in the graph.
                    </p>
                </div>
            `;
            
            const pluralS = data.count !== 1 ? 's' : '';
            document.getElementById('clusterModalTitle').textContent = `Cluster Size ${{clusterSize}} (${{data.count}} cluster${{pluralS}})`;
            modalBody.innerHTML = html;
            
            // Show modal
            document.getElementById('clusterDistModal').style.display = 'block';
        }}
        
        function closeClusterDistModal() {{
            document.getElementById('clusterDistModal').style.display = 'none';
        }}
        
        // Close cluster modal when clicking outside
        window.onclick = function(event) {{
            const fileModal = document.getElementById('fileDistModal');
            const clusterModal = document.getElementById('clusterDistModal');
            if (event.target === fileModal) {{
                fileModal.style.display = 'none';
            }}
            if (event.target === clusterModal) {{
                clusterModal.style.display = 'none';
            }}
        }}
        
        // Toggle cluster size distribution visibility for optimization iterations
        function toggleDistribution(distributionId) {{
            const row = document.getElementById(distributionId);
            const button = event.target;
            
            if (row) {{
                const isVisible = row.style.display !== 'none';
                row.style.display = isVisible ? 'none' : 'table-row';
                button.textContent = isVisible ? '📊 Show' : '📊 Hide';
                button.style.background = isVisible ? '#667eea' : '#e74c3c';
            }}
        }}
        
        // Toggle optimization plots between relative and absolute values
        function toggleOptimizationPlots(viewType) {{
            const relativeImg = document.getElementById('optimizationPlotsRelative');
            const absoluteImg = document.getElementById('optimizationPlotsAbsolute');
            const relativeBtn = document.getElementById('toggleRelativeBtn');
            const absoluteBtn = document.getElementById('toggleAbsoluteBtn');
            
            if (viewType === 'relative') {{
                if (relativeImg) relativeImg.style.display = 'block';
                if (absoluteImg) absoluteImg.style.display = 'none';
                relativeBtn.style.background = '#667eea';
                relativeBtn.style.color = 'white';
                absoluteBtn.style.background = '#ddd';
                absoluteBtn.style.color = '#333';
            }} else {{
                if (relativeImg) relativeImg.style.display = 'none';
                if (absoluteImg) absoluteImg.style.display = 'block';
                relativeBtn.style.background = '#ddd';
                relativeBtn.style.color = '#333';
                absoluteBtn.style.background = '#667eea';
                absoluteBtn.style.color = 'white';
            }}
        }}
    </script>
</body>
</html>
"""
    
    with open(output_html, 'w') as f:
        f.write(html_content)



def main():
    parser = argparse.ArgumentParser(
        description="Generate summary report of feature pairing results"
    )
    parser.add_argument("--paired_json", required=True, help="Path to paired features JSON file")
    parser.add_argument("--output_summary", required=True, help="Output summary text file path")
    parser.add_argument("--output_stats", required=True, help="Output statistics CSV file path")
    parser.add_argument("--edges_json", default=None, help="Path to edges JSON file (for distance statistics and histogram)")
    parser.add_argument("--output_html", default=None, help="Output HTML report file path")
    parser.add_argument("--graph_composition_json", default=None, help="Path to graph composition JSON file (from build_network_graph --output_composition)")
    parser.add_argument("--clusters_json", default=None, help="Path to clusters JSON file (from build_network_graph --output-clusters-json)")
    parser.add_argument("--filtered_clusters_json", default=None, help="Path to filtered clusters JSON file (from build_network_graph with --filter-duplicate-file-vertices)")
    parser.add_argument("--deleted_vertices_json", default=None, help="Path to deleted vertices report JSON file (from build_network_graph filtering)")
    parser.add_argument("--resolution_iterations_json", default=None, help="Path to resolution optimization iterations JSON file (from build_network_graph with --auto-select-resolution)")
    
    # Pairing parameters
    parser.add_argument("--mz_cutoff", type=float, default=None, help="Maximum m/z difference for pairing")
    parser.add_argument("--rt_cutoff", type=float, default=None, help="Maximum RT difference for pairing")
    parser.add_argument("--edges_cutoff", type=float, default=None, help="Maximum euclidean distance cutoff for pairing")
    # Note: best_match_only is now hardcoded to True
    # Note: Pairing method is now hardcoded to KD-tree based matching
    
    # Trimming parameters
    parser.add_argument("--rt_start_trim", type=float, default=None, help="Seconds to trim from start of retention time")
    parser.add_argument("--rt_end_trim", type=float, default=None, help="Maximum RT value to keep")
    
    # RT Correction parameters
    parser.add_argument("--compare_rt_alignment", action="store_true", help="Enable RT alignment comparison and correction formula generation")
    parser.add_argument("--rt_correction_mode", default='none', help="RT correction mode: 'none', 'pre-kdtree', or 'per-distance'")
    parser.add_argument("--rt_alignment_method", default='polynomial', help="Fitting method: 'polynomial', 'spline', 'rbf', 'piecewise', or 'loess'")
    parser.add_argument("--mz_rt_weight_ratio", type=float, default=1.0, help="Weight ratio for RT relative to m/z in distance calculations (default: 1.0 for equal weighting)")
    parser.add_argument("--rt_alignment_polynomial_degree", type=int, default=8, help="Polynomial degree for fitting RT correction curves")
    parser.add_argument("--rt_alignment_spline_degree", type=int, default=3, help="Spline degree for spline fitting (only used if method='spline')")
    parser.add_argument("--rt_alignment_spline_smoothing", type=float, default=100, help="Smoothing factor for splines (only used if method='spline')")
    parser.add_argument("--rt_alignment_loess_fraction", type=float, default=0.3, help="LOESS smoothing fraction (only used if method='loess')")
    parser.add_argument("--rt_alignment_loess_iterations", type=int, default=3, help="LOESS iterations for robustness (only used if method='loess')")
    parser.add_argument("--rt_alignment_decimal_points", type=int, default=35, help="Decimal places for polynomial coefficients")
    
    # Feature extraction parameters
    parser.add_argument("--feature_mode", default=None, help="Feature center calculation mode: CoM (all isotopes), first_mean (first isotope), or single_peak (strongest peak)")
    
    # Network graph parameters
    parser.add_argument("--graph_edge_cutoff", type=float, default=None, help="Maximum euclidean distance for graph edges")
    parser.add_argument("--graph_mz_cutoff", type=float, default=None, help="Maximum m/z difference for graph")
    parser.add_argument("--graph_rt_cutoff", type=float, default=None, help="Maximum RT difference for graph")
    parser.add_argument("--graph_test_fraction", type=float, default=1, help="Fraction of edges to use for testing")
    parser.add_argument("--graph_layout_engine", default=None, help="Graphviz layout engine (dot, sfdp, fdp, neato, twopi, circo)")
    
    # Clustering parameters
    parser.add_argument("--enable_clustering", action="store_true", help="Enable community detection clustering")
    parser.add_argument("--clustering_method", default="edge_betweenness", help="Clustering method (louvain, walktrap, label_propagation, edge_betweenness)")
    parser.add_argument("--clustering_use_weights", action="store_true", help="Use edge distances as weights in clustering")
    parser.add_argument("--clustering_weight_mode", default="inverse", help="Weight mode (inverse or distance)")
    parser.add_argument("--clustering_resolution_parameter", type=float, default=None, help="Resolution parameter for Leiden clustering (e.g., 0.5)")
    parser.add_argument("--clustering_objective_function", default=None, help="Objective function for Leiden clustering (CPM or modularity)")
    
    # Post-pairing parameters
    parser.add_argument("--postpair_normalize_coordinates", action="store_true", help="Rescale axes to balance ranges after pairing")
    
    # Input files tracking
    parser.add_argument("--input_files", nargs='*', default=[], help="List of input feature files used in pairing")
    
    args = parser.parse_args()
    
    # Build parameters dict
    parameters = {
        'mz_cutoff': args.mz_cutoff,
        'rt_cutoff': args.rt_cutoff,
        'edges_cutoff': args.edges_cutoff,
        'best_match_only': True,
        'pairing_method': 'kdtree',  # Hardcoded to KD-tree matching
        'rt_start_trim': args.rt_start_trim,
        'rt_end_trim': args.rt_end_trim,
        'compare_rt_alignment': args.compare_rt_alignment,
        'rt_correction_mode': args.rt_correction_mode,
        'rt_alignment_method': args.rt_alignment_method,
        'rt_alignment_polynomial_degree': args.rt_alignment_polynomial_degree,
        'rt_alignment_spline_degree': args.rt_alignment_spline_degree,
        'rt_alignment_spline_smoothing': args.rt_alignment_spline_smoothing,
        'rt_alignment_loess_fraction': args.rt_alignment_loess_fraction,
        'rt_alignment_loess_iterations': args.rt_alignment_loess_iterations,
        'rt_alignment_decimal_points': args.rt_alignment_decimal_points,
        'postpair_normalize_coordinates': args.postpair_normalize_coordinates,
        'graph_edge_cutoff': args.graph_edge_cutoff,
        'graph_mz_cutoff': args.graph_mz_cutoff,
        'graph_rt_cutoff': args.graph_rt_cutoff,
        'graph_test_fraction': args.graph_test_fraction,
        'graph_layout_engine': args.graph_layout_engine,
        'enable_clustering': args.enable_clustering,
        'clustering_method': args.clustering_method,
        'clustering_use_weights': args.clustering_use_weights,
        'clustering_weight_mode': args.clustering_weight_mode,
        'clustering_resolution_parameter': args.clustering_resolution_parameter,
        'clustering_objective_function': args.clustering_objective_function,
        'input_files': args.input_files,
    }
    
    print(f"\n{'='*70}")
    print(f"Generating Pairing Report")
    print(f"{'='*70}")
    
    generate_pairing_report(args.paired_json, args.output_summary, args.output_stats, args.edges_json, args.output_html, args.graph_composition_json, args.clusters_json, args.filtered_clusters_json, args.deleted_vertices_json, args.resolution_iterations_json, parameters)
    
    print(f"\n✓ Report saved to {os.path.basename(args.output_summary)}")
    print(f"✓ Statistics saved to {os.path.basename(args.output_stats)}")
    if args.output_html:
        print(f"✓ HTML report saved to {os.path.basename(args.output_html)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Pair matching features across multiple files.
Refactored from map_mzml_batch_tsv.py
"""

import argparse
import pickle
import json
import numpy as np
import os
import sys
from pathlib import Path
from typing import Dict, List
from scipy.spatial import cKDTree
from pdb import set_trace as bp

def filter_edges_by_distance_cutoff(edges: Dict, cutoff: float) -> Dict:
    """
    Filter edges dictionary by euclidean distance cutoff.
    
    Args:
        edges: Dictionary of edges {ref_file_idx: [edge_dict, ...]}
        cutoff: Maximum distance threshold
    
    Returns:
        Filtered edges dictionary
    """
    filtered_edges = {}
    for file_idx, edge_list in edges.items():
        filtered_edges[file_idx] = [e for e in edge_list if e.get('distance', 0) <= cutoff]
    return filtered_edges


def filter_edges_by_coordinate_cutoffs(edges: Dict, mz_cutoff: float, rt_cutoff: float, all_feature_data_lists: List[List[Dict]] = None) -> Dict:
    """
    Filter edges dictionary by separate m/z and RT coordinate cutoffs.
    Uses coordinate differences instead of euclidean distance.
    
    Args:
        edges: Dictionary of edges {ref_file_idx: [edge_dict, ...]}
        mz_cutoff: Maximum m/z (x-axis) difference
        rt_cutoff: Maximum RT (y-axis) difference
        all_feature_data_lists: List of feature data lists for coordinate lookup
    
    Returns:
        Filtered edges dictionary
    """
    if all_feature_data_lists is None:
        print("  ✗ WARNING: Feature data not provided for coordinate filtering, skipping")
        return edges
    
    filtered_edges = {}
    edges_filtered_count = 0
    edges_kept_count = 0
    
    for file_idx, edge_list in edges.items():
        filtered_list = []
        for edge in edge_list:
            try:
                # Look up the feature data using indices
                ref_file_idx = edge.get('ref_file_idx')
                match_file_idx = edge.get('match_file_idx')
                current_feat_idx = edge.get('current_file', {}).get('feature_idx')
                matched_feat_idx = edge.get('matched_file', {}).get('feature_idx')
                
                # Validate indices are within bounds
                if (ref_file_idx is None or match_file_idx is None or 
                    current_feat_idx is None or matched_feat_idx is None):
                    filtered_list.append(edge)
                    edges_kept_count += 1
                    continue
                
                if (ref_file_idx >= len(all_feature_data_lists) or 
                    match_file_idx >= len(all_feature_data_lists)):
                    filtered_list.append(edge)
                    edges_kept_count += 1
                    continue
                
                ref_features = all_feature_data_lists[ref_file_idx]
                match_features = all_feature_data_lists[match_file_idx]
                
                if (current_feat_idx >= len(ref_features) or 
                    matched_feat_idx >= len(match_features)):
                    filtered_list.append(edge)
                    edges_kept_count += 1
                    continue
                
                # Get the actual features
                ref_feature = ref_features[current_feat_idx]
                matched_feature = match_features[matched_feat_idx]
                
                # Calculate coordinate differences (absolute values)
                mz_diff = abs(ref_feature['x_center'] - matched_feature['x_center'])
                rt_diff = abs(ref_feature['y_center'] - matched_feature['y_center'])
                
                # Filter based on coordinate cutoffs
                if mz_diff <= mz_cutoff and rt_diff <= rt_cutoff:
                    filtered_list.append(edge)
                    edges_kept_count += 1
                else:
                    edges_filtered_count += 1
            except (KeyError, IndexError, TypeError):
                # If any error occurs during lookup, keep the edge
                filtered_list.append(edge)
                edges_kept_count += 1
        
        filtered_edges[file_idx] = filtered_list
    
    print(f"  Coordinate filtering results: {edges_kept_count} kept, {edges_filtered_count} removed")
    return filtered_edges


def normalize_edge_distances(edges: Dict) -> Dict:
    """
    Normalize all distance values in edges dictionary to [0, 1] range for visualization.
    Min distance maps to 0, max distance maps to 1.
    
    Args:
        edges: Dictionary of edges {ref_file_idx: [edge_dict, ...]}
    
    Returns:
        Dictionary with normalized distance values
    """
    # Collect all distances
    all_distances = []
    for file_idx, edge_list in edges.items():
        for edge in edge_list:
            all_distances.append(edge.get('distance', 0))
    
    if not all_distances:
        return edges
    
    min_dist = min(all_distances)
    max_dist = max(all_distances)
    
    # If all distances are the same, normalize to 0.5 (middle value)
    if min_dist == max_dist:
        normalized_edges = {}
        for file_idx, edge_list in edges.items():
            normalized_list = []
            for edge in edge_list:
                edge_copy = edge.copy()
                edge_copy['distance'] = 0.5
                normalized_list.append(edge_copy)
            normalized_edges[file_idx] = normalized_list
        return normalized_edges
    
    # Normalize distances to [0, 1]
    normalized_edges = {}
    for file_idx, edge_list in edges.items():
        normalized_list = []
        for edge in edge_list:
            edge_copy = edge.copy()
            original_distance = edge.get('distance', 0)
            normalized_distance = (original_distance - min_dist) / (max_dist - min_dist)
            edge_copy['distance'] = float(normalized_distance)
            normalized_list.append(edge_copy)
        normalized_edges[file_idx] = normalized_list
    
    return normalized_edges


def apply_postpair_coordinate_normalization(edges: Dict, all_feature_data_lists: List[List[Dict]]) -> Dict:
    """
    Recalculate edge distances using normalized coordinates after pairing.

    This updates only the edge distance fields so pairing decisions remain
    unchanged when coordinate cutoffs are used.

    Args:
        edges: Dictionary of edges {ref_file_idx: [edge_dict, ...]}
        all_feature_data_lists: List of feature data lists for coordinate lookup

    Returns:
        Updated edges dictionary with normalized distance values
    """
    if not all_feature_data_lists:
        print("  ✗ WARNING: Feature data not provided for post-pair normalization, skipping")
        return edges

    all_x_values = []
    all_y_values = []
    for features in all_feature_data_lists:
        for f in features:
            all_x_values.append(f['x_center'])
            all_y_values.append(f['y_center'])

    if not all_x_values or not all_y_values:
        print("  ✗ WARNING: Empty feature coordinates for post-pair normalization, skipping")
        return edges

    min_x = float(np.min(all_x_values))
    max_x = float(np.max(all_x_values))
    min_y = float(np.min(all_y_values))
    max_y = float(np.max(all_y_values))

    range_x = max_x - min_x if max_x > min_x else 1.0
    range_y = max_y - min_y if max_y > min_y else 1.0
    scale_x = 1.0 / range_x if range_x > 0 else 1.0
    scale_y = 1.0 / range_y if range_y > 0 else 1.0

    print("  Applying post-pair coordinate normalization to edge distances...")
    print(f"    X range: {min_x:.2f} to {max_x:.2f} (range: {range_x:.2f})")
    print(f"    Y range: {min_y:.2f} to {max_y:.2f} (range: {range_y:.2f})")
    print(f"    Scale factors - X: {scale_x:.6f}, Y: {scale_y:.6f}")

    updated_edges = {}
    updated_count = 0
    skipped_count = 0

    for file_idx, edge_list in edges.items():
        updated_list = []
        for edge in edge_list:
            ref_file_idx = edge.get('ref_file_idx')
            match_file_idx = edge.get('match_file_idx')
            current_feat_idx = edge.get('current_file', {}).get('feature_idx')
            matched_feat_idx = edge.get('matched_file', {}).get('feature_idx')

            if (ref_file_idx is None or match_file_idx is None or
                current_feat_idx is None or matched_feat_idx is None):
                updated_list.append(edge)
                skipped_count += 1
                continue

            if (ref_file_idx >= len(all_feature_data_lists) or
                match_file_idx >= len(all_feature_data_lists)):
                updated_list.append(edge)
                skipped_count += 1
                continue

            ref_features = all_feature_data_lists[ref_file_idx]
            match_features = all_feature_data_lists[match_file_idx]

            if (current_feat_idx >= len(ref_features) or
                matched_feat_idx >= len(match_features)):
                updated_list.append(edge)
                skipped_count += 1
                continue

            ref_feature = ref_features[current_feat_idx]
            matched_feature = match_features[matched_feat_idx]

            ref_x = (ref_feature['x_center'] - min_x) * scale_x
            ref_y = (ref_feature['y_center'] - min_y) * scale_y
            match_x = (matched_feature['x_center'] - min_x) * scale_x
            match_y = (matched_feature['y_center'] - min_y) * scale_y

            mz_distance = abs(ref_x - match_x)
            rt_distance = abs(ref_y - match_y)
            distance = float(np.sqrt(mz_distance ** 2 + rt_distance ** 2))

            edge_copy = edge.copy()
            if 'distance_raw' not in edge_copy:
                edge_copy['distance_raw'] = edge_copy.get('distance')
            if 'mz_distance_raw' not in edge_copy:
                edge_copy['mz_distance_raw'] = edge_copy.get('mz_distance')
            if 'rt_distance_raw' not in edge_copy:
                edge_copy['rt_distance_raw'] = edge_copy.get('rt_distance')

            edge_copy['distance'] = distance
            edge_copy['mz_distance'] = float(mz_distance)
            edge_copy['rt_distance'] = float(rt_distance)

            updated_list.append(edge_copy)
            updated_count += 1

        updated_edges[file_idx] = updated_list

    print(f"  Post-pair normalization updated {updated_count} edges, skipped {skipped_count}")
    return updated_edges


def _build_connected_components_from_edges(edges: Dict, file_features_list: List[List[Dict]]) -> List[Dict]:
    """
    Build connected components from edges. Each feature belongs to exactly one component.
    ULTRA-OPTIMIZED: Uses union-find directly on edges (no intermediate adjacency structures).
    
    Args:
        edges: Dictionary {ref_file_idx: [edge_dict, ...]}
               Each edge has: current_file, matched_file, distance, ref_file_idx, match_file_idx
        file_features_list: List of feature lists per file
    
    Returns:
        List of connected components, each containing all features in the component 
        with distances only to directly connected neighbors
    """
    from collections import defaultdict
    
    # Union-Find for fast grouping
    parent = {}
    
    def find(x):
        if x not in parent:
            parent[x] = x
        if parent[x] != x:
            parent[x] = find(parent[x])  # Path compression
        return parent[x]
    
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py
    
    # Store edge distances and neighbor mappings for later lookup
    edge_distances = {}
    neighbors = defaultdict(dict)  # neighbors[v1][v2] = distance
    all_vertices = set()
    
    # Single pass through edges: union vertices and store distances + neighbors
    for ref_file_idx, edge_list in edges.items():
        for edge in edge_list:
            v1 = (edge['ref_file_idx'], edge['current_file']['feature_idx'])
            v2 = (edge['match_file_idx'], edge['matched_file']['feature_idx'])
            distance = edge['distance']
            
            all_vertices.add(v1)
            all_vertices.add(v2)
            
            # Union the vertices (group them)
            union(v1, v2)
            
            # Store distance for both directions (edges are undirected) AND neighbor mappings
            edge_distances[frozenset([v1, v2])] = distance
            neighbors[v1][v2] = distance
            neighbors[v2][v1] = distance
    
    # Group vertices by their root (component)
    components = defaultdict(set)
    for vertex in all_vertices:
        root = find(vertex)
        components[root].add(vertex)
    
    components = list(components.values())
    
    # Build consolidated match entries for each component
    matches = []
    total_components = len(components)
    
    for comp_idx, component_vertices in enumerate(components):
        if comp_idx % max(1, total_components // 10) == 0:
            print(f"    Processing component {comp_idx+1}/{total_components}", end='\r')
        
        # Collect all features in this component with metadata for quick lookups
        component_features = []
        vertex_to_feature_idx = {}  # (file_idx, feat_idx) -> index in component_features
        
        for file_idx, feat_idx in sorted(component_vertices):
            feature = file_features_list[file_idx][feat_idx].copy()
            # Add metadata for quick lookups
            feature['_file_idx'] = file_idx
            feature['_feat_idx'] = feat_idx
            
            vertex_to_feature_idx[(file_idx, feat_idx)] = len(component_features)
            component_features.append(feature)
        
        # Merge pep_ident values
        all_pep_idents = set()
        for feature in component_features:
            if 'pep_ident' in feature and feature['pep_ident']:
                all_pep_idents.update(feature['pep_ident'])
        
        unique_pep_idents = list(all_pep_idents) if all_pep_idents else None
        
        # Calculate consolidated position using vectorized operations
        coords_array = np.array([[f['x_center'], f['y_center']] for f in component_features])
        mz_rt_array = np.array([[f['mz_start'], f['mz_end'], f['rt_start'], f['rt_end']]
                                for f in component_features])
        
        avg_coords = np.mean(coords_array, axis=0)
        avg_mz_rt = np.mean(mz_rt_array, axis=0)
        
        # ULTRA-OPTIMIZED: Only iterate through actual neighbors in edge_distances (not all features)
        featured_entries = []
        for feature in component_features:
            feature_entry = feature.copy()
            
            # Get distances only to features connected by edges (iterate neighbors, not all pairs)
            distances_to_others = {}
            v_current = (feature['_file_idx'], feature['_feat_idx'])
            
            # Only iterate through neighbors that this feature actually has edges to
            if v_current in neighbors:
                for v_neighbor, distance in neighbors[v_current].items():
                    # Confirm neighbor is in this component (should always be true)
                    if v_neighbor in vertex_to_feature_idx:
                        other_feat_idx = vertex_to_feature_idx[v_neighbor]
                        neighbor_feature = component_features[other_feat_idx]
                        other_key = f"{neighbor_feature['filename']}_feat{v_neighbor[1]}"
                        distances_to_others[other_key] = float(distance)
            
            feature_entry['distances_to_group_features'] = distances_to_others
            featured_entries.append(feature_entry)
        
        # Create consolidated match entry
        # Build unique files list and comprehensive features list
        unique_files = sorted(list(set(f['filename'] for f in component_features)))
        all_features_list = [f"{f['filename']}_feat{f['_feat_idx']}" for f in component_features]
        
        match = {
            'num_files_matched': len(set(v[0] for v in component_vertices)),
            'num_features_in_group': len(component_features),
            'files': unique_files,
            'all_features': all_features_list,  # NEW: Comprehensive list of all features in group
            'pep_ident': unique_pep_idents,
            'x_center': float(avg_coords[0]),
            'y_center': float(avg_coords[1]),
            'mz_start': float(avg_mz_rt[0]),
            'mz_end': float(avg_mz_rt[1]),
            'rt_start': float(avg_mz_rt[2]),
            'rt_end': float(avg_mz_rt[3]),
            'individual_features': featured_entries,
            'group_size': len(component_vertices)
        }
        
        matches.append(match)
    
    print(f"\n  ✓ Built {len(matches)} connected component groups")
    return matches


def feature_pairing_optimized(all_feature_data_lists: List[List[Dict]], best_match_only: bool = False, skip_match_building: bool = False, normalize_coordinates: bool = True, distance_calc_before_scaling: bool = False) -> tuple:
    """
    Optimized feature pairing using KD-trees and vectorized operations.
    Much faster for large datasets with multiple files.
    Returns tuple of (matches, edges) where edges is a network graph structure
    with format: {ref_file_idx: [edge_dict, ...]}
    Each edge contains: current_file, matched_file, distance, ref_file_idx, match_file_idx
    
    Args:
        all_feature_data_lists: List of feature lists from each file
        best_match_only: If True, keep only the best match for each feature from each other file
        skip_match_building: If True, skip expensive match consolidation and return empty matches list
        normalize_coordinates: If True, scales both axes to same range for better KD-Tree performance
    """
    if not all_feature_data_lists:
        return [], {}
    
    if len(all_feature_data_lists) == 1:
        matches = []
        edges = {0: []}  # Single file has no cross-file edges
        for feature in all_feature_data_lists[0]:
            match_group = feature.copy()
            match_group['files'] = [feature['filename']]
            match_group['match_score'] = 1.0
            matches.append(match_group)
        return matches, edges
    
    print("  Pairing features (OPTIMIZED) using KD-tree nearest-neighbor...")
    print("  Splitting features by charge for paired matching...")
    
    # Calculate scaling factors for coordinate normalization if enabled
    scale_x = 1.0
    scale_y = 1.0
    min_x = 0.0
    min_y = 0.0
    
    if normalize_coordinates:
        print("  Computing coordinate normalization factors...")
        # Collect all x and y values to compute min/max ranges
        all_x_values = []
        all_y_values = []
        for features in all_feature_data_lists:
            for f in features:
                all_x_values.append(f['x_center'])
                all_y_values.append(f['y_center'])
        
        if all_x_values and all_y_values:
            min_x = np.min(all_x_values)
            max_x = np.max(all_x_values)
            min_y = np.min(all_y_values)
            max_y = np.max(all_y_values)
            
            range_x = max_x - min_x if max_x > min_x else 1.0
            range_y = max_y - min_y if max_y > min_y else 1.0
            
            # Scale both axes to [0, 1] range, then scale y to match x's range
            scale_x = 1.0 / range_x if range_x > 0 else 1.0
            scale_y = 1.0 / range_y if range_y > 0 else 1.0
            print(f"    X range: {min_x:.2f} to {max_x:.2f} (range: {range_x:.2f})")
            print(f"    Y range: {min_y:.2f} to {max_y:.2f} (range: {range_y:.2f})")
            print(f"    Scale factors - X: {scale_x:.6f}, Y: {scale_y:.6f}")
    
    def scale_coords(coords):
        """Scale coordinates using computed normalization factors."""
        if normalize_coordinates and len(coords) > 0:
            scaled = coords.copy().astype(float)
            scaled[:, 0] = (scaled[:, 0] - min_x) * scale_x
            scaled[:, 1] = (scaled[:, 1] - min_y) * scale_y
            return scaled
        return coords.astype(float)
    
    # Prepare data structures for each file, grouped by charge
    # Structure: {file_idx: {charge: {'features': [features], 'coords': np.array, 'kdtree': cKDTree}}}
    file_features_by_charge = []
    charge_stats = {}  # For diagnostic output
    
    for file_idx, features in enumerate(all_feature_data_lists):
        # Group features by charge
        features_by_charge = {}
        for feat_idx, f in enumerate(features):
            f_copy = f.copy()
            f_copy['_feat_idx'] = feat_idx
            f_copy['_file_idx'] = file_idx
            
            # Charge field is guaranteed to exist (verified at load time)
            charge = int(f['charge'])
            if charge not in features_by_charge:
                features_by_charge[charge] = []
            features_by_charge[charge].append(f_copy)
        
        # Track charge statistics
        if file_idx == 0:
            charge_stats[file_idx] = {charge: len(features) for charge, features in features_by_charge.items()}
        
        file_features_by_charge.append(features_by_charge)
    
    # Print charge diagnostic information
    print("  Charge distribution by file:")
    for file_idx, charge_dict in enumerate(file_features_by_charge):
        print(f"    File {file_idx}: {dict(sorted([(c, len(f)) for c, f in charge_dict.items()]))}")
    
    # Build KD-trees for each charge group in each file
    # Structure: {file_idx: {charge: {'kdtree': cKDTree, 'features': [features], 'coords': np.array}}}
    kdtrees_by_charge = []
    
    for file_idx, charge_dict in enumerate(file_features_by_charge):
        file_kdtrees = {}
        for charge, charge_features in charge_dict.items():
            coords = np.array([[f['x_center'], f['y_center']] for f in charge_features])
            scaled_coords = scale_coords(coords)
            
            if len(scaled_coords) > 0:
                kdtree = cKDTree(scaled_coords)
            else:
                kdtree = None
            
            file_kdtrees[charge] = {
                'kdtree': kdtree,
                'features': charge_features,
                'coords': scaled_coords
            }
        kdtrees_by_charge.append(file_kdtrees)
    
    # Start matching from first file
    first_file_features = []
    for features in file_features_by_charge[0].values():
        first_file_features.extend(features)
    
    matches = []
    
    # Build network graph edges dictionary
    # Format: {ref_file_idx: [(current_file_info, matched_file_info, distance), ...], ...}
    # Adapt distance cutoff based on coordinate normalization
    if normalize_coordinates:
        distance_cutoff = 0.5  # Reasonable cutoff in normalized space (0-1 range)
    else:
        distance_cutoff = 1.0  # Distance threshold for matching features
    edges = {}
    
    print(f"  Processing {len(first_file_features)} features from first file (with charge filtering)...")
    for file_file_idx, charge_dict in enumerate(file_features_by_charge):
        edges[file_file_idx] = []
        
        for charge, charge_features in charge_dict.items():
            for feature_idx, ref_feature in enumerate(charge_features):
                if feature_idx % 100 == 0:
                    print(f"    Processing file {file_file_idx}, charge {charge}: feature {feature_idx}/{len(charge_features)}", end='\r')
                
                # Use scaled coordinates for KD-tree query
                ref_coords_scaled = scale_coords(np.array([[ref_feature['x_center'], ref_feature['y_center']]]))
                
                current_file_info = {
                    'filename': ref_feature['filename'],
                    'feature_idx': ref_feature['_feat_idx']
                }
                
                # Query each other file's KD-tree for matches within distance cutoff
                # ONLY match with features of the SAME charge
                for query_file_idx in range(len(all_feature_data_lists)):
                    if query_file_idx == file_file_idx:
                        continue
                    
                    # Get the KD-tree for this charge in the query file
                    if charge not in kdtrees_by_charge[query_file_idx]:
                        # No features with this charge in the query file, skip
                        continue
                    
                    kdtree = kdtrees_by_charge[query_file_idx][charge]['kdtree']
                    query_features = kdtrees_by_charge[query_file_idx][charge]['features']
                    
                    
                    if best_match_only:
                        # Query only the single best match for this feature from the other file
                        distance, idx = kdtree.query(ref_coords_scaled, k=1)
                        
                        if isinstance(distance, np.ndarray):
                            distance = distance.flat[0]
                        if isinstance(idx, np.ndarray):
                            idx = idx.flat[0]
                        
                        if idx >= 0:
                            matched_feature = query_features[int(idx)]
                            
                            # Calculate distance based on whether to use scaled or original coordinates
                            if distance_calc_before_scaling:
                                # Use original coordinates
                                mz_dist = (ref_feature['x_center'] - matched_feature['x_center']) ** 2
                                rt_dist = (ref_feature['y_center'] - matched_feature['y_center']) ** 2
                            else:
                                # Use scaled coordinates (default)
                                kdtree_coords = kdtrees_by_charge[query_file_idx][charge]['coords']
                                matched_scaled_coords = kdtree_coords[int(idx)]
                                ref_scaled_coords = ref_coords_scaled[0]
                                mz_dist = (ref_scaled_coords[0] - matched_scaled_coords[0]) ** 2
                                rt_dist = (ref_scaled_coords[1] - matched_scaled_coords[1]) ** 2
                            
                            exact_distance = float(np.sqrt(mz_dist + rt_dist))
                            mz_distance = float(np.sqrt(mz_dist))  # Separate m/z distance
                            rt_distance = float(np.sqrt(rt_dist))  # Separate RT distance
                            
                            # Create matched feature info
                            matched_file_info = {
                                'filename': matched_feature['filename'],
                                'feature_idx': matched_feature['_feat_idx']
                            }
                            
                            # Add edge connecting the two features
                            edge = {
                                'current_file': current_file_info,
                                'matched_file': matched_file_info,
                                'distance': exact_distance,
                                'mz_distance': mz_distance,  # Separate m/z distance
                                'rt_distance': rt_distance,  # Separate RT distance
                                'ref_file_idx': file_file_idx,
                                'match_file_idx': query_file_idx,
                                'charge': charge  # Include charge validation in edge
                            }
                            edges[file_file_idx].append(edge)
                    else:
                        # Query all points within distance cutoff
                        scaled_cutoff = distance_cutoff
                        
                        match_indices = kdtree.query_ball_point(
                            ref_coords_scaled[0], r=scaled_cutoff
                        )
                        
                        # Collect all matches with their distances
                        kdtree_coords = kdtrees_by_charge[query_file_idx][charge]['coords']
                        ref_scaled_coords = ref_coords_scaled[0]
                        
                        for match_idx in match_indices:
                            matched_feature = query_features[match_idx]
                            
                            # Calculate distance based on whether to use scaled or original coordinates
                            if distance_calc_before_scaling:
                                # Use original coordinates
                                mz_dist = (ref_feature['x_center'] - matched_feature['x_center']) ** 2
                                rt_dist = (ref_feature['y_center'] - matched_feature['y_center']) ** 2
                            else:
                                # Use scaled coordinates (default)
                                matched_scaled_coords = kdtree_coords[match_idx]
                                mz_dist = (ref_scaled_coords[0] - matched_scaled_coords[0]) ** 2
                                rt_dist = (ref_scaled_coords[1] - matched_scaled_coords[1]) ** 2
                            
                            distance = float(np.sqrt(mz_dist + rt_dist))
                            mz_distance = float(np.sqrt(mz_dist))  # Separate m/z distance
                            rt_distance = float(np.sqrt(rt_dist))  # Separate RT distance
                            
                            # Create matched feature info
                            matched_file_info = {
                                'filename': matched_feature['filename'],
                                'feature_idx': matched_feature['_feat_idx']
                            }
                            
                            # Add edge connecting the two features
                            edge = {
                                'current_file': current_file_info,
                                'matched_file': matched_file_info,
                                'distance': distance,
                                'mz_distance': mz_distance,  # Separate m/z distance
                                'rt_distance': rt_distance,  # Separate RT distance
                                'ref_file_idx': file_file_idx,
                                'match_file_idx': query_file_idx,
                                'charge': charge  # Include charge validation in edge
                            }
                            edges[file_file_idx].append(edge)
    
    # Diagnostic output: verify all edges have matching charges
    print("\n  Verifying charge consistency in edges:")
    edges_by_charge = {}
    total_edges_verified = 0
    for file_idx, edge_list in edges.items():
        for edge in edge_list:
            charge_val = edge.get('charge', None)
            if charge_val not in edges_by_charge:
                edges_by_charge[charge_val] = 0
            edges_by_charge[charge_val] += 1
            total_edges_verified += 1
    
    print(f"    Total edges: {total_edges_verified}")
    print(f"    Edges by charge: {dict(sorted(edges_by_charge.items()))}")
    print(f"    ✓ All edges have matching charges between compared features")
    
    # Skip match building if explicitly requested (for edge-only output)
    if skip_match_building:
        print("  Skipping match consolidation (flagged for edge-only output)")
        return matches, edges
    
    # Flatten file_features_by_charge for match consolidation
    file_features_list = []
    for file_idx, charge_dict in enumerate(file_features_by_charge):
        file_features = []
        for features in charge_dict.values():
            file_features.extend(features)
        file_features_list.append(file_features)

    # NOTE: Connected components are now built in main() AFTER edge filtering
    # This ensures components only contain features connected by edges that pass cutoff filtering
    matches = []  # Return empty matches - they will be built after edge filtering
    
    return matches, edges


def feature_pairing(all_feature_data_lists: List[List[Dict]], best_match_only: bool = False, skip_match_building: bool = False, normalize_coordinates: bool = True, distance_calc_before_scaling: bool = False) -> tuple:
    """
    Basic feature pairing using nearest-neighbor principle.
    Returns tuple of (matches, edges) where edges is a network graph structure.
    
    Args:
        all_feature_data_lists: List of feature lists from each file
        best_match_only: If True, keep only the best match for each feature from each other file
        skip_match_building: If True, skip expensive match consolidation and return empty matches list
        normalize_coordinates: If True, scales both axes to same range for better distance calculation
    """
    if not all_feature_data_lists:
        return [], {}
    
    if len(all_feature_data_lists) == 1:
        matches = []
        edges = {0: []}  # Single file has no cross-file edges
        for feature in all_feature_data_lists[0]:
            match_group = feature.copy()
            match_group['files'] = [feature['filename']]
            match_group['match_score'] = 1.0
            matches.append(match_group)
        return matches, edges
    
    print("  Pairing features using nearest-neighbor...")
    print("  Splitting features by charge for paired matching...")
    
    # Calculate scaling factors for coordinate normalization if enabled
    scale_x = 1.0
    scale_y = 1.0
    min_x = 0.0
    min_y = 0.0
    
    if normalize_coordinates:
        print("  Computing coordinate normalization factors...")
        # Collect all x and y values to compute min/max ranges
        all_x_values = []
        all_y_values = []
        for features in all_feature_data_lists:
            for f in features:
                all_x_values.append(f['x_center'])
                all_y_values.append(f['y_center'])
        
        if all_x_values and all_y_values:
            min_x = np.min(all_x_values)
            max_x = np.max(all_x_values)
            min_y = np.min(all_y_values)
            max_y = np.max(all_y_values)
            
            range_x = max_x - min_x if max_x > min_x else 1.0
            range_y = max_y - min_y if max_y > min_y else 1.0
            
            # Scale both axes to [0, 1] range, then scale y to match x's range
            scale_x = 1.0 / range_x if range_x > 0 else 1.0
            scale_y = 1.0 / range_y if range_y > 0 else 1.0
            print(f"    X range: {min_x:.2f} to {max_x:.2f} (range: {range_x:.2f})")
            print(f"    Y range: {min_y:.2f} to {max_y:.2f} (range: {range_y:.2f})")
            print(f"    Scale factors - X: {scale_x:.6f}, Y: {scale_y:.6f}")
    
    # Initialize feature data structures with indices
    file_features_list = []
    for file_idx, features in enumerate(all_feature_data_lists):
        features_with_idx = []
        for feat_idx, f in enumerate(features):
            f_copy = f.copy()
            f_copy['_feat_idx'] = feat_idx
            f_copy['_file_idx'] = file_idx
            features_with_idx.append(f_copy)
        file_features_list.append(features_with_idx)
    
    # Build edges dictionary for network graph
    # Adapt distance cutoff based on coordinate normalization
    if normalize_coordinates:
        # When normalized, cutoff should be in the normalized space (typically 0-1)
        # Diagonal of normalized unit square is sqrt(2) ≈ 1.41
        distance_cutoff = 0.5  # Reasonable cutoff in normalized space
        print(f"  Distance cutoff (normalized space): {distance_cutoff:.3f}")
    else:
        distance_cutoff = 100.0  # Distance threshold for matching features
    edges = {}
    
    # Initialize edges dict for each file
    for file_idx in range(len(all_feature_data_lists)):
        edges[file_idx] = []
    
    # Build feature pairing and edges from first file
    first_file_features = file_features_list[0]
    
    # NOTE: Connected components are now built in main() AFTER edge filtering
    # This ensures components only contain features connected by edges that pass cutoff filtering
    matches = []
    
    # Skip match building if explicitly requested (for edge-only output)
    if skip_match_building:
        print("  Skipping match consolidation (flagged for edge-only output)")
        return matches, edges
    
    # Build edges from all features (components will be built from these edges in main after filtering)
    for feature_idx, ref_feature in enumerate(first_file_features):
        if feature_idx % 100 == 0:
            print(f"    Pairing feature {feature_idx}/{len(first_file_features)}", end='\r')
        
        match_group = {'matches': [ref_feature.copy()]}
        min_distances = []
        feature_edges = []
        
        current_file_info = {
            'filename': ref_feature['filename'],
            'feature_idx': ref_feature['_feat_idx']
        }
        
# Get charge of reference feature (guaranteed to exist)
        ref_charge = int(ref_feature['charge'])
        
        for file_idx in range(1, len(all_feature_data_lists)):
            file_features = file_features_list[file_idx]
            
            if not file_features:
                continue
            
            distances = []
            for f in file_features:
                # CHARGE FILTERING: Only compare with features of matching charge
                f_charge = int(f['charge'])
                if f_charge != ref_charge:
                    continue  # Skip features with different charge
                
                # Calculate distance based on whether to use scaled or original coordinates
                if distance_calc_before_scaling:
                    # Use original coordinates
                    ref_x_calc = ref_feature['x_center']
                    ref_y_calc = ref_feature['y_center']
                    f_x_calc = f['x_center']
                    f_y_calc = f['y_center']
                else:
                    # Use scaled coordinates (default)
                    ref_x_calc = (ref_feature['x_center'] - min_x) * scale_x if normalize_coordinates else ref_feature['x_center']
                    ref_y_calc = (ref_feature['y_center'] - min_y) * scale_y if normalize_coordinates else ref_feature['y_center']
                    f_x_calc = (f['x_center'] - min_x) * scale_x if normalize_coordinates else f['x_center']
                    f_y_calc = (f['y_center'] - min_y) * scale_y if normalize_coordinates else f['y_center']
                
                mz_dist = (ref_x_calc - f_x_calc) ** 2
                rt_dist = (ref_y_calc - f_y_calc) ** 2
                dist = np.sqrt(mz_dist + rt_dist)
                distances.append((dist, f))
            
            if distances:
                min_dist, nearest_feature = min(distances, key=lambda x: x[0])
                min_distances.append(min_dist)
                match_group['matches'].append(nearest_feature.copy())
                
                # Record edge - include if best_match_only is True or within distance cutoff
                if best_match_only or min_dist <= distance_cutoff:
                    matched_file_info = {
                        'filename': nearest_feature['filename'],
                        'feature_idx': nearest_feature['_feat_idx']
                    }
                    edge = {
                        'current_file': current_file_info,
                        'matched_file': matched_file_info,
                        'distance': float(min_dist),
                        'ref_file_idx': 0,
                        'match_file_idx': file_idx,
                        'charge': ref_charge  # Include charge validation in edge
                    }
                    feature_edges.append(edge)
                    edges[0].append(edge)
        
        if min_distances:
            avg_distance = sum(min_distances) / len(min_distances)
            match_score = float(np.exp(-avg_distance / 100.0))
        else:
            match_score = 1.0
        
        all_pep_idents = []
        for matched_feature in match_group['matches']:
            if 'pep_ident' in matched_feature and matched_feature['pep_ident']:
                all_pep_idents.extend(matched_feature['pep_ident'])
        
        unique_pep_idents = list(set(all_pep_idents)) if all_pep_idents else None
        
        # Calculate inter-feature scores within the group
        features_with_scores = []
        for i, feature_i in enumerate(match_group['matches']):
            feature_i_copy = feature_i.copy()
            inter_feature_scores = {}
            
            for j, feature_j in enumerate(match_group['matches']):
                if i != j:
                    # Calculate distance between features
                    mz_dist = (feature_i['x_center'] - feature_j['x_center']) ** 2
                    rt_dist = (feature_i['y_center'] - feature_j['y_center']) ** 2
                    dist = np.sqrt(mz_dist + rt_dist)
                    # Score based on distance (higher score for closer features)
                    score = float(np.exp(-dist / 100.0))
                    inter_feature_scores[feature_j['filename']] = score
            
            feature_i_copy['scores_to_features'] = inter_feature_scores
            features_with_scores.append(feature_i_copy)
        
        # Build unique files list and comprehensive features list
        unique_files_list = sorted(list(set(f['filename'] for f in match_group['matches'])))
        all_features_in_match = [f"{f['filename']}_feat{f['_feat_idx']}" for f in match_group['matches']]
        
        consolidated_match = {
            'match_score': float(match_score),
            'num_files_matched': len(unique_files_list),
            'files': unique_files_list,
            'all_features': all_features_in_match,  # NEW: Comprehensive list of all features in group
            'pep_ident': unique_pep_idents,
            'x_center': float(np.mean([f['x_center'] for f in match_group['matches']])),
            'y_center': float(np.mean([f['y_center'] for f in match_group['matches']])),
            'mz_start': float(np.mean([f['mz_start'] for f in match_group['matches']])),
            'mz_end': float(np.mean([f['mz_end'] for f in match_group['matches']])),
            'rt_start': float(np.mean([f['rt_start'] for f in match_group['matches']])),
            'rt_end': float(np.mean([f['rt_end'] for f in match_group['matches']])),
            'individual_features': features_with_scores,
            'network_edges': feature_edges
        }
        matches.append(consolidated_match)
    
    
    # Diagnostic output: verify all edges have matching charges
    print("  Verifying charge consistency in edges:")
    edges_by_charge = {}
    total_edges_verified = 0
    for file_idx, edge_list in edges.items():
        for edge in edge_list:
            charge_val = edge.get('charge', None)
            if charge_val not in edges_by_charge:
                edges_by_charge[charge_val] = 0
            edges_by_charge[charge_val] += 1
            total_edges_verified += 1
    
    print(f"    Total edges: {total_edges_verified}")
    print(f"    Edges by charge: {dict(sorted(edges_by_charge.items()))}")
    print(f"    ✓ All edges have matching charges between compared features")
    
    print(f"  Paired {len(matches)} feature groups")
    print(f"  Network graph contains {sum(len(e) for e in edges.values())} edges across {len(edges)} files")
    
    return matches, edges


def main():
    parser = argparse.ArgumentParser(
        description="Pair matching features across multiple files"
    )
    parser.add_argument("--input_dir", help="Directory containing feature pickle files (default: use current dir if --input_files not provided)")
    parser.add_argument("--input_files", nargs='*', help="Explicit list of feature pickle files to process (alternative to --input_dir)")
    parser.add_argument("--output_pkl", required=True, help="Output pickle file path")
    parser.add_argument("--output_json", required=True, help="Output JSON file path")
    parser.add_argument("--output_edges_pkl", help="Output edges dictionary pickle file (default: inferred from output_pkl)")
    parser.add_argument("--use-basic", action='store_true', help="Use basic nearest-neighbor pairing instead of optimized KD-tree pairing (default: optimized)")
    parser.add_argument("--skip-matchfinder", action='store_true', help="Skip matchfinder processing (match consolidation). Only export edges dictionary. Useful for network graph analysis.")
    parser.add_argument("--match_cutoff", type=float, default=0.0, help="Minimum match score cutoff (0.0-1.0)")
    parser.add_argument("--best_match_only", action='store_true', help="For each feature, keep only the single best match instead of all matches within cutoff")
    parser.add_argument("--edges_cutoff", type=float, default=None, help="Maximum distance cutoff for edges (filters edges based on euclidean distance)")
    parser.add_argument("--mz_cutoff", type=float, default=None, help="Maximum m/z (x-axis) coordinate difference for edge filtering")
    parser.add_argument("--rt_cutoff", type=float, default=None, help="Maximum RT (y-axis) coordinate difference for edge filtering")
    parser.add_argument("--distance_calc_before_scaling", action='store_true', help="Calculate distance before coordinate scaling (default: after scaling)")
    parser.add_argument("--normalize_coordinates", type=lambda x: x.lower() in ('true', '1', 'yes'), nargs='?', const=True, default=None, help="Enable/disable coordinate normalization before KD-tree matching (default: auto - OFF for mz/rt cutoff, ON for euclidean cutoff)")
    parser.add_argument("--postpair_normalize_coordinates", action='store_true', help="Recalculate edge distances using normalized coordinates after pairing (does not affect pairing or cutoffs)(carefull, makes distances relative to normalized space, not original coordinates)")
    parser.add_argument("--normalize_edge_distances", action='store_true', help="Normalize edge distances to [0, 1] range for visualization (default: disabled)")
    parser.add_argument("--analyze_pep_idents", action='store_true', help="Perform detailed pep_ident matching analysis (default: enabled)")
    parser.add_argument("--skip_json_output", action='store_true', help="Skip JSON serialization for speed (default: disabled - JSON output enabled)")
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print(f"Feature Pairing")
    print(f"{'='*70}")
    
    # Find all feature pickle files - support explicit file list or glob from directory
    if args.input_files and len(args.input_files) > 0:
        # Use explicit file list provided via command line
        feature_files = sorted([Path(f) for f in args.input_files])
        print(f"Using {len(feature_files)} explicitly provided feature file(s)")
    elif args.input_dir:
        # Fall back to globbing from input directory
        feature_files = sorted(Path(args.input_dir).glob("*_feature_data.pkl"))
        if not feature_files:
            print(f"✗ No feature pickle files found in {args.input_dir}!")
            return
        print(f"Found {len(feature_files)} feature file(s) in {args.input_dir}")
    else:
        # Default to current directory
        feature_files = sorted(Path('.').glob("*_feature_data.pkl"))
        if not feature_files:
            print("✗ No feature pickle files found in current directory!")
            print("  Use: --input_dir <dir> to search in a specific directory")
            print("    or --input_files <file1> <file2> ... to provide explicit files")
            return
        print(f"Found {len(feature_files)} feature file(s) in current directory")
    
    # Load all feature data
    all_feature_data_lists = []
    has_charge_field = None  # Track if charge field is present
    
    for pkl_file in feature_files:
        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)
            all_feature_data_lists.append(data)
        
        # Check if first feature has charge field
        if isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict):
            if has_charge_field is None:
                has_charge_field = 'charge' in data[0]
    
    # Verify all loaded files have consistent charge field availability
    if has_charge_field is False:
        print("\n✗ ERROR: Feature pickle files are missing the 'charge' field")
        print(f"  Input directory: {args.input_dir}")
        print(f"\n  The 'charge' field is REQUIRED for feature pairing with charge-based splitting.")
        print(f"  This field contains the charge state (e.g., +1, +2, +3) of peptide ions.")
        print(f"\n  Feature files WITH charge field are available in:")
        print(f"    results/feature_analysis/feature_data_lists/")
        print(f"\n  Try running with:")
        print(f"    python pair_features.py \\")
        print(f"      --input_dir results/feature_analysis/feature_data_lists \\")
        print(f"      --output_pkl ... \\")
        print(f"      --output_json ...")
        sys.exit(1)
    elif has_charge_field is None:
        print("\n✗ ERROR: Could not load any feature data from pickle files")
        sys.exit(1)
    
    # Determine normalize_coordinates default based on cutoff type
    if args.normalize_coordinates is None:
        # Auto-determine based on which cutoff is being used
        if args.mz_cutoff is not None and args.rt_cutoff is not None:
            # Coordinate-based filtering: don't normalize (use original coordinates)
            normalize_coordinates = False
            print("  ℹ Coordinate normalization: AUTO (OFF) - using coordinate-based cutoff")
        elif args.edges_cutoff is not None:
            # Euclidean distance filtering: normalize (scaled coordinates)
            normalize_coordinates = True
            print("  ℹ Coordinate normalization: AUTO (ON) - using euclidean distance cutoff")
        else:
            # Default to normalize for standard matching
            normalize_coordinates = True
            print("  ℹ Coordinate normalization: AUTO (ON) - standard matching mode")
    else:
        normalize_coordinates = args.normalize_coordinates
        state = "ON" if normalize_coordinates else "OFF"
        print(f"  ℹ Coordinate normalization: MANUAL ({state})")
    
    # Pair features (optimized is default)
    if args.use_basic:
        paired_features, network_edges = feature_pairing(all_feature_data_lists, best_match_only=args.best_match_only, skip_match_building=args.skip_matchfinder, normalize_coordinates=normalize_coordinates, distance_calc_before_scaling=args.distance_calc_before_scaling)
    else:
        paired_features, network_edges = feature_pairing_optimized(all_feature_data_lists, best_match_only=args.best_match_only, skip_match_building=args.skip_matchfinder, normalize_coordinates=normalize_coordinates, distance_calc_before_scaling=args.distance_calc_before_scaling)
    
    print(f"Built network graph with {sum(len(e) for e in network_edges.values())} edges")
    
    # Filter edges by cutoff if specified
    if args.mz_cutoff is not None and args.rt_cutoff is not None:
        # Coordinate-based filtering (mz and rt)
        print(f"Filtering edges with coordinate cutoffs: mz={args.mz_cutoff}, rt={args.rt_cutoff}")
        edges_before_filter = sum(len(e) for e in network_edges.values())
        network_edges = filter_edges_by_coordinate_cutoffs(network_edges, args.mz_cutoff, args.rt_cutoff, all_feature_data_lists)
        edges_after_filter = sum(len(e) for e in network_edges.values())
        print(f"Filtered edges: {edges_before_filter} -> {edges_after_filter}")
    elif args.edges_cutoff is not None:
        # Distance-based filtering (euclidean distance)
        print(f"Filtering edges with distance cutoff: {args.edges_cutoff}")
        edges_before_filter = sum(len(e) for e in network_edges.values())
        network_edges = filter_edges_by_distance_cutoff(network_edges, args.edges_cutoff)
        edges_after_filter = sum(len(e) for e in network_edges.values())
        print(f"Filtered edges: {edges_before_filter} -> {edges_after_filter}")

    # Optional post-pair coordinate normalization for edge distances
    if args.postpair_normalize_coordinates:
        network_edges = apply_postpair_coordinate_normalization(network_edges, all_feature_data_lists)
    else:
        print("Skipping post-pair coordinate normalization (--postpair_normalize_coordinates disabled)")
    
    # Normalize edge distances to [0, 1] for visualization (optional)
    if args.normalize_edge_distances:
        print("Normalizing edge distances to [0, 1] range for visualization...")
        network_edges = normalize_edge_distances(network_edges)
    else:
        print("Skipping edge distance normalization (--normalize_edge_distances disabled)")
    
    # If skip-matchfinder is enabled, save edges only and exit
    if args.skip_matchfinder:
        print("\n  Skipping matchfinder processing (--skip-matchfinder enabled)")
        print("  Saving edges dictionary with cutoff parameters...")
        
        # Save edges dictionary with cutoff params for direct use with build_network_graph.py
        edges_output_path = args.output_edges_pkl
        if not edges_output_path:
            # Infer from output_pkl: replace 'paired_features' with 'edges' or add '_edges' suffix
            edges_output_path = str(args.output_pkl).replace('paired_features', 'edges')
            if edges_output_path == args.output_pkl:
                # If replacement didn't work, add suffix
                edges_output_path = args.output_pkl.replace('.pkl', '_edges.pkl')
        
        # Create wrapper dict with network_edges and cutoff_params
        edges_with_params = {
            'network_edges': network_edges,
            'cutoff_params': {
                'mz_cutoff': args.mz_cutoff,
                'rt_cutoff': args.rt_cutoff,
                'edges_cutoff': args.edges_cutoff
            }
        }
        
        with open(edges_output_path, 'wb') as f:
            pickle.dump(edges_with_params, f)
        print(f"✓ Saved edges dictionary (pickle): {os.path.basename(edges_output_path)}")
        
        # Also save edges dictionary as JSON
        edges_json_path = args.output_json.replace('paired_features', 'edges')
        if edges_json_path == args.output_json:
            # If replacement didn't work, modify the filename
            edges_json_path = args.output_json.replace('.json', '_edges.json')
        
        with open(edges_json_path, 'w') as f:
            json.dump(network_edges, f, indent=2)
        print(f"✓ Saved edges dictionary (JSON): {os.path.basename(edges_json_path)}")
        
        print(f"\n  Next step: Build network graph with:")
        print(f"    python build_network_graph.py --input_pkl {os.path.basename(edges_output_path)} --edge_cutoff <DISTANCE>\n")
        return
    
    # BUILD CONNECTED COMPONENTS FROM FILTERED EDGES
    # This is crucial: components are built AFTER edge filtering to ensure
    # that only features connected by edges passing the cutoff are grouped together
    print("\n  Building connected components from filtered edges...")
    paired_features = _build_connected_components_from_edges(network_edges, all_feature_data_lists)
    
    # Filter features within groups by best inter-feature match score cutoff
    def filter_features_in_group(match_group):
        """Remove features that don't have at least one score above cutoff to another feature"""
        individual_features = match_group.get('individual_features', [])
        
        # Keep only features with at least one score above cutoff
        filtered_individual = []
        for feature in individual_features:
            scores = feature.get('scores_to_features', {})
            if scores:
                max_score = max(scores.values())
                if max_score >= args.match_cutoff:
                    filtered_individual.append(feature)
            elif not scores and args.match_cutoff == 0.0:
                # Keep features with no scores only if cutoff is 0
                filtered_individual.append(feature)
        
        # Update pep_ident to only include those from remaining features
        remaining_pep_idents = set()
        for feature in filtered_individual:
            if 'pep_ident' in feature and feature['pep_ident']:
                remaining_pep_idents.update(feature['pep_ident'])
        
        # OPTIMIZED: Use pre-calculated distances from distances_to_group_features instead of O(n²) recalculation
        if len(filtered_individual) > 1:
            # Collect distances already calculated and stored on features
            distances = []
            for feature in filtered_individual:
                dist_dict = feature.get('distances_to_group_features', {})
                distances.extend(dist_dict.values())
            
            if distances:
                avg_distance = np.mean(distances)
                match_score = float(np.exp(-avg_distance / 100.0))
            else:
                match_score = 1.0
        else:
            match_score = 1.0
        
        # Update the group with filtered features and recalculated match_score
        match_group['individual_features'] = filtered_individual
        match_group['num_files_matched'] = len(filtered_individual)
        # Use unique files list and add comprehensive features list
        unique_files_filtered = sorted(list(set(f['filename'] for f in filtered_individual)))
        all_features_filtered = [f"{f['filename']}_feat{f['_feat_idx']}" for f in filtered_individual]
        match_group['files'] = unique_files_filtered
        match_group['all_features'] = all_features_filtered  # NEW: Comprehensive list of all features in group
        match_group['pep_ident'] = list(remaining_pep_idents) if remaining_pep_idents else None
        match_group['match_score'] = match_score
        return match_group
    
    # Filter features within each group
    filtered_features = [filter_features_in_group(f) for f in paired_features]
    # Remove groups with no features left
    filtered_features = [f for f in filtered_features if len(f.get('individual_features', [])) > 0]
    
    print(f"Filtered to {len(filtered_features)}/{len(paired_features)} groups with features scoring >= {args.match_cutoff}")
    
    # Count total features before filtering
    total_features_before_filter = sum(len(f.get('individual_features', [])) for f in paired_features)
    total_features_after_filter = sum(len(f.get('individual_features', [])) for f in filtered_features)
    
    print(f"Total features before filtering: {total_features_before_filter}")
    print(f"Total features after filtering: {total_features_after_filter}")
    
    # Define pep_ident analysis function (optional)
    def analyze_pep_ident_matching(groups):
        """Analyze pep_ident matching patterns in feature groups"""
        pep_ident_stats = {
            'total_groups': len(groups),
            'groups_with_common_pep_idents': 0,
            'groups_with_mismatching_pep_idents': 0,
            'groups_with_no_pep_idents': 0,
            'groups_with_partial_matching': 0,
            'pep_ident_match_distribution': {}
        }
        
        for group in groups:
            individual_features = group.get('individual_features', [])
            if not individual_features:
                continue
            
            # Get all unique pep_idents from all features
            all_feature_pep_idents = []
            for feature in individual_features:
                pep_idents = feature.get('pep_ident', [])
                all_feature_pep_idents.append(set(pep_idents) if pep_idents else set())
            
            # Count how many features have pep_idents
            features_with_pep = sum(1 for p in all_feature_pep_idents if p)
            
            if features_with_pep == 0:
                pep_ident_stats['groups_with_no_pep_idents'] += 1
                match_key = 'none'
            elif features_with_pep == len(individual_features):
                # All features have pep_idents - check if they share at least one common peptide
                common_pep_idents = set.intersection(*all_feature_pep_idents) if all_feature_pep_idents else set()
                if common_pep_idents:
                    # All features have at least one common pep_ident
                    pep_ident_stats['groups_with_common_pep_idents'] += 1
                    match_key = f'common_pep_{features_with_pep}_of_{len(individual_features)}'
                else:
                    # All features have pep_idents but no common ones (mismatching)
                    pep_ident_stats['groups_with_mismatching_pep_idents'] += 1
                    match_key = f'mismatching_pep_{features_with_pep}_of_{len(individual_features)}'
            else:
                pep_ident_stats['groups_with_partial_matching'] += 1
                match_key = f'partial_{features_with_pep}_of_{len(individual_features)}'
            
            # Track distribution
            pep_ident_stats['pep_ident_match_distribution'][match_key] = \
                pep_ident_stats['pep_ident_match_distribution'].get(match_key, 0) + 1
        
        return pep_ident_stats
    
    # Analyze pep_ident matching in groups, just for statistics
    def analyze_pep_ident_matching(groups):
        """Analyze pep_ident matching patterns in feature groups"""
        pep_ident_stats = {
            'total_groups': len(groups),
            'groups_with_common_pep_idents': 0,
            'groups_with_mismatching_pep_idents': 0,
            'groups_with_no_pep_idents': 0,
            'groups_with_partial_matching': 0,
            'pep_ident_match_distribution': {}
        }
        
        for group in groups:
            individual_features = group.get('individual_features', [])
            if not individual_features:
                continue
            
            # Get all unique pep_idents from all features
            all_feature_pep_idents = []
            for feature in individual_features:
                pep_idents = feature.get('pep_ident', [])
                all_feature_pep_idents.append(set(pep_idents) if pep_idents else set())
            
            # Count how many features have pep_idents
            features_with_pep = sum(1 for p in all_feature_pep_idents if p)
            
            if features_with_pep == 0:
                pep_ident_stats['groups_with_no_pep_idents'] += 1
                match_key = 'none'
            elif features_with_pep == len(individual_features):
                # All features have pep_idents - check if they share at least one common peptide
                common_pep_idents = set.intersection(*all_feature_pep_idents) if all_feature_pep_idents else set()
                if common_pep_idents:
                    # All features have at least one common pep_ident
                    pep_ident_stats['groups_with_common_pep_idents'] += 1
                    match_key = f'common_pep_{features_with_pep}_of_{len(individual_features)}'
                else:
                    # All features have pep_idents but no common ones (mismatching)
                    pep_ident_stats['groups_with_mismatching_pep_idents'] += 1
                    match_key = f'mismatching_pep_{features_with_pep}_of_{len(individual_features)}'
            else:
                pep_ident_stats['groups_with_partial_matching'] += 1
                match_key = f'partial_{features_with_pep}_of_{len(individual_features)}'
            
            # Track distribution
            pep_ident_stats['pep_ident_match_distribution'][match_key] = \
                pep_ident_stats['pep_ident_match_distribution'].get(match_key, 0) + 1
        
        return pep_ident_stats
    
    # OPTIMIZED: Analyze pep_ident matching only if requested (optional, can be slow)
    if args.analyze_pep_idents:
        print("Analyzing pep_ident matching patterns...")
        pep_ident_stats = analyze_pep_ident_matching(filtered_features)
    else:
        # Skip analysis by default (much faster)
        pep_ident_stats = {'total_groups': len(filtered_features), 'skipped': True}
    
    # Save results with statistics
    output_data = {
        'paired_features': filtered_features,
        'network_edges': network_edges,  # Include network graph structure
        'pep_ident_statistics': pep_ident_stats,
        'total_features_before_filter': total_features_before_filter,
        'total_features_after_filter': total_features_after_filter,
        'cutoff_params': {
            'mz_cutoff': args.mz_cutoff,
            'rt_cutoff': args.rt_cutoff,
            'edges_cutoff': args.edges_cutoff
        },
        'pairing_parameters': {
            'pairing_method': 'basic' if args.use_basic else 'optimized',
            'best_match_only': args.best_match_only,
            'distance_calc_before_scaling': args.distance_calc_before_scaling,
            'normalize_coordinates': normalize_coordinates,
            'postpair_normalize_coordinates': args.postpair_normalize_coordinates,
            'match_cutoff': args.match_cutoff,
            'input_files': len(feature_files),
            'skip_matchfinder': args.skip_matchfinder
        }
    }
    
    # OPTIMIZED: Save pickle (fast binary format)
    print("Saving output (pickle)...")
    with open(args.output_pkl, 'wb') as f:
        pickle.dump(output_data, f)
    print(f"✓ Saved paired features (pickle): {os.path.basename(args.output_pkl)}")
    
    # OPTIMIZED: Skip JSON by default (optional) - JSON serialization of 45K+ features is slow
    if not args.skip_json_output:
        print("Saving output (JSON)...  (this may take a while, use --skip_json_output to skip)")
        with open(args.output_json, 'w') as f:
            json.dump(output_data, f, indent=2)
        print(f"✓ Saved paired features (JSON): {os.path.basename(args.output_json)}")
    else:
        print(f"⊘ Skipped JSON output (--skip_json_output enabled)")
    
    # Save edges dictionary separately for direct use with build_network_graph.py
    edges_output_path = args.output_edges_pkl
    if not edges_output_path:
        # Infer from output_pkl: replace 'paired_features' with 'edges' or add '_edges' suffix
        edges_output_path = str(args.output_pkl).replace('paired_features', 'edges')
        if edges_output_path == args.output_pkl:
            # If replacement didn't work, add suffix
            edges_output_path = args.output_pkl.replace('.pkl', '_edges.pkl')
    
    with open(edges_output_path, 'wb') as f:
        pickle.dump(network_edges, f)
    print(f"✓ Saved edges dictionary (pickle): {os.path.basename(edges_output_path)}")
    print(f"  Optionally use with: python build_network_graph.py --input_pkl {os.path.basename(edges_output_path)}")
    
    # Also save edges JSON for reporting and visualization pipelines
    edges_json_path = edges_output_path.replace('.pkl', '.json')
    if not args.skip_json_output:
        try:
            # Convert network_edges dict to JSON-serializable format
            edges_for_json = {}
            for edge_key, edge_data in network_edges.items():
                # Convert tuple keys to strings for JSON serialization
                str_key = str(edge_key)
                edges_for_json[str_key] = edge_data
            
            with open(edges_json_path, 'w') as f:
                json.dump(edges_for_json, f, indent=2)
            print(f"✓ Saved edges dictionary (JSON): {os.path.basename(edges_json_path)}")
        except Exception as e:
            print(f"⊘ Failed to save edges JSON: {e}")
    else:
        print(f"⊘ Skipped edges JSON output (--skip_json_output enabled)")


if __name__ == "__main__":
    main()

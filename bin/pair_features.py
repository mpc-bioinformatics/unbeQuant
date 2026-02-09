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
from pathlib import Path
from typing import Dict, List
from scipy.spatial import cKDTree
from pdb import set_trace as bp

def feature_pairing_optimized(all_feature_data_lists: List[List[Dict]]) -> List[Dict]:
    """
    Optimized feature pairing using KD-trees and vectorized operations.
    Much faster for large datasets with multiple files.
    """
    if not all_feature_data_lists:
        return []
    
    if len(all_feature_data_lists) == 1:
        matches = []
        for feature in all_feature_data_lists[0]:
            match_group = feature.copy()
            match_group['files'] = [feature['filename']]
            match_group['match_score'] = 1.0
            matches.append(match_group)
        return matches
    
    print("  Pairing features (OPTIMIZED) using KD-tree nearest-neighbor...")
    
    # Prepare data structures for each file
    file_features_list = []
    kdtrees = []
    feature_coords = []
    
    for file_idx, features in enumerate(all_feature_data_lists):
        coords = np.array([[f['x_center'], f['y_center']] for f in features])
        feature_coords.append(coords)
        
        if len(coords) > 0:
            kdtrees.append(cKDTree(coords))
        else:
            kdtrees.append(None)
        
        features_with_idx = []
        for feat_idx, f in enumerate(features):
            f_copy = f.copy()
            f_copy['_feat_idx'] = feat_idx
            f_copy['_file_idx'] = file_idx
            features_with_idx.append(f_copy)
        file_features_list.append(features_with_idx)
    
    # Start matching from first file
    first_file_features = file_features_list[0]
    matches = []
    
    print(f"  Processing {len(first_file_features)} features from first file...")
    
    for feature_idx, ref_feature in enumerate(first_file_features):
        if feature_idx % 100 == 0:
            print(f"    Paired {feature_idx}/{len(first_file_features)}", end='\r')
        #TODO: multiply m/z-values to bring it to the same scala as the rt-values for better KD-Tree performance
        ref_coords = np.array([[ref_feature['x_center'], ref_feature['y_center']]])
        bp()
        match_group = {'matches': [ref_feature.copy()]}
        min_distances = []
        
        for file_idx in range(1, len(all_feature_data_lists)):
            if kdtrees[file_idx] is None or len(file_features_list[file_idx]) == 0:
                continue
            
            distance, idx = kdtrees[file_idx].query(ref_coords, k=1)
            
            if isinstance(distance, np.ndarray):
                distance = distance.flat[0]
            if isinstance(idx, np.ndarray):
                idx = idx.flat[0]
            
            if idx >= 0:
                nearest_feature = file_features_list[file_idx][int(idx)].copy()
                min_distances.append(distance)
                match_group['matches'].append(nearest_feature)
                
        
        # Calculate match score
        #TODO: Check if it makes sense!!!
        if min_distances:
            avg_distance = np.mean(min_distances)
            match_score = float(np.exp(-avg_distance / 100.0))
        else:
            match_score = 1.0
        
        # Merge pep_ident values
        all_pep_idents = set()
        for matched_feature in match_group['matches']:
            if 'pep_ident' in matched_feature and matched_feature['pep_ident']:
                all_pep_idents.update(matched_feature['pep_ident'])
        
        unique_pep_idents = list(all_pep_idents) if all_pep_idents else None
        
        # Vectorized averaging of positions
        match_coords = np.array([[f['x_center'], f['y_center']] for f in match_group['matches']])
        match_mz_rt = np.array([[f['mz_start'], f['mz_end'], f['rt_start'], f['rt_end']] 
                                 for f in match_group['matches']])
        
        avg_coords = np.mean(match_coords, axis=0)
        avg_mz_rt = np.mean(match_mz_rt, axis=0)
        
        # Calculate per-feature scores based on distance from consolidated center
        scored_features = []
        for i, f in enumerate(match_group['matches']):
            feature_dist = np.sqrt((f['x_center'] - avg_coords[0])**2 + (f['y_center'] - avg_coords[1])**2)
            feature_score = float(np.exp(-feature_dist / 50.0))
            f_scored = f.copy()
            f_scored['feature_score'] = feature_score
            
            # Calculate scores to other features in the group
            inter_feature_scores = {}
            for j, other_f in enumerate(match_group['matches']):
                if i != j:
                    mz_dist = (f['x_center'] - other_f['x_center']) ** 2
                    rt_dist = (f['y_center'] - other_f['y_center']) ** 2
                    dist = np.sqrt(mz_dist + rt_dist)
                    score = float(np.exp(-dist / 100.0))
                    inter_feature_scores[other_f['filename']] = score
            
            f_scored['scores_to_features'] = inter_feature_scores
            scored_features.append(f_scored)
        
        # Create consolidated match entry
        consolidated_match = {
            'match_score': match_score,
            'num_files_matched': len(match_group['matches']),
            'files': [f['filename'] for f in match_group['matches']],
            'pep_ident': unique_pep_idents,
            'x_center': float(avg_coords[0]),
            'y_center': float(avg_coords[1]),
            'mz_start': float(avg_mz_rt[0]),
            'mz_end': float(avg_mz_rt[1]),
            'rt_start': float(avg_mz_rt[2]),
            'rt_end': float(avg_mz_rt[3]),
            'individual_features': scored_features
        }
        
        matches.append(consolidated_match)
    #TODO: IGraph to build data structure
    #First KD-Tree on each file to find nearest neighbors - cutoff direkt auf euclidische/manhatten-distanz Distanz
    #manhatten-Distanz für einzelne berechnung der einzelnen Achsen
    #TODO: Maybe visualize in Gravis
    print(f"  Paired {len(matches)} feature groups")
    return matches


def feature_pairing(all_feature_data_lists: List[List[Dict]]) -> List[Dict]:
    """
    Basic feature pairing using nearest-neighbor principle.
    """
    if not all_feature_data_lists:
        return []
    
    if len(all_feature_data_lists) == 1:
        matches = []
        for feature in all_feature_data_lists[0]:
            match_group = feature.copy()
            match_group['files'] = [feature['filename']]
            match_group['match_score'] = 1.0
            matches.append(match_group)
        return matches
    
    print("  Pairing features using nearest-neighbor...")
    
    all_features = []
    for file_idx, features in enumerate(all_feature_data_lists):
        for feature in features:
            feature_with_file = feature.copy()
            feature_with_file['_file_idx'] = file_idx
            all_features.append(feature_with_file)
    
    first_file_features = [f for f in all_features if f['_file_idx'] == 0]
    matches = []
    
    for feature_idx, ref_feature in enumerate(first_file_features):
        if feature_idx % 100 == 0:
            print(f"    Pairing feature {feature_idx}/{len(first_file_features)}", end='\r')
        
        match_group = {'matches': [ref_feature.copy()]}
        min_distances = []
        
        for file_idx in range(1, len(all_feature_data_lists)):
            file_features = [f for f in all_features if f['_file_idx'] == file_idx]
            
            if not file_features:
                continue
            
            distances = []
            for f in file_features:
                mz_dist = (ref_feature['x_center'] - f['x_center']) ** 2
                rt_dist = (ref_feature['y_center'] - f['y_center']) ** 2
                dist = (mz_dist + rt_dist) ** 0.5
                distances.append((dist, f))
            
            if distances:
                min_dist, nearest_feature = min(distances, key=lambda x: x[0])
                min_distances.append(min_dist)
                match_group['matches'].append(nearest_feature.copy())
        
        if min_distances:
            avg_distance = sum(min_distances) / len(min_distances)
            match_score = np.exp(-avg_distance / 100.0)
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
        
        consolidated_match = {
            'match_score': float(match_score),
            'num_files_matched': len(match_group['matches']),
            'files': [f['filename'] for f in match_group['matches']],
            'pep_ident': unique_pep_idents,
            'x_center': float(np.mean([f['x_center'] for f in match_group['matches']])),
            'y_center': float(np.mean([f['y_center'] for f in match_group['matches']])),
            'mz_start': float(np.mean([f['mz_start'] for f in match_group['matches']])),
            'mz_end': float(np.mean([f['mz_end'] for f in match_group['matches']])),
            'rt_start': float(np.mean([f['rt_start'] for f in match_group['matches']])),
            'rt_end': float(np.mean([f['rt_end'] for f in match_group['matches']])),
            'individual_features': features_with_scores
        }
        matches.append(consolidated_match)
    
    print(f"  Paired {len(matches)} feature groups")
    return matches


def main():
    parser = argparse.ArgumentParser(
        description="Pair matching features across multiple files"
    )
    parser.add_argument("--input_dir", required=True, help="Directory containing feature pickle files")
    parser.add_argument("--output_pkl", required=True, help="Output pickle file path")
    parser.add_argument("--output_json", required=True, help="Output JSON file path")
    parser.add_argument("--optimize", action='store_true', help="Use optimized KD-tree based pairing")
    parser.add_argument("--match_cutoff", type=float, default=0.0, help="Minimum match score cutoff (0.0-1.0)")
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print(f"Feature Pairing")
    print(f"{'='*70}")
    
    # Find all feature pickle files
    feature_files = sorted(Path(args.input_dir).glob("*_feature_data.pkl"))
    
    if not feature_files:
        print("✗ No feature pickle files found!")
        return
    
    print(f"Found {len(feature_files)} feature file(s)")
    
    # Load all feature data
    all_feature_data_lists = []
    for pkl_file in feature_files:
        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)
            all_feature_data_lists.append(data)
    
    # Pair features
    if args.optimize:
        paired_features = feature_pairing_optimized(all_feature_data_lists)
    else:
        paired_features = feature_pairing(all_feature_data_lists)
    
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
        
        # Recalculate match_score based on remaining features only
        if len(filtered_individual) > 1:
            # Calculate average distance between remaining features
            distances = []
            for i, f1 in enumerate(filtered_individual):
                for j, f2 in enumerate(filtered_individual):
                    if i < j:
                        mz_dist = (f1['x_center'] - f2['x_center']) ** 2
                        rt_dist = (f1['y_center'] - f2['y_center']) ** 2
                        dist = np.sqrt(mz_dist + rt_dist)
                        distances.append(dist)
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
        match_group['files'] = [f['filename'] for f in filtered_individual]
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
    
    pep_ident_stats = analyze_pep_ident_matching(filtered_features)
    
    # Save results with statistics
    output_data = {
        'paired_features': filtered_features,
        'pep_ident_statistics': pep_ident_stats,
        'total_features_before_filter': total_features_before_filter,
        'total_features_after_filter': total_features_after_filter
    }
    
    with open(args.output_pkl, 'wb') as f:
        pickle.dump(output_data, f)
    print(f"✓ Saved paired features (pickle): {os.path.basename(args.output_pkl)}")
    
    with open(args.output_json, 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"✓ Saved paired features (JSON): {os.path.basename(args.output_json)}")


if __name__ == "__main__":
    main()

import pyopenms as oms
from pdb import set_trace as bp
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
import pickle
import time
import psutil
import os
import gc
import ast
from pathlib import Path
from typing import Dict, List, Tuple
import json
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.optimize import linear_sum_assignment

# ============================================================================
# Configuration
# ============================================================================
round_up_to = 2  # Number of decimal places to round m/z values
testing = False  # If True, only process first spectrum for testing
log_scale = False  # If True, apply logarithmic scaling to intensities
scale_colors = False  # If True, scale colors based on min/max intensities, only for visualization
invert_colors = True  # If True, invert color scale (high intensity = dark)
feature_mode = "CoM"  # Options: "rectangle", "CoM" (Center of Mass)
optimize_pairing = True  # If True, use optimized KD-tree based pairing (faster for large datasets)
feature_center_diagnostic = True  # If True, generate diagnostic plots for feature centers

# Directory paths - Update these as needed
mzML_directory = "/workspaces/unbeQuant/results/mgfs_mzmls/"
tsv_directory = "/workspaces/unbeQuant/results/quantifications/features_with_annotated_identifications/"
output_directory = "/workspaces/unbeQuant/results/feature_data_lists/"

# Feature visualization settings
box_color = "red"
box_width = 2
box_alpha = 0.7

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# ============================================================================
# Global RGB Dictionary (computed once at startup)
# ============================================================================
rgb_intensity_dict = None

def initialize_rgb_dict():
    """Precompute RGB dictionary once for all files."""
    global rgb_intensity_dict
    
    if scale_colors:
        print("Precomputing RGB intensity dictionary...")
        
        def intensity_to_rgb(i):
            r = 0
            g = 0
            b = 0
            
            # Distribute 24 bits round-robin: bits 0,3,6,9... → B, bits 1,4,7,10... → G, bits 2,5,8,11... → R
            for bit_pos in range(24):
                if i & (1 << bit_pos):
                    channel = bit_pos % 3
                    bit_index = bit_pos // 3
                    
                    if channel == 0:
                        b |= (1 << bit_index)
                    elif channel == 1:
                        g |= (1 << bit_index)
                    else:
                        r |= (1 << bit_index)
            
            return (r, g, b)
        
        rgb_intensity_dict = {i: intensity_to_rgb(i) for i in range(16777216)}
        print("RGB dictionary precomputed and ready")
    else:
        # Simple function for non-scaled colors
        rgb_intensity_dict = {}

def intensity_to_rgb_lookup(intensity_int):
    """Look up RGB values using precomputed dictionary."""
    if scale_colors and rgb_intensity_dict:
        return rgb_intensity_dict.get(intensity_int, (0, 0, 0))
    else:
        return rgb_int2tuple(intensity_int)

def rgb_int2tuple(rgbint):
    """Convert RGB integer to tuple for non-scaled colors."""
    return (rgbint // 256 // 256 % 256, rgbint // 256 % 256, rgbint % 256)

# ============================================================================
# Helper Functions
# ============================================================================

def parse_list(s):
    """Parse list strings from TSV - handles numeric values."""
    if isinstance(s, str):
        try:
            parsed = ast.literal_eval(s)
            return [float(x) for x in parsed]
        except:
            return []
    return []

def parse_list_flexible(s):
    """Parse list strings from TSV - handles both numeric and string values."""
    if isinstance(s, str):
        try:
            parsed = ast.literal_eval(s)
            # Try to convert to float, if it fails keep as string
            result = []
            for x in parsed:
                try:
                    result.append(float(x))
                except (ValueError, TypeError):
                    result.append(str(x))
            return result
        except:
            return []
    return []

def find_file_pairs(mzml_dir: str, tsv_dir: str) -> List[Tuple[str, str, str]]:
    """
    Find matching mzML and TSV file pairs.
    Returns list of tuples: (mzml_path, tsv_path, common_name)
    """
    pairs = []
    
    # Get all mzML files
    mzml_files = sorted(Path(mzml_dir).glob("*.mzML"))
    
    for mzml_path in mzml_files:
        filename_base = mzml_path.stem  # Name without extension
        
        # Look for corresponding TSV file
        tsv_path = Path(tsv_dir) / f"{filename_base}.tsv"
        
        if tsv_path.exists():
            pairs.append((str(mzml_path), str(tsv_path), filename_base))
            print(f"✓ Found pair: {filename_base}")
        else:
            print(f"✗ No TSV found for {filename_base}")
    
    return pairs

def process_mzml_file(mzml_path: str) -> Tuple[np.ndarray, List[float], List[float], Dict, Dict]:
    """
    Process a single mzML file.
    Returns: (spec_total, RTINSECONDS_arr, mz_total_arr, mz_dict, rt_dict)
    """
    print(f"\n  Loading mzML: {os.path.basename(mzml_path)}")
    
    exp = oms.MSExperiment()
    oms.MzMLFile().load(mzml_path, exp)
    
    spec_total = []
    RTINSECONDS_arr = []
    mz_arr = []
    counter = 0
    
    for spec in exp:
        if spec.getMSLevel() == 1:
            RT = spec.getRT()
            mz, intensity = spec.get_peaks()
            counter += 1
            RTINSECONDS_arr.append(round(RT, 8))
            
            if counter % 100 == 0:
                print(f"    Processed {counter} spectra", end='\r')
            
            if len(mz) == len(intensity) and len(mz) > 0:
                for i in range(len(mz)):
                    mz_i = mz[i]
                    spec_curr = [round(RT, 8), round(mz_i, round_up_to), intensity[i]]
                    spec_total.append(spec_curr)
                    mz_arr.append(mz_i)
                    bp()
            else:
                print(f"\n  Warning: Mismatch in spectrum at RT={RT}")
    
    # Convert to arrays and process
    spec_total = np.array(spec_total)
    mz_arr = np.array(mz_arr)
    mz_arr = np.round(mz_arr, round_up_to)
    
    # Create m/z range array
    mz_total_arr = np.arange(min(mz_arr), max(mz_arr) + (1/(10**round_up_to)), (1/(10**round_up_to)))
    mz_total_arr = np.round(mz_total_arr, round_up_to)
    mz_total_arr = list(mz_total_arr)
    
    # Create lookup dictionaries for O(1) access
    mz_dict = {v: i for i, v in enumerate(mz_total_arr)}
    rt_dict = {v: i for i, v in enumerate(RTINSECONDS_arr)}
    
    print(f"    Loaded {len(spec_total)} spectrum points, {len(RTINSECONDS_arr)} RT values")
    print(f"    m/z range: {min(mz_arr):.2f} - {max(mz_arr):.2f}")
    
    return spec_total, RTINSECONDS_arr, mz_total_arr, mz_dict, rt_dict

def create_heatmap_image(spec_total: np.ndarray, RTINSECONDS_arr: List[float], 
                         mz_total_arr: List[float], mz_dict: Dict, rt_dict: Dict,
                         output_path: str) -> np.ndarray:
    """
    Create and save heatmap image, return raw intensity array.
    Returns: img_raw (raw intensity values for center of mass calculation)
    """
    print(f"  Creating heatmap image...")
    
    img_height = int(len(RTINSECONDS_arr) + 10)
    mz_total_num = (max(mz_total_arr) - min(mz_total_arr)) * (10 ** round_up_to)
    img_width = int(mz_total_num) + 10
    
    # Check RAM requirements
    img_size_bytes = img_height * img_width * 3
    available_ram_bytes = psutil.virtual_memory().available
    use_memmap = img_size_bytes > (available_ram_bytes * 0.8)
    
    print(f"    Image size: {img_size_bytes / (1024**3):.2f} GB")
    print(f"    Available RAM: {available_ram_bytes / (1024**3):.2f} GB")
    print(f"    Using memmap: {use_memmap}")
    
    # Create arrays
    if use_memmap:
        temp_filepath = output_path.replace('.png', '_temp.npy')
        img_frame = np.memmap(temp_filepath, dtype=np.uint8, mode='w+', 
                              shape=(img_height, img_width, 3))
        img_frame[:] = 255
    else:
        img_frame = np.zeros((img_height, img_width, 3), dtype=np.uint8)
        img_frame[:] = 255
    
    img_raw = np.zeros((img_height, img_width), dtype=np.float32)
    
    # Normalize intensity values
    min_spec = spec_total[:, 2].min()
    max_spec = spec_total[:, 2].max()
    
    if log_scale:
        log_min = np.log1p(min_spec)
        log_max = np.log1p(max_spec)
        log_diff = log_max - log_min
    else:
        log_min = min_spec
        log_max = max_spec
        log_diff = max_spec - min_spec
    
    # Map spectrum points to image
    print(f"  Mapping {len(spec_total)} spectrum points to image...")
    for i in range(len(spec_total)):
        if i % 50000 == 0:
            print(f"    Mapped {i}/{len(spec_total)}", end='\r')
        
        try:
            index_mz = mz_dict[round(spec_total[i][1], round_up_to)]
            index_sec = rt_dict[round(spec_total[i][0], 8)]
            intensity = spec_total[i][2]
            bp()
            if log_scale:
                log_intensity = np.log1p(intensity)
            else:
                log_intensity = intensity
            
            if invert_colors:
                norm_log_intensity = 1 - ((log_intensity - log_min) / log_diff)
            else:
                norm_log_intensity = ((log_intensity - log_min) / log_diff)
            
            rgb_tuple = intensity_to_rgb_lookup(int(round(norm_log_intensity * 16777215)))
            img_frame[index_sec][index_mz] = rgb_tuple
            img_raw[index_sec][index_mz] = intensity  # Store full precision intensity
        except (KeyError, IndexError):
            pass
    
    # Flush memmap if used
    if use_memmap:
        img_frame.flush()
        # Load into memory for image saving
        img_frame_data = np.array(img_frame)
        del img_frame
        gc.collect()
        img_frame = img_frame_data
    
    # Save as PNG
    print(f"  Saving heatmap to {os.path.basename(output_path)}")
    img_pil = Image.fromarray(img_frame, mode='RGB')
    img_pil.save(output_path, format='PNG')
    img_pil = None
    
    gc.collect()
    
    return img_raw

def extract_feature_data(features_tsv: str, img_raw: np.ndarray, mz_dict: Dict, 
                        rt_dict: Dict, mz_total_arr: List[float], 
                        RTINSECONDS_arr: List[float], filename_base: str) -> List[Dict]:
    """
    Extract feature data and create overlays with boxes/CoM markers.
    """
    print(f"  Loading features from TSV...")
    features_data = pd.read_csv(features_tsv, sep='\t')
    print(f"    Loaded {len(features_data)} features")
    
    img_height = len(RTINSECONDS_arr) + 10
    mz_total_num = (max(mz_total_arr) - min(mz_total_arr)) * (10 ** round_up_to)
    img_width = int(mz_total_num) + 10
    
    feature_data_list = []
    box_count = 0
    skipped_count = 0
    
    # Color mapping for overlay
    color_map = {
        'red': (255, 0, 0, int(255 * box_alpha)),
        'blue': (0, 0, 255, int(255 * box_alpha)),
        'green': (0, 255, 0, int(255 * box_alpha)),
        'yellow': (255, 255, 0, int(255 * box_alpha)),
        'cyan': (0, 255, 255, int(255 * box_alpha)),
        'magenta': (255, 0, 255, int(255 * box_alpha)),
        'white': (255, 255, 255, int(255 * box_alpha)),
        'black': (0, 0, 0, int(255 * box_alpha))
    }
    box_color_rgba = color_map.get(box_color.lower(), (255, 0, 0, int(255 * box_alpha)))
    
    print(f"  Processing {len(features_data)} features...")
    
    for idx, row in features_data.iterrows():
        if idx % 100 == 0:
            print(f"    Processing feature {idx}/{len(features_data)}", end='\r')
        
        mz_starts = parse_list(row['l_mz_start'])
        mz_ends = parse_list(row['l_mz_end'])
        rt_starts = parse_list(row['l_rt_start'])
        rt_ends = parse_list(row['l_rt_end'])
        pep_ident = parse_list_flexible(row['l_pep_ident'])
        #bp()
        if mz_starts and mz_ends and rt_starts and rt_ends:
            mz_start = min(mz_starts + mz_ends)
            mz_end = max(mz_starts + mz_ends)
            rt_start = min(rt_starts + rt_ends)
            rt_end = max(rt_starts + rt_ends)
            
            mz_start_rounded = round(mz_start, round_up_to)
            mz_end_rounded = round(mz_end, round_up_to)
            
            try:
                x_start = mz_dict[mz_start_rounded]
                x_end = mz_dict[mz_end_rounded]
                y_start = rt_dict[round(rt_start, 8)]
                y_end = rt_dict[round(rt_end, 8)]
                
                # Calculate center of mass
                x_center_geo = (x_start + x_end) / 2
                y_center_geo = (y_start + y_end) / 2
                x_center = x_center_geo
                y_center = y_center_geo
                
                # Calculate intensity-weighted center of mass
                intensity_sum = 0
                for x in range(int(x_start), int(x_end) + 1):
                    for y in range(int(y_start), int(y_end) + 1):
                        if 0 <= x < img_raw.shape[1] and 0 <= y < img_raw.shape[0]:
                            intensity_sum += img_raw[y, x]
                
                if intensity_sum > 0:
                    for x in range(int(x_start), int(x_end) + 1):
                        for y in range(int(y_start), int(y_end) + 1):
                            if 0 <= x < img_raw.shape[1] and 0 <= y < img_raw.shape[0]:
                                x_vec = (x - x_center_geo) * (img_raw[y, x] / intensity_sum)
                                y_vec = (y - y_center_geo) * (img_raw[y, x] / intensity_sum)
                                x_center += x_vec
                                y_center += y_vec
                                
                if feature_center_diagnostic == True:             
                    plt.subplot(2, 2, 2)
                    plt.scatter([x_center], [y_center], color='red', s=100, marker='x', label='Center of Gravity')
                    plt.scatter([(x_start + x_end) / 2], [(y_start + y_end) / 2], color='blue', marker='o', label='Geometric Center')
                    plt.xlim(x_start - 2, x_end + 2)
                    plt.ylim(y_start - 2, y_end + 2)
                    plt.legend()
                    plt.title(f'Feature {idx} Centers')
                    plt.xlabel('X pixels')
                    plt.ylabel('Y pixels')
                    plt.gca().invert_yaxis()

                    plt.subplot(2, 1, 2)
                    heatmap_data = img_raw[y_start:y_end+1, x_start:x_end+1]
                    im = plt.imshow(heatmap_data, cmap='hot', aspect='auto', interpolation='bilinear')
                    plt.colorbar(im, label='Intensity')
                    plt.title(f'Feature {idx} Heatmap')
                    plt.xlabel('X pixels')
                    plt.ylabel('Y pixels')

                    plt.tight_layout()
                    plt.savefig(f'/workspaces/unbeQuant/results/feature_{idx}_analysis.png', dpi=100, bbox_inches='tight')
                    plt.close()
                    bp()
                
                feature_data_list.append({
                    'idx': idx,
                    'x_center_geo': float(x_center_geo),
                    'y_center_geo': float(y_center_geo),
                    'mz_start': float(mz_start),
                    'mz_end': float(mz_end),
                    'rt_start': float(rt_start),
                    'rt_end': float(rt_end),
                    'x_center': float(x_center),
                    'y_center': float(y_center),
                    'pep_ident': pep_ident,
                    'filename': filename_base
                })
                
                box_count += 1
                
            except KeyError:
                skipped_count += 1
    
    print(f"  Extracted {box_count} features ({skipped_count} skipped)")
    
    return feature_data_list

def feature_pairing(all_feature_data_lists: List[List[Dict]]) -> List[Dict]:
    """
    Pair matching features across multiple files using nearest-neighbor principle.
    
    Args:
        all_feature_data_lists: List of feature_data_lists from different files
    
    Returns:
        List of matched feature groups with scores
    """
    if not all_feature_data_lists:
        return []
    
    if len(all_feature_data_lists) == 1:
        # Single file - just add match score of 1.0
        matches = []
        for feature in all_feature_data_lists[0]:
            match_group = feature.copy()
            match_group['files'] = [feature['filename']]
            match_group['match_score'] = 1.0
            matches.append(match_group)
        return matches
    
    print("  Pairing features across files using nearest-neighbor...")
    
    # Flatten all features with file index
    all_features = []
    for file_idx, features in enumerate(all_feature_data_lists):
        for feature in features:
            feature_with_file = feature.copy()
            feature_with_file['_file_idx'] = file_idx
            all_features.append(feature_with_file)
    
    # Start with features from first file
    first_file_features = [f for f in all_features if f['_file_idx'] == 0]
    matches = []
    
    for feature_idx, ref_feature in enumerate(first_file_features):
        if feature_idx % 100 == 0:
            print(f"    Pairing feature {feature_idx}/{len(first_file_features)}", end='\r')
        
        # Find nearest neighbor in each other file
        match_group = {'matches': [ref_feature.copy()]}
        min_distances = []
        
        for file_idx in range(1, len(all_feature_data_lists)):
            file_features = [f for f in all_features if f['_file_idx'] == file_idx]
            
            if not file_features:
                continue
            
            # Calculate distances to all features in this file
            distances = []
            for f in file_features:
                # Euclidean distance in (m/z, RT) space
                mz_dist = (ref_feature['mz_start'] - f['mz_start']) ** 2
                rt_dist = (ref_feature['rt_start'] - f['rt_start']) ** 2
                dist = (mz_dist + rt_dist) ** 0.5
                distances.append((dist, f))
            
            # Find nearest
            if distances:
                min_dist, nearest_feature = min(distances, key=lambda x: x[0])
                min_distances.append(min_dist)
                match_group['matches'].append(nearest_feature.copy())
        
        # Calculate match score (inverse of average distance)
        if min_distances:
            avg_distance = sum(min_distances) / len(min_distances)
            # Normalize to 0-1 range using exponential decay
            match_score = np.exp(-avg_distance / 100.0)
        else:
            match_score = 1.0
        
        # Merge pep_ident values
        all_pep_idents = []
        for matched_feature in match_group['matches']:
            if 'pep_ident' in matched_feature and matched_feature['pep_ident']:
                all_pep_idents.extend(matched_feature['pep_ident'])
        
        # Remove duplicates and keep unique identifications
        unique_pep_idents = list(set(all_pep_idents)) if all_pep_idents else None
        
        # Create consolidated match entry
        consolidated_match = {
            'match_score': float(match_score),
            'num_files_matched': len(match_group['matches']),
            'files': [f['filename'] for f in match_group['matches']],
            'pep_ident': unique_pep_idents,
            # Take the average position across all matches
            'x_center': float(np.mean([f['x_center'] for f in match_group['matches']])),
            'y_center': float(np.mean([f['y_center'] for f in match_group['matches']])),
            # Take the average m/z and RT
            'mz_start': float(np.mean([f['mz_start'] for f in match_group['matches']])),
            'mz_end': float(np.mean([f['mz_end'] for f in match_group['matches']])),
            'rt_start': float(np.mean([f['rt_start'] for f in match_group['matches']])),
            'rt_end': float(np.mean([f['rt_end'] for f in match_group['matches']])),
            # Store individual feature data for reference
            'individual_features': match_group['matches']
        }
        
        matches.append(consolidated_match)
    
    print(f"  Paired {len(matches)} feature groups")
    
    return matches

def feature_pairing_optimized(all_feature_data_lists: List[List[Dict]]) -> List[Dict]:
    """
    Optimized feature pairing using KD-trees and vectorized operations.
    Much faster for large datasets with multiple files.
    
    Args:
        all_feature_data_lists: List of feature_data_lists from different files
    
    Returns:
        List of matched feature groups with scores
    """
    if not all_feature_data_lists:
        return []
    
    if len(all_feature_data_lists) == 1:
        # Single file - just add match score of 1.0
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
        # Extract coordinates for KD-tree (m/z, RT)
        coords = np.array([[f['mz_start'], f['rt_start']] for f in features])
        feature_coords.append(coords)
        
        # Build KD-tree for this file
        if len(coords) > 0:
            kdtrees.append(cKDTree(coords))
        else:
            kdtrees.append(None)
        
        # Add file index to each feature for reference
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
        
        ref_coords = np.array([[ref_feature['mz_start'], ref_feature['rt_start']]])
        match_group = {'matches': [ref_feature.copy()]}
        min_distances = []
        
        # Find nearest neighbor in each other file using KD-tree
        for file_idx in range(1, len(all_feature_data_lists)):
            if kdtrees[file_idx] is None or len(file_features_list[file_idx]) == 0:
                continue
            
            # Query KD-tree for nearest neighbor
            distance, idx = kdtrees[file_idx].query(ref_coords, k=1)
            
            # Handle both scalar and array returns from query
            if isinstance(distance, np.ndarray):
                distance = distance.flat[0]
            if isinstance(idx, np.ndarray):
                idx = idx.flat[0]
            
            if idx >= 0:
                nearest_feature = file_features_list[file_idx][int(idx)].copy()
                min_distances.append(distance)
                match_group['matches'].append(nearest_feature)
        
        # Calculate match score using vectorized operations
        if min_distances:
            avg_distance = np.mean(min_distances)
            match_score = float(np.exp(-avg_distance / 100.0))
        else:
            match_score = 1.0
        
        # Merge pep_ident values using set for O(1) lookup
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
            'individual_features': match_group['matches']
        }
        
        matches.append(consolidated_match)
    
    print(f"  Paired {len(matches)} feature groups")
    
    return matches

# ============================================================================
# Main Processing
# ============================================================================

def process_single_file_pair(mzml_path: str, tsv_path: str, filename_base: str) -> Tuple[bool, List[Dict]]:
    """
    Process a single mzML/TSV pair. 
    Returns: (success: bool, feature_data_list: List[Dict])
    """
    try:
        print(f"\n{'='*70}")
        print(f"Processing: {filename_base}")
        print(f"{'='*70}")
        
        # Step 1: Load and process mzML
        spec_total, RTINSECONDS_arr, mz_total_arr, mz_dict, rt_dict = process_mzml_file(mzml_path)
        
        # Step 2: Create heatmap
        heatmap_output = os.path.join(output_directory, f"{filename_base}_heatmap.png")
        img_raw = create_heatmap_image(spec_total, RTINSECONDS_arr, mz_total_arr, 
                                       mz_dict, rt_dict, heatmap_output)
        
        # Step 3: Extract feature data
        feature_data_list = extract_feature_data(tsv_path, img_raw, mz_dict, rt_dict,
                                                  mz_total_arr, RTINSECONDS_arr, filename_base)
        
        # Step 4: Export feature data list
        output_pickle = os.path.join(output_directory, f"{filename_base}_feature_data.pkl")
        output_json = os.path.join(output_directory, f"{filename_base}_feature_data.json")
        
        with open(output_pickle, 'wb') as f:
            pickle.dump(feature_data_list, f)
        print(f"  ✓ Saved pickle: {os.path.basename(output_pickle)}")
        
        # Also save as JSON for readability
        with open(output_json, 'w') as f:
            json.dump(feature_data_list, f, indent=2)
        print(f"  ✓ Saved JSON: {os.path.basename(output_json)}")
        
        # Clean up memory
        del spec_total, RTINSECONDS_arr, mz_total_arr, mz_dict, rt_dict, img_raw
        gc.collect()
        
        return True, feature_data_list
        
    except Exception as e:
        print(f"  ✗ Error processing {filename_base}: {e}")
        import traceback
        traceback.print_exc()
        return False, []

def main():
    """
    Main entry point.
    Uses the global optimize_pairing configuration variable.
    """
    print("\n" + "="*70)
    print("Multi-File mzML to Feature Data Converter")
    print("="*70)
    if optimize_pairing:
        print("(Using OPTIMIZED pairing algorithm)")
    
    # Initialize RGB dictionary once
    initialize_rgb_dict()
    
    # Find file pairs
    print(f"\nSearching for file pairs...")
    print(f"  mzML directory: {mzML_directory}")
    print(f"  TSV directory: {tsv_directory}")
    
    file_pairs = find_file_pairs(mzML_directory, tsv_directory)
    
    if not file_pairs:
        print("\n✗ No matching mzML/TSV pairs found!")
        return
    
    print(f"\n✓ Found {len(file_pairs)} file pair(s)")
    
    # Process each pair and collect feature data lists
    start_time = time.time()
    successful = 0
    failed = 0
    all_feature_data_lists = []
    
    for mzml_path, tsv_path, filename_base in file_pairs:
        success, feature_data_list = process_single_file_pair(mzml_path, tsv_path, filename_base)
        if success:
            successful += 1
            all_feature_data_lists.append(feature_data_list)
        else:
            failed += 1
    
    # Cross-file feature pairing
    if all_feature_data_lists:
        print(f"\n{'='*70}")
        print("CROSS-FILE FEATURE PAIRING")
        print(f"{'='*70}")
        
        # Choose pairing algorithm based on global config
        if optimize_pairing:
            paired_features = feature_pairing_optimized(all_feature_data_lists)
        else:
            paired_features = feature_pairing(all_feature_data_lists)
        
        # Export paired features
        paired_output_pickle = os.path.join(output_directory, "paired_features.pkl")
        paired_output_json = os.path.join(output_directory, "paired_features.json")
        
        with open(paired_output_pickle, 'wb') as f:
            pickle.dump(paired_features, f)
        print(f"\n✓ Saved paired features (pickle): {os.path.basename(paired_output_pickle)}")
        
        with open(paired_output_json, 'w') as f:
            json.dump(paired_features, f, indent=2)
        print(f"✓ Saved paired features (JSON): {os.path.basename(paired_output_json)}")
    
    # Summary
    elapsed_time = time.time() - start_time
    print("\n" + "="*70)
    print("PROCESSING SUMMARY")
    print("="*70)
    print(f"Files processed: {successful}/{successful + failed}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    if all_feature_data_lists:
        print(f"Paired feature groups: {len(paired_features)}")
    print(f"Total time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
    print(f"Output directory: {output_directory}")
    print("="*70 + "\n")

if __name__ == "__main__":
    main()

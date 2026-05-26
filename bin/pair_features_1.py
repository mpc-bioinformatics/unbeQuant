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
import gc
import tracemalloc
import psutil
from pathlib import Path
from typing import Dict, List, Tuple
from scipy.spatial import cKDTree
from pdb import set_trace as bp

# Import model classes from compare_retention_time_alignment for unpickling
sys.path.insert(0, os.path.dirname(__file__))
try:
    from compare_retention_time_alignment import LOESSInterpolator, PiecewisePolynomial
except ImportError:
    # Define fallback classes if import fails
    class LOESSInterpolator:
        """LOESS interpolator that can be pickled."""
        def __init__(self, x_data, y_data):
            self.x_data = np.array(x_data)
            self.y_data = np.array(y_data)
        
        def __call__(self, x_val):
            if isinstance(x_val, (list, np.ndarray)):
                return np.array([self._eval_scalar(x) for x in np.atleast_1d(x_val)])
            else:
                return self._eval_scalar(x_val)
        
        def _eval_scalar(self, x_val):
            if x_val <= self.x_data[0]:
                return self.y_data[0]
            if x_val >= self.x_data[-1]:
                return self.y_data[-1]
            idx = np.searchsorted(self.x_data, x_val)
            x0, x1 = self.x_data[idx-1], self.x_data[idx]
            y0, y1 = self.y_data[idx-1], self.y_data[idx]
            return y0 + (y1 - y0) * (x_val - x0) / (x1 - x0)
    
    class PiecewisePolynomial:
        """Piecewise polynomial model that can be pickled (fallback)."""
        pass


def get_memory_info(label="", print_details=False):
    """
    Get comprehensive memory usage info using multiple methods.
    Returns RSS in MB, but prints all available metrics.
    """
    # Start tracemalloc if not already started
    if not tracemalloc.is_tracing():
        tracemalloc.start()
    
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    
    # Get all memory metrics
    rss_mb = mem_info.rss / (1024 * 1024)  # Resident Set Size - actual physical memory
    vms_mb = mem_info.vms / (1024 * 1024)  # Virtual Memory Size - all mapped memory
    
    # Get tracemalloc peak (Python object allocations only)
    current, peak = tracemalloc.get_traced_memory()
    current_mb = current / (1024 * 1024)
    peak_mb = peak / (1024 * 1024)
    
    if print_details:
        print(f"[COMPREHENSIVE MEMORY TRACKING] {label}")
        print(f"  RSS (physical memory in use):       {rss_mb:8.1f} MB")
        print(f"  VSZ (virtual memory mapped):        {vms_mb:8.1f} MB")
        print(f"  Tracemalloc current (Python objs):  {current_mb:8.1f} MB")
        print(f"  Tracemalloc peak (Python objs):     {peak_mb:8.1f} MB")
        try:
            percent = process.memory_percent()
            print(f"  System memory percent:              {percent:8.1f} %")
        except:
            pass
        # Show difference from baseline
        print()
    
    return rss_mb


def get_object_size_recursive(obj, seen=None, depth=0, max_depth=3):
    """
    Recursively estimate size of a Python object including nested structures.
    Useful for estimating dictionary/list memory usage.
    """
    if seen is None:
        seen = set()
    
    obj_id = id(obj)
    if obj_id in seen or depth > max_depth:
        return 0
    
    seen.add(obj_id)
    size = sys.getsizeof(obj)
    
    if isinstance(obj, dict):
        for k, v in obj.items():
            size += get_object_size_recursive(k, seen, depth + 1, max_depth)
            size += get_object_size_recursive(v, seen, depth + 1, max_depth)
    elif isinstance(obj, (list, tuple)):
        for item in obj:
            size += get_object_size_recursive(item, seen, depth + 1, max_depth)
    elif isinstance(obj, np.ndarray):
        size += obj.nbytes  # Direct array data
    
    return size


def print_structure_sizes(label, structures_dict):
    """Print estimated memory sizes of key data structures."""
    print(f"\n[STRUCTURE SIZES] {label}")
    total_size = 0
    for name, obj in structures_dict.items():
        if obj is None:
            print(f"  {name:30s}: Not allocated")
            continue
        try:
            if isinstance(obj, dict):
                # For dictionaries, estimate based on size of content
                if isinstance(obj, dict) and len(obj) > 0:
                    size = get_object_size_recursive(obj)
                    size_mb = size / (1024 * 1024)
                    print(f"  {name:30s}: {size_mb:8.1f} MB ({len(obj)} entries)")
                    total_size += size
                else:
                    print(f"  {name:30s}: <1 MB (empty)")
            elif isinstance(obj, list):
                size = get_object_size_recursive(obj)
                size_mb = size / (1024 * 1024)
                print(f"  {name:30s}: {size_mb:8.1f} MB ({len(obj)} entries)")
                total_size += size
            elif isinstance(obj, np.ndarray):
                size_mb = obj.nbytes / (1024 * 1024)
                print(f"  {name:30s}: {size_mb:8.1f} MB (shape {obj.shape})")
                total_size += obj.nbytes
            else:
                size = sys.getsizeof(obj)
                size_mb = size / (1024 * 1024)
                print(f"  {name:30s}: {size_mb:8.1f} MB")
                total_size += size
        except Exception as e:
            print(f"  {name:30s}: Error measuring - {e}")
    
    print(f"  {'TOTAL ESTIMATED':30s}: {total_size / (1024 * 1024):8.1f} MB")
    print()


def build_trafoxmls_fid_dict(trafoxmls_dir: str) -> Dict[str, float]:
    """
    Build a dictionary mapping feature IDs (f_<node>) to transformed RT values (to=)
    from trafoXML files.
    
    The trafoXML files contain <Pair note="<node>" to="<transformed_rt>"/> elements.
    The node attribute is converted to openms_fid format: f_<node>
    
    Args:
        trafoxmls_dir: Directory containing trafoXML files
    
    Returns:
        Dict mapping f_<node> (string) -> transformed_rt (float)
    """
    import xml.etree.ElementTree as ET
    from pathlib import Path
    import os
    
    fid_to_rt_dict = {}
    
    # Convert to Path object and resolve to absolute path
    dir_path = Path(trafoxmls_dir).resolve()
    
    print(f"[DEBUG] TrafoXML directory resolution:")
    print(f"  Input path: {repr(trafoxmls_dir)}")
    print(f"  Resolved path: {repr(str(dir_path))}")
    print(f"  Path exists: {dir_path.exists()}")
    print(f"  Path is dir: {dir_path.is_dir()}")
    
    if not dir_path.exists():
        raise FileNotFoundError(f"TrafoXML directory does not exist: {trafoxmls_dir}")
    
    if not dir_path.is_dir():
        raise NotADirectoryError(f"TrafoXML path is not a directory: {trafoxmls_dir}")
    
    # List all files ending in .trafoXML
    trafo_files = sorted(list(dir_path.glob("**/*.trafoXML")))
    
    print(f"[DEBUG] Glob search results:")
    print(f"  Pattern: {dir_path}/**/*.trafoXML")
    print(f"  Found {len(trafo_files)} trafoXML files")
    
    if not trafo_files:
        # Try non-recursive glob as well
        trafo_files_nonrec = sorted(list(dir_path.glob("*.trafoXML")))
        print(f"  Non-recursive glob (*.trafoXML): {len(trafo_files_nonrec)} files")
        if trafo_files_nonrec:
            trafo_files = trafo_files_nonrec
        else:
            # List what files are actually in the directory
            all_files = list(dir_path.glob("*"))
            print(f"  All files in directory ({len(all_files)} total):")
            for f in sorted(all_files)[:10]:
                print(f"    - {f.name}")
            if len(all_files) > 10:
                print(f"    ... and {len(all_files) - 10} more")
            raise FileNotFoundError(f"No .trafoXML files found in: {trafoxmls_dir}")
    
    print(f"[PROGRESS] Building feature ID -> transformed RT mappings from {len(trafo_files)} trafoXML files")
    
    for trafo_file in trafo_files:
        try:
            tree = ET.parse(str(trafo_file))
            root = tree.getroot()
            
            transformation = root.find('Transformation')
            if transformation is None:
                continue
                
            pairs_elem = transformation.find('Pairs')
            if pairs_elem is None:
                continue
            
            pairs = pairs_elem.findall('Pair')
            for pair in pairs:
                node_str = pair.get('note')  # Node ID from trafoXML
                to_rt_str = pair.get('to')   # Transformed RT
                
                if node_str and to_rt_str:
                    # Convert node ID to openms_fid format by prepending "f_"
                    fid = f"f_{node_str}"
                    to_rt = float(to_rt_str)
                    fid_to_rt_dict[fid] = to_rt
            
            print(f"  ✓ {trafo_file.name}: {len(pairs)} feature transformations loaded")
            
        except Exception as e:
            print(f"  ✗ Warning: Failed to parse {trafo_file.name}: {e}")
    
    print(f"  ✓ Total feature ID mappings: {len(fid_to_rt_dict)}")
    return fid_to_rt_dict


def apply_aligntree_rt_transformations(all_feature_data_lists: List[List[Dict]], 
                                       fid_to_rt_dict: Dict[str, float]) -> List[List[Dict]]:
    """
    Apply AlignTree RT transformations by matching feature IDs (openms_fid).
    Updates y_center to the transformed RT value from trafoXML files.
    
    Args:
        all_feature_data_lists: List of feature lists
        fid_to_rt_dict: Dict mapping openms_fid (f_<node>) -> transformed_rt
    
    Returns:
        Transformed feature lists with updated y_center values
    """
    transformed_lists = []
    total_matched = 0
    total_features = 0
    
    # DEBUG: Show sample from fid_to_rt_dict
    print(f"[DEBUG] fid_to_rt_dict sample (first 3 items, format: f_<node> -> transformed_rt):")
    for (fid, rt) in list(fid_to_rt_dict.items())[:3]:
        print(f"  {repr(fid)}: {rt:.2f}")
    
    for file_idx, feature_list in enumerate(all_feature_data_lists):
        transformed_features = []
        matched = 0
        
        for feature in feature_list:
            total_features += 1
            feature_copy = feature.copy()
            
            # Get feature ID (openms_fid field)
            fid = feature.get('openms_fid')
            
            if fid and fid in fid_to_rt_dict:
                # Apply transformation
                new_rt = fid_to_rt_dict[fid]
                feature_copy['y_center'] = new_rt
                matched += 1
                total_matched += 1
            
            transformed_features.append(feature_copy)
        
        transformed_lists.append(transformed_features)
        print(f"  File {file_idx + 1}: {matched}/{len(feature_list)} features transformed (y_center updated)")
    
    print(f"\n[PROGRESS] Total features transformed: {total_matched}/{total_features}")
    return transformed_lists
# ============================================================================
# RT Drift Correction Functions
# ============================================================================

def evaluate_polynomial(x_value: float, coefficients: List[float]) -> float:
    """
    Evaluate polynomial using coefficients from numpy.poly1d.
    
    Coefficients are stored in descending order of powers:
    [a_n, a_{n-1}, ..., a_1, a_0] represents:
    a_n*x^n + a_{n-1}*x^{n-1} + ... + a_1*x + a_0
    
    Args:
        x_value: Point at which to evaluate
        coefficients: List of polynomial coefficients (highest degree first)
    
    Returns:
        Evaluated polynomial value
    """
    if not coefficients:
        return 0.0
    
    # Use numpy.poly1d for evaluation
    poly = np.poly1d(coefficients)
    return float(poly(x_value))


def load_rt_corrections_json(json_path: str) -> Dict[Tuple[str, str], Dict]:
    """
    Load RT corrections from JSON file (supports all fitting methods).
    
    Args:
        json_path: Path to JSON file with fitted corrections
    
    Returns:
        Dictionary with (file1_name, file2_name) tuples as keys and correction dicts as values
    
    Raises:
        FileNotFoundError: If JSON file doesn't exist
        ValueError: If JSON format is invalid
    """
    import json
    import pickle
    import base64
    
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"RT correction JSON file not found: {json_path}")
    
    # Check if file is empty or only contains whitespace
    file_size = os.path.getsize(json_path)
    if file_size == 0:
        raise ValueError(f"RT correction JSON file is empty: {json_path}")
    
    with open(json_path, 'r') as f:
        content = f.read().strip()
        if not content:
            raise ValueError(f"RT correction JSON file contains no data: {json_path}")
        
        data = json.loads(content)
    
    if not isinstance(data, dict) or 'corrections' not in data:
        raise ValueError(f"JSON must have 'corrections' key")
    
    corrections = data['corrections']
    if not corrections:
        raise ValueError(f"RT correction JSON has no corrections (empty 'corrections' dict)")
    
    result = {}
    
    for pair_key, correction_data in corrections.items():
        # pair_key format: "file1_vs_file2"
        parts = pair_key.split('_vs_')
        if len(parts) != 2:
            raise ValueError(f"Invalid pair key format: {pair_key}, expected 'file1_vs_file2'")
        
        file1, file2 = parts
        # Remove .mzML extension if present
        file1 = file1.replace('.mzML', '')
        file2 = file2.replace('.mzML', '')
        
        result[(file1, file2)] = correction_data
    
    return result


def apply_rt_correction_json(rt_value: float, correction_dict: Dict) -> float:
    """
    Calculate RT correction amount using loaded JSON correction data.
    
    Supports multiple fitting methods (polynomial, spline, rbf, piecewise, loess).
    All methods use RT values as input (x-axis is RT from file 1).
    
    Args:
        rt_value: Retention time value (from reference file) used as input to the correction model
        correction_dict: Correction data dict with 'method' and 'serialization' keys
    
    Returns:
        RT correction amount (to be subtracted from original RT)
    
    Raises:
        ValueError: If method is unsupported or deserialization fails
    """
    import json
    import pickle
    import base64
    import numpy as np
    
    method = correction_dict.get('method', 'polynomial')
    serialization = correction_dict.get('serialization', {})
    
    try:
        if method == 'polynomial':
            # Extract coefficients and evaluate
            coefficients = serialization.get('coefficients', [])
            if not coefficients:
                return 0.0  # No correction available
            
            rt_diff = evaluate_polynomial(rt_value, coefficients)
            return rt_diff
        
        elif method in ['rbf', 'piecewise', 'loess', 'spline']:
            # Deserialize pickle-based models
            pickle_str = serialization.get('model_pickle')
            if not pickle_str:
                return 0.0
            
            model_bytes = base64.b64decode(pickle_str)
            model = pickle.loads(model_bytes)
            
            # Model is callable: rt_diff = model(rt_value)
            rt_diff = model(rt_value)
            return rt_diff
        
        else:
            raise ValueError(f"Unsupported fitting method: {method}")
    
    except Exception as e:
        raise ValueError(f"Failed to apply correction (method={method}): {e}")


def apply_rt_corrections_json_to_features(features: List[Dict], correction_dict: Dict, rt_source: str = 'y_center', mz_source: str = 'x_center') -> List[Dict]:
    """
    Apply RT correction to a list of features using JSON model.
    
    The correction works on RT values but is predicted from m/z:
    corrected_rt = rt - model(mz)
    
    Args:
        features: List of feature dictionaries
        correction_dict: Correction data dict with serialization info
        rt_source: Which RT field to correct
        mz_source: Which m/z field to use for prediction
    
    Returns:
        List of features with corrected RT coordinates
    """
    corrected_features = []
    
    for feature in features:
        feature_copy = feature.copy()
        
        # Get m/z (used as input to model)
        if mz_source not in feature:
            continue  # Skip if m/z not available
        
        mz_value = feature[mz_source]
        
        # Get original RT
        if rt_source not in feature:
            continue  # Skip if RT not available
        
        rt_value = feature[rt_source]
        
        # Apply correction
        try:
            # Equations are RT vs RT-Diff, so input is RT value, not m/z
            rt_diff = apply_rt_correction_json(rt_value, correction_dict)
            corrected_rt = rt_value - rt_diff
            
            # Backup original and update
            if 'y_center_original' not in feature_copy:
                feature_copy['y_center_original'] = feature_copy['y_center']
            
            feature_copy['y_center'] = corrected_rt
            corrected_features.append(feature_copy)
        except Exception as e:
            # If correction fails, keep original
            corrected_features.append(feature_copy)
    
    return corrected_features


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
                # Handle both index-based format (during pairing) and restored format (after JSON restoration)
                ref_file_idx = edge.get('ref_file_idx')
                match_file_idx = edge.get('match_file_idx')
                current_feat_idx = edge.get('current_feature_idx', edge.get('current_file', {}).get('feature_idx'))
                matched_feat_idx = edge.get('matched_feature_idx', edge.get('matched_file', {}).get('feature_idx'))
                
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
    Rescale edge distances using balanced coordinate axes.
    
    This rescales the smaller axis to match the larger axis's range, ensuring
    both dimensions contribute equally to distance calculations. Then calculates
    euclidean distance in this rescaled space (no further normalization to [0,1]).

    Preserves original coordinate distances while adding rescaled distance fields.
    
    Original distance is kept in: distance
    Rescaled distance is added to: distance_rescaled
    
    This updates only the edge distance fields so pairing decisions remain
    unchanged when coordinate cutoffs are used.

    Args:
        edges: Dictionary of edges {ref_file_idx: [edge_dict, ...]}
        all_feature_data_lists: List of feature data lists for coordinate lookup

    Returns:
        Updated edges dictionary with original and rescaled distance values
    """
    if not all_feature_data_lists:
        print("  ✗ WARNING: Feature data not provided for post-pair rescaling, skipping")
        return edges

    all_x_values = []
    all_y_values = []
    for features in all_feature_data_lists:
        for f in features:
            all_x_values.append(f['x_center'])
            all_y_values.append(f['y_center'])

    if not all_x_values or not all_y_values:
        print("  ✗ WARNING: Empty feature coordinates for post-pair rescaling, skipping")
        return edges

    min_x = float(np.min(all_x_values))
    max_x = float(np.max(all_x_values))
    min_y = float(np.min(all_y_values))
    max_y = float(np.max(all_y_values))

    range_x = max_x - min_x if max_x > min_x else 1.0
    range_y = max_y - min_y if max_y > min_y else 1.0
    
    # Find the larger range and scale both axes to it
    max_range = max(range_x, range_y)
    scale_x = max_range / range_x if range_x > 0 else 1.0
    scale_y = max_range / range_y if range_y > 0 else 1.0

    print("  Applying post-pair coordinate rescaling to edge distances...")
    print(f"    X range: {min_x:.2f} to {max_x:.2f} (range: {range_x:.2f}, scale: {scale_x:.4f})")
    print(f"    Y range: {min_y:.2f} to {max_y:.2f} (range: {range_y:.2f}, scale: {scale_y:.4f})")
    print(f"    Both axes rescaled to range: {max_range:.2f}")

    updated_edges = {}
    updated_count = 0
    skipped_count = 0

    for file_idx, edge_list in edges.items():
        updated_list = []
        for edge in edge_list:
            ref_file_idx = edge.get('ref_file_idx')
            match_file_idx = edge.get('match_file_idx')
            # Handle both index-based format (during pairing) and restored format (after JSON restoration)
            current_feat_idx = edge.get('current_feature_idx', edge.get('current_file', {}).get('feature_idx'))
            matched_feat_idx = edge.get('matched_feature_idx', edge.get('matched_file', {}).get('feature_idx'))

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

            # Keep original coordinate differences
            mz_distance = abs(ref_feature['x_center'] - matched_feature['x_center'])
            rt_distance = abs(ref_feature['y_center'] - matched_feature['y_center'])
            distance_original = float(np.sqrt(mz_distance ** 2 + rt_distance ** 2))
            
            # Calculate rescaled coordinates (offset + scale, but keeping range proportional)
            ref_x_rescaled = (ref_feature['x_center'] - min_x) * scale_x
            ref_y_rescaled = (ref_feature['y_center'] - min_y) * scale_y
            match_x_rescaled = (matched_feature['x_center'] - min_x) * scale_x
            match_y_rescaled = (matched_feature['y_center'] - min_y) * scale_y
            
            # Calculate rescaled coordinate differences
            mz_distance_rescaled = abs(ref_x_rescaled - match_x_rescaled)
            rt_distance_rescaled = abs(ref_y_rescaled - match_y_rescaled)
            
            # Calculate euclidean distance in rescaled space
            distance_rescaled = float(np.sqrt(mz_distance_rescaled ** 2 + rt_distance_rescaled ** 2))

            edge_copy = edge.copy()
            if 'distance_raw' not in edge_copy:
                edge_copy['distance_raw'] = edge_copy.get('distance')

            # Update distance to original euclidean distance
            edge_copy['distance'] = distance_original
            
            # Add rescaled versions (when axes need balancing)
            edge_copy['distance_rescaled'] = distance_rescaled  # Euclidean in rescaled space

            updated_list.append(edge_copy)
            updated_count += 1

        updated_edges[file_idx] = updated_list

    print(f"  Post-pair rescaling updated {updated_count} edges, skipped {skipped_count}")
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
    
    print("[DEBUG_CC] Entering _build_connected_components_from_edges...")
    mem_cc_start = get_memory_info("CC START", print_details=True)
    
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
    
    print("[DEBUG_CC] Starting edge processing loop...")
    # Single pass through edges: union vertices and store distances + neighbors
    for ref_file_idx, edge_list in edges.items():
        print(f"[DEBUG_CC]   Processing file_idx {ref_file_idx}, {len(edge_list)} edges")
        for edge_num, edge in enumerate(edge_list):
            if edge_num % 100000 == 0 and edge_num > 0:
                print(f"[DEBUG_CC]     Processed {edge_num}/{len(edge_list)} edges...")
            # Handle both index-based format (during pairing) and restored format (after JSON restoration)
            # Index-based: 'current_feature_idx', restored: 'current_file'['feature_idx']
            v1_file = edge.get('current_file_idx', edge.get('ref_file_idx'))
            v1_feat = edge.get('current_feature_idx', edge.get('current_file', {}).get('feature_idx'))
            v2_file = edge.get('matched_file_idx', edge.get('match_file_idx'))
            v2_feat = edge.get('matched_feature_idx', edge.get('matched_file', {}).get('feature_idx'))
            
            v1 = (v1_file, v1_feat)
            v2 = (v2_file, v2_feat)
            distance = edge.get('distance', 0)
            
            all_vertices.add(v1)
            all_vertices.add(v2)
            
            # Union the vertices (group them)
            union(v1, v2)
            
            # Store distance for both directions (edges are undirected) AND neighbor mappings
            edge_distances[frozenset([v1, v2])] = distance
            neighbors[v1][v2] = distance
            neighbors[v2][v1] = distance
    
    print(f"[DEBUG_CC] Edge processing complete. all_vertices: {len(all_vertices)}, neighbors entries: {len(neighbors)}")
    mem_after_edge_proc = get_memory_info("CC AFTER EDGE PROCESSING", print_details=True)
    print(f"[DEBUG_CC] Memory delta after edge processing: {mem_after_edge_proc - mem_cc_start:.1f} MB")
    
    # Group vertices by their root (component)
    print("[DEBUG_CC] Grouping vertices by component...")
    components = defaultdict(set)
    for vertex in all_vertices:
        root = find(vertex)
        components[root].add(vertex)

    components = list(components.values())
    print(f"[DEBUG_CC] Component grouping complete. Total components: {len(components)}")
    mem_after_grouping = get_memory_info("CC AFTER GROUPING", print_details=True)
    print(f"[DEBUG_CC] Memory delta after grouping: {mem_after_grouping - mem_after_edge_proc:.1f} MB")
    
    # Build consolidated match entries for each component
    matches = []
    total_components = len(components)
    print("[DEBUG_CC] Starting component consolidation...")
    mem_before_consolidation = get_memory_info("CC BEFORE CONSOLIDATION", print_details=True)
    
    for comp_idx, component_vertices in enumerate(components):
        if comp_idx % max(1, total_components // 10) == 0:
            if comp_idx > 0:
                mem_now = get_memory_info(f"CC CONSOLIDATION {comp_idx}/{total_components}", print_details=False)
                print(f"    Processing component {comp_idx+1}/{total_components} (mem: {mem_now:.1f} MB)", end='\r', flush=True)
            else:
                print(f"    Processing component {comp_idx+1}/{total_components}", end='\r', flush=True)
        
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
                for v_neighbor, neighbor_distance in neighbors[v_current].items():
                    # Confirm neighbor is in this component (should always be true)
                    if v_neighbor in vertex_to_feature_idx:
                        other_feat_idx = vertex_to_feature_idx[v_neighbor]
                        neighbor_feature = component_features[other_feat_idx]
                        other_key = f"{neighbor_feature['filename']}_feat{v_neighbor[1]}"
                        distances_to_others[other_key] = float(neighbor_distance)
            
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
    mem_cc_end = get_memory_info("CC END", print_details=True)
    print(f"[DEBUG_CC] Total memory delta in CC function: {mem_cc_end - mem_cc_start:.1f} MB")
    return matches


def extract_and_save_unpaired_vertices(network_edges: Dict, all_feature_data_lists: List[List[Dict]], output_path: str = None) -> List[Dict]:
    """
    Extract vertices that have no edges in the network (unpaired), and save to JSON.
    
    Compares all vertices in all_feature_data_lists with vertices appearing in network_edges.
    Any vertex NOT present in network_edges is marked as unpaired.
    
    Args:
        network_edges: Dictionary {ref_file_idx: [edge_dict, ...]}
        all_feature_data_lists: List of feature lists from each file
        output_path: Optional path to save unpaired vertices JSON (if None, only returns list)
    
    Returns:
        List of unpaired vertex dictionaries
    """
    
    # Build set of all vertices that appear in edges (fast O(n) operation)
    # This is the ONLY way to determine if a vertex is paired or not
    paired_vertex_ids = set()
    
    for file_idx, edge_list in network_edges.items():
        for edge in edge_list:
            # Handle both index-based format (during pairing) and restored format (after JSON restoration)
            current_file_idx = edge.get('current_file_idx', edge.get('ref_file_idx'))
            current_feat_idx = edge.get('current_feature_idx', edge.get('current_file', {}).get('feature_idx'))
            paired_vertex_ids.add((current_file_idx, current_feat_idx))
            
            match_file_idx = edge.get('matched_file_idx', edge.get('match_file_idx'))
            match_feat_idx = edge.get('matched_feature_idx', edge.get('matched_file', {}).get('feature_idx'))
            paired_vertex_ids.add((match_file_idx, match_feat_idx))
    
    # Find all unpaired vertices by comparing with all_feature_data_lists
    unpaired_vertices = []
    
    for file_idx, features in enumerate(all_feature_data_lists):
        for feat_idx, feature in enumerate(features):
            vertex_id = (file_idx, feat_idx)
            
            # Include if this vertex does NOT appear in any edge in network_edges
            if vertex_id not in paired_vertex_ids:
                unpaired_vertex = {
                    'vertex_id': f"{feature.get('filename', f'file_{file_idx}')}_{feat_idx}",
                    'filename': feature.get('filename', f'file_{file_idx}'),
                    'feature_idx': feat_idx,
                    'file_idx': file_idx,
                    'x_center': float(feature.get('x_center', 0)),
                    'y_center': float(feature.get('y_center', 0)),
                    'mz_start': float(feature.get('mz_start', 0)),
                    'mz_end': float(feature.get('mz_end', 0)),
                    'rt_start': float(feature.get('rt_start', 0)),
                    'rt_end': float(feature.get('rt_end', 0)),
                    'charge': int(feature.get('charge', 0)),
                    'intensity': float(feature.get('intensity', 0)),
                    'pep_ident': feature.get('pep_ident', []),
                    'prot_ident': feature.get('prot_ident', []),
                    'openms_fid': feature.get('openms_fid', ''),
                    'ms2_scans': feature.get('ms2_scans', '')
                }
                unpaired_vertices.append(unpaired_vertex)
    
    # Save to JSON (only if output_path is provided)
    if output_path and unpaired_vertices:
        try:
            with open(output_path, 'w') as f:
                json.dump({
                    'total_unpaired': len(unpaired_vertices),
                    'unpaired_vertices': unpaired_vertices
                }, f, indent=2)
            print(f"\n[PROGRESS] Saved {len(unpaired_vertices)} unpaired vertices to: {os.path.basename(output_path)}")
            print(f"  ℹ Unpaired vertices are features that appear in all_feature_data but have no edges in network_edges")
        except Exception as e:
            print(f"✗ Warning: Could not save unpaired vertices to {output_path}: {e}")
    elif not unpaired_vertices:
        print(f"\n[PROGRESS] No unpaired vertices found (all features have at least one edge)")
    
    return unpaired_vertices


def feature_pairing_optimized(all_feature_data_lists: List[List[Dict]], best_match_only: bool = False, skip_match_building: bool = False, normalize_coordinates: bool = True, distance_calc_before_scaling: bool = False, rt_correction_json: str = None, rt_correction_mode: str = 'none', rt_source: str = 'rt_start', mz_rt_weight_ratio: float = 1.0, mz_cutoff: float = None, rt_cutoff: float = None, edges_cutoff: float = None) -> tuple:
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
        distance_calc_before_scaling: If True, use original coordinates for distance calculation
        rt_correction_json: Path to JSON file with fitted RT correction models (optional)
        rt_correction_mode: 'none', 'pre-kdtree', or 'per-distance'
            - 'none': No RT correction applied
            - 'pre-kdtree': Apply correction before building KD-tree for each pair (Option A)
            - 'per-distance': Apply correction during distance calculation (Option B)
        rt_source: Which RT field to use for correction ('rt_start', 'rt_center', 'rt_geo_center')
        mz_rt_weight_ratio: Weight ratio for RT relative to m/z (default: 1.0 for equal weighting)
            - < 1.0: RT weighs less than m/z (m/z more important)
            - 1.0: Equal weighting (default)
            - > 1.0: RT weighs more than m/z (RT more important)
    """
    if not all_feature_data_lists:
        return [], {}, {}
    
    # Build file index → filename mapping for memory optimization
    file_index_map = {}
    for file_idx, features in enumerate(all_feature_data_lists):
        if features and len(features) > 0:
            file_index_map[file_idx] = features[0].get('filename', f'file_{file_idx}')
    
    if len(all_feature_data_lists) == 1:
        matches = []
        edges = {0: []}  # Single file has no cross-file edges
        for feature in all_feature_data_lists[0]:
            match_group = feature.copy()
            match_group['files'] = [feature['filename']]
            match_group['match_score'] = 1.0
            matches.append(match_group)
        return matches, edges, file_index_map
    
    print("\n[PROGRESS] Starting feature pairing (optimized)")
    print("  Pairing features (OPTIMIZED) using KD-tree nearest-neighbor...")
    print("  Splitting features by charge for paired matching...")
    
    # Load RT correction models from JSON if provided
    rt_corrections = {}
    print(f"\n[PROGRESS] Step 1: Loading RT corrections...")
    
    if rt_correction_mode and rt_correction_mode != 'none':
        json_str = str(rt_correction_json).strip() if rt_correction_json else ""
        
        try:
            if json_str and os.path.exists(json_str):
                print(f"  Loading RT corrections from: {json_str}")
                rt_corrections = load_rt_corrections_json(json_str)
                print(f"  ✓ Loaded {len(rt_corrections)} RT correction pairs")
                print(f"  ✓ Mode: {rt_correction_mode} | RT source: {rt_source}")
                
                # Log correction keys (already using real filenames)
                if len(rt_corrections) > 0:
                    print(f"\n  Sample correction keys:")
                    for (ref_key, query_key) in list(rt_corrections.keys())[:3]:
                        print(f"    {repr(ref_key)} → {repr(query_key)}")
                    if len(rt_corrections) > 3:
                        print(f"    ... and {len(rt_corrections) - 3} more pairs")
            else:
                print(f"  ⚠ No RT correction JSON file provided - proceeding without corrections")
        except Exception as e:
            print(f"  ✗ Error loading RT corrections: {e}")
            raise ValueError(f"Failed to load RT corrections: {e}")
    else:
        print(f"  ⊘ RT correction disabled (mode={repr(rt_correction_mode)})")
    
    # Calculate scaling factors for coordinate normalization if enabled
    scale_x = 1.0
    scale_y = 1.0
    min_x = 0.0
    min_y = 0.0
    
    if normalize_coordinates:
        print("\n[PROGRESS] Step 2: Computing coordinate normalization factors...")
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
            
            # Scale both axes to [0, 1] range, then apply optional RT/m/z weighting
            scale_x = 1.0 / range_x if range_x > 0 else 1.0
            scale_y = 1.0 / range_y if range_y > 0 else 1.0
            
            # Apply mz_rt_weight_ratio to emphasize one dimension over the other
            # ratio > 1.0 means RT weighs more, < 1.0 means m/z weighs more
            scale_y_weighted = scale_y * mz_rt_weight_ratio
            
            print(f"    X range: {min_x:.2f} to {max_x:.2f} (range: {range_x:.2f})")
            print(f"    Y range: {min_y:.2f} to {max_y:.2f} (range: {range_y:.2f})")
            print(f"    Scale factors - X (m/z): {scale_x:.6f}, Y (RT): {scale_y:.6f}")
            if mz_rt_weight_ratio != 1.0:
                print(f"    Dimension weighting - m/z:RT ratio: 1.0:{mz_rt_weight_ratio:.4f}")
                print(f"    Weighted scale factor - Y (RT): {scale_y_weighted:.6f}")
    
    def scale_coords(coords):
        """Scale coordinates using computed normalization factors and optional dimension weighting."""
        if normalize_coordinates and len(coords) > 0:
            scaled = coords.copy().astype(float)
            scaled[:, 0] = (scaled[:, 0] - min_x) * scale_x
            scaled[:, 1] = (scaled[:, 1] - min_y) * scale_y_weighted
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
    
    # Apply JSON-based RT corrections if available
    print(f"\n[PROGRESS] Step 2.5: RT Correction Application Status:")
    print(f"  - rt_correction_mode: {repr(rt_correction_mode)}")
    print(f"  - rt_corrections loaded: {len(rt_corrections)} pairs")
    if len(rt_corrections) > 0:
        print(f"  - Available correction keys: {list(rt_corrections.keys())[:5]}")  # Show first 5
        # Show details of first correction
        first_key = list(rt_corrections.keys())[0]
        first_correction = rt_corrections[first_key]
        print(f"  - First correction key: {repr(first_key)}")
        print(f"    - method: {first_correction.get('method', 'unknown')}")
        print(f"    - has serialization: {'serialization' in first_correction}")
    
    if len(rt_corrections) > 0 and rt_correction_mode == 'pre-kdtree':
        print("\n[PROGRESS] Step 2.6: Applying JSON-based RT corrections to features...")
        
        # Get reference file name from first file's features
        ref_file_name = None
        if file_features_by_charge and file_features_by_charge[0]:
            first_charge_features = list(file_features_by_charge[0].values())[0]
            if first_charge_features:
                ref_file_name = first_charge_features[0].get('filename', None)
        
        if ref_file_name:
            print(f"  Reference file: {repr(ref_file_name)}")
            
            # Apply corrections to each query file
            corrections_applied = 0
            corrections_skipped = 0
            corrections_found_key = 0
            corrections_not_found_key = 0
            
            for file_idx in range(1, len(file_features_by_charge)):
                # Get query file name from first feature
                query_file_name = None
                query_features_list = list(file_features_by_charge[file_idx].values())
                if query_features_list:
                    first_query_features = query_features_list[0]
                    if first_query_features:
                        query_file_name = first_query_features[0].get('filename', None)
                
                if not query_file_name:
                    print(f"    File {file_idx}: Could not determine filename")
                    continue
                
                # Look up correction for this file pair
                correction_key = (ref_file_name, query_file_name)
                print(f"    File {file_idx}: Checking correction key: {repr(correction_key)}")
                if correction_key in rt_corrections:
                    print(f"      ✓ Correction key found! Applying...")
                    correction_dict = rt_corrections[correction_key]
                    print(f"      Method: {correction_dict.get('method', 'unknown')}")
                    
                    # Apply correction to all features in this file
                    for charge, charge_features in file_features_by_charge[file_idx].items():
                        for feature in charge_features:
                            mz_value = feature.get('x_center', 0.0)
                            rt_value = feature.get('y_center', 0.0)
                            
                            try:
                                # Get the RT difference correction
                                # Equations are RT vs RT-Diff, so input is RT value, not m/z
                                rt_diff = apply_rt_correction_json(rt_value, correction_dict)
                                # Apply correction: corrected_rt = rt - rt_diff
                                corrected_rt = rt_value - rt_diff
                                
                                # Store original RT and update with corrected value
                                if 'y_center_original' not in feature:
                                    feature['y_center_original'] = rt_value
                                feature['y_center'] = corrected_rt
                                corrections_applied += 1
                            except Exception as e:
                                print(f"      WARNING: Failed to apply correction to feature: {e}")
                                corrections_skipped += 1
                else:
                    # Try reverse direction
                    correction_key_rev = (query_file_name, ref_file_name)
                    print(f"      ✗ Correction key NOT found, trying reverse: {repr(correction_key_rev)}")
                    if correction_key_rev in rt_corrections:
                        correction_dict = rt_corrections[correction_key_rev]
                        print(f"      ✓ Reverse correction key found! Applying...")
                        print(f"      Method: {correction_dict.get('method', 'unknown')}")
                        
                        # Apply reverse correction to all features in this file
                        for charge, charge_features in file_features_by_charge[file_idx].items():
                            for feature in charge_features:
                                mz_value = feature.get('x_center', 0.0)
                                rt_value = feature.get('y_center', 0.0)
                                
                                try:
                                    # Get the RT difference correction (this is for ref vs query)
                                    # Equations are RT vs RT-Diff, so input is RT value, not m/z
                                    # So to reverse it: corrected_rt = rt + rt_diff
                                    rt_diff = apply_rt_correction_json(rt_value, correction_dict)
                                    corrected_rt = rt_value + rt_diff
                                    
                                    # Store original RT and update with corrected value
                                    if 'y_center_original' not in feature:
                                        feature['y_center_original'] = rt_value
                                    feature['y_center'] = corrected_rt
                                    corrections_applied += 1
                                except Exception as e:
                                    print(f"      WARNING: Failed to apply reverse correction to feature: {e}")
                                    corrections_skipped += 1
                    else:
                        print(f"      ✗ Neither forward nor reverse correction key found, skipping")
                        corrections_skipped += 1
            
            print(f"\n  RT Correction Application Summary:")
            print(f"  - Applied corrections to {corrections_applied} features")
            print(f"  - Skipped/not found: {corrections_skipped} files")
        else:
            print(f"  ✗ Could not determine reference file name, skipping RT corrections")
    elif len(rt_corrections) == 0:
        print(f"  ⊘ No RT corrections loaded (rt_corrections dict is empty)")
    elif rt_correction_mode != 'pre-kdtree':
        print(f"  ⊘ RT correction mode is {repr(rt_correction_mode)}, not 'pre-kdtree' - not applying corrections")
    
    # Print charge diagnostic information
    print("  Charge distribution by file:")
    for file_idx, charge_dict in enumerate(file_features_by_charge):
        print(f"    File {file_idx}: {dict(sorted([(c, len(f)) for c, f in charge_dict.items()]))}")
    
    # Build KD-trees for each charge group in each file
    # Structure: {file_idx: {charge: {'kdtree': cKDTree, 'features': [features], 'coords': np.array}}}
    print("\n[PROGRESS] Step 3: Building KD-trees for each file/charge combination...")
    get_memory_info("BEFORE KD-TREE BUILDING", print_details=True)
    
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
    
    print(f"  ✓ KD-trees built successfully")
    get_memory_info("AFTER KD-TREE BUILDING", print_details=True)
    
    # Pre-compute corrected KD-trees for all file pairs (if RT corrections enabled)
    # Note: Pre-kdtree computation is skipped for JSON-based corrections
    # JSON corrections are applied directly during feature matching
    kdtrees_corrected_by_pair = {}
    if rt_correction_mode == 'pre-kdtree' and len(rt_corrections) > 0:
        print("\n[PROGRESS] Step 3.5: Pre-computing corrected KD-trees for all file pairs...")
        print(f"  WARNING: Pre-kdtree computation is only supported for TSV-based corrections.")
        print(f"  Received JSON-based corrections ({len(rt_corrections)} pairs).")
        print(f"  JSON corrections will be applied during matching instead.")
        print(f"  ⊘ Skipping pre-kdtree computation")
    else:
        print(f"\n[PROGRESS] Step 3.5: RT correction mode analysis...")
        if rt_correction_mode != 'pre-kdtree':
            print(f"  ⊘ Pre-kdtree correction not requested (mode: {rt_correction_mode})")
        elif len(rt_corrections) == 0:
            print(f"  ⊘ No RT corrections loaded")
    
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
    
    print("\n[PROGRESS] Step 4: Performing pairwise feature matching...")
    get_memory_info("BEFORE PAIRWISE MATCHING", print_details=True)
    num_files = len(all_feature_data_lists)
    total_pairwise_comparisons = 0
    processed_pairwise_comparisons = 0
    # Calculate expected number of comparisons for progress tracking
    for file_file_idx, charge_dict in enumerate(file_features_by_charge):
        for charge, charge_features in charge_dict.items():
            for ref_feature in charge_features:
                for query_file_idx in range(len(all_feature_data_lists)):
                    if query_file_idx != file_file_idx:
                        total_pairwise_comparisons += 1
    print(f"  Expected total pairwise comparisons: {total_pairwise_comparisons}")
    print(f"  Processing {len(first_file_features)} features from first file (with charge filtering)...")
    
    for file_file_idx, charge_dict in enumerate(file_features_by_charge):
        edges[file_file_idx] = []
        
        for charge, charge_features in charge_dict.items():
            for feature_idx, ref_feature in enumerate(charge_features):
                if feature_idx % 100 == 0:
                    print(f"    Processing file {file_file_idx}, charge {charge}: feature {feature_idx}/{len(charge_features)}", end='\r', flush=True)
                
                # Use scaled coordinates for KD-tree query
                ref_coords_scaled = scale_coords(np.array([[ref_feature['x_center'], ref_feature['y_center']]]))
                
                # [MEMORY OPTIMIZATION] Store file/feature indices instead of filename strings
                current_file_idx = ref_feature['_file_idx']
                current_feature_idx = ref_feature['_feat_idx']
                
                # Query each other file's KD-tree for matches within distance cutoff
                # ONLY match with features of the SAME charge
                for query_file_idx in range(len(all_feature_data_lists)):
                    processed_pairwise_comparisons += 1
                    if processed_pairwise_comparisons % max(1, total_pairwise_comparisons // 20) == 0:
                        pct = int(100 * processed_pairwise_comparisons / total_pairwise_comparisons) if total_pairwise_comparisons > 0 else 0
                        print(f"    [Pairwise comparisons: {processed_pairwise_comparisons}/{total_pairwise_comparisons} ({pct}%)]", end='\r', flush=True)
                    
                    if query_file_idx == file_file_idx:
                        continue
                    
                    # [CRITICAL CHECK] Skip if KD-trees for this file have been freed (set to None)
                    # This can happen if we're comparing against a file that was already processed
                    if kdtrees_by_charge[query_file_idx] is None:
                        continue
                    
                    # Get the KD-tree for this charge in the query file
                    if charge not in kdtrees_by_charge[query_file_idx]:
                        # No features with this charge in the query file, skip
                        continue
                    
                    # Use pre-computed corrected KD-tree if available, otherwise use uncorrected
                    if (rt_correction_mode == 'pre-kdtree' and 
                        file_file_idx in kdtrees_corrected_by_pair and 
                        query_file_idx in kdtrees_corrected_by_pair[file_file_idx] and
                        charge in kdtrees_corrected_by_pair[file_file_idx][query_file_idx]):
                        # Use pre-computed corrected KD-tree
                        kdtree = kdtrees_corrected_by_pair[file_file_idx][query_file_idx][charge]['kdtree']
                        query_features = kdtrees_corrected_by_pair[file_file_idx][query_file_idx][charge]['features']
                        kdtree_coords = kdtrees_corrected_by_pair[file_file_idx][query_file_idx][charge]['coords']
                    else:
                        # Use uncorrected KD-tree
                        kdtree = kdtrees_by_charge[query_file_idx][charge]['kdtree']
                        query_features = kdtrees_by_charge[query_file_idx][charge]['features']
                        kdtree_coords = kdtrees_by_charge[query_file_idx][charge]['coords']
                    
                    # Query only the single best match for this feature from the other file
                        # (best_match_only is now hardcoded to True)
                        distance, idx = kdtree.query(ref_coords_scaled, k=1)
                        
                        if isinstance(distance, np.ndarray):
                            distance = distance.flat[0]
                        if isinstance(idx, np.ndarray):
                            idx = idx.flat[0]
                        
                        if idx >= 0:
                            matched_feature = query_features[int(idx)]
                            
                            # Calculate distance with optional RT correction
                            if distance_calc_before_scaling:
                                # Use original coordinates
                                mz_dist = (ref_feature['x_center'] - matched_feature['x_center']) ** 2
                                rt_val_query = matched_feature['y_center']
                                
                                # Note: Per-distance RT correction was previously TSV-based only
                                # JSON corrections are applied during pre-kdtree or static mode
                                
                                rt_dist = (ref_feature['y_center'] - rt_val_query) ** 2
                            else:
                                # Use scaled coordinates (default)
                                matched_scaled_coords = kdtree_coords[int(idx)]
                                ref_scaled_coords = ref_coords_scaled[0]
                                mz_dist = (ref_scaled_coords[0] - matched_scaled_coords[0]) ** 2
                                rt_dist = (ref_scaled_coords[1] - matched_scaled_coords[1]) ** 2
                            
                            # Calculate distances in SCALED space (for KD-tree efficiency/visualization)
                            exact_distance = float(np.sqrt(mz_dist + rt_dist))
                            mz_distance = float(np.sqrt(mz_dist))  # Separate m/z distance in scaled space
                            rt_distance = float(np.sqrt(rt_dist))  # Separate RT distance in scaled space
                            
                            # [CRITICAL FIX] Calculate distances in ORIGINAL coordinate space for cutoff comparison
                            # This ensures cutoff parameters (which are in original space) are compared fairly
                            original_mz_dist = (ref_feature['x_center'] - matched_feature['x_center']) ** 2
                            original_rt_dist = (ref_feature['y_center'] - matched_feature['y_center']) ** 2
                            original_mz_distance = float(np.sqrt(original_mz_dist))
                            original_rt_distance = float(np.sqrt(original_rt_dist))
                            original_exact_distance = float(np.sqrt(original_mz_dist + original_rt_dist))
                            
                            # [MEMORY OPTIMIZATION] Filter edges at creation time to avoid storing rejected edges
                            # Apply cutoff filters DURING loop instead of post-processing
                            # ALL COMPARISONS USE ORIGINAL-SPACE DISTANCES FOR CONSISTENCY
                            if mz_cutoff is not None and rt_cutoff is not None:
                                # Coordinate-based filtering: check both m/z and RT cutoffs (in original space)
                                if original_mz_distance > mz_cutoff or original_rt_distance > rt_cutoff:
                                    continue  # Skip this edge, don't add to dict
                            elif edges_cutoff is not None:
                                # Distance-based filtering: check euclidean distance cutoff (in original space)
                                if original_exact_distance > edges_cutoff:
                                    continue  # Skip this edge, don't add to dict
                            
                            # [MEMORY OPTIMIZATION] Store indices instead of filename dicts
                            matched_feature_idx = matched_feature['_feat_idx']
                            
                            # Add edge connecting the two features (indices only, no filenames)
                            edge = {
                                'current_feature_idx': current_feature_idx,
                                'matched_feature_idx': matched_feature_idx,
                                'distance': exact_distance,
                                'ref_file_idx': file_file_idx,      # Reference file being looped over
                                'match_file_idx': query_file_idx,   # Query file for matching
                                'charge': charge  # Include charge validation in edge
                            }
                            edges[file_file_idx].append(edge)
        
        # [MEMORY OPTIMIZATION] After completing all charges/features for this file, force garbage collection
        # This prevents edges[0], edges[1], ... edges[N] from all accumulating simultaneously
        # Reduces memory bloat as the loop processes multiple files sequentially
        print(f"    File {file_file_idx} complete: {len(edges[file_file_idx])} edges. Collecting garbage...")
        gc.collect()
                            
    
    # Diagnostic output: verify all edges have matching charges
    print("\n  Verifying charge consistency in edges:")
    get_memory_info("AFTER PAIRWISE MATCHING (before consolidation)", print_details=True)
    
    # Print structure sizes before consolidation
    print_structure_sizes("STRUCTURE SIZES AFTER MATCHING", {
        'edges (total)': edges,
        'file_features_by_charge': file_features_by_charge,
        'all_feature_data_lists': all_feature_data_lists
    })
    
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
        
        # [AGGRESSIVE CLEANUP] Free large intermediate data structures before returning
        print("  ✓ Performing aggressive memory cleanup...")
        mem_before_cleanup = get_memory_info("BEFORE CLEANUP", print_details=True)
        
        # Delete large data structures that are no longer needed
        del file_features_by_charge
        del kdtrees_by_charge
        del scale_coords
        gc.collect()
        
        mem_after_cleanup = get_memory_info("AFTER CLEANUP", print_details=True)
        print(f"    Memory freed: {mem_before_cleanup - mem_after_cleanup:.1f} MB")
        
        return matches, edges, file_index_map
    
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
    
    print(f"\n[PROGRESS] Step 4 complete: Pairwise matching done ({processed_pairwise_comparisons} comparisons)")
    print(f"[PROGRESS] Feature pairing completed successfully")
    
    # [DELETION POINT #3] rt_corrections loaded from JSON is no longer used after pairing
    if 'rt_corrections' in locals():
        del rt_corrections
    
    return matches, edges, file_index_map


def restore_filenames_in_edges_json(json_path: str, file_index_map: Dict[int, str]) -> None:
    """
    Post-process edges JSON to restore filenames from indices and clean up edge structure.
    During pairing, filenames are replaced with indices to save memory (~1.26GB).
    This function restores filenames and removes index fields, creating clean nested structures.
    
    Converts from:
    - ref_file_idx, match_file_idx, current_feature_idx, matched_feature_idx
    
    To clean nested structure for build_network_graph.py:
    - current_file: {filename, feature_idx}
    - matched_file: {filename, feature_idx}
    
    Args:
        json_path: Path to JSON file with edges using indices
        file_index_map: Mapping {file_idx: filename}
    """
    if not os.path.exists(json_path):
        print(f"⊘ Edges JSON not found: {json_path}")
        return
    
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        print(f"Restoring filenames in edges JSON and cleaning structure...")
        edges_updated = 0
        edges_missing_filenames = 0
        
        # Handle both direct edges dict and wrapped format
        if isinstance(data, dict) and 'network_edges' in data:
            # Wrapped format from edges save
            edges_dict = data['network_edges']
        else:
            # Direct edges dict
            edges_dict = data
        
        # Restore filenames in edges and remove index fields
        for file_idx_str, edge_list in edges_dict.items():
            for edge in edge_list:
                # Build current_file nested structure
                ref_file_idx = edge.get('ref_file_idx')
                if ref_file_idx is not None and ref_file_idx in file_index_map:
                    edge['current_file'] = {
                        'filename': file_index_map[ref_file_idx],
                        'feature_idx': edge.get('current_feature_idx', -1)
                    }
                else:
                    edges_missing_filenames += 1
                
                # Build matched_file nested structure
                match_file_idx = edge.get('match_file_idx')
                if match_file_idx is not None and match_file_idx in file_index_map:
                    edge['matched_file'] = {
                        'filename': file_index_map[match_file_idx],
                        'feature_idx': edge.get('matched_feature_idx', -1)
                    }
                else:
                    edges_missing_filenames += 1
                
                # Remove index fields (keep only nested structures)
                for key in ['ref_file_idx', 'match_file_idx', 'current_feature_idx', 'matched_feature_idx']:
                    if key in edge:
                        del edge[key]
                
                edges_updated += 1
        
        # Save cleaned JSON
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"  ✓ Restored and cleaned {edges_updated} edges in {os.path.basename(json_path)}")
        if edges_missing_filenames > 0:
            print(f"  ⚠ {edges_missing_filenames} edge endpoints missing from file_index_map (may cause issues downstream)")
    
    except Exception as e:
        print(f"✗ Warning: Failed to restore filenames in JSON: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Pair matching features across multiple files"
    )
    parser.add_argument("--input_dir", help="Directory containing feature pickle files (default: use current dir if --input_files not provided)")
    parser.add_argument("--input_files", nargs='*', help="Explicit list of feature pickle files to process (alternative to --input_dir)")
    parser.add_argument("--output_pkl", required=True, help="Output pickle file path")
    parser.add_argument("--output_json", required=True, help="Output JSON file path")
    parser.add_argument("--output_edges_pkl", help="Output edges dictionary pickle file (default: inferred from output_pkl)")
    parser.add_argument("--output_unpaired_json", help="Output JSON file for unpaired vertices (features with no edges)")

    parser.add_argument("--skip-matchfinder", action='store_true', help="Skip matchfinder processing (match consolidation). Only export edges dictionary. Useful for network graph analysis.")
    parser.add_argument("--match_cutoff", type=float, default=None, help="Minimum match score cutoff (0.0-1.0)")
    # Note: best_match_only is now hardcoded to True
    parser.add_argument("--edges_cutoff", type=float, default=None, help="Maximum distance cutoff for edges (filters edges based on euclidean distance)")
    parser.add_argument("--mz_cutoff", type=float, default=None, help="Maximum m/z (x-axis) coordinate difference for edge filtering")
    parser.add_argument("--rt_cutoff", type=float, default=None, help="Maximum RT (y-axis) coordinate difference for edge filtering")
    parser.add_argument("--distance_calc_before_scaling", action='store_true', help="Calculate distance before coordinate scaling (default: after scaling)")
    parser.add_argument("--normalize_coordinates", type=lambda x: x.lower() in ('true', '1', 'yes'), nargs='?', const=True, default=None, help="Enable/disable coordinate normalization before KD-tree matching (default: auto - OFF for mz/rt cutoff, ON for euclidean cutoff)")
    parser.add_argument("--postpair_normalize_coordinates", action='store_true', help="Rescale axes to balance coordinate ranges and recalculate distances (adds rescaled distance fields; original distances preserved)")
    parser.add_argument("--normalize_edge_distances", action='store_true', help="Normalize edge distances to [0, 1] range for visualization (default: disabled)")
    parser.add_argument("--analyze_pep_idents", action='store_true', help="Perform detailed pep_ident matching analysis (default: enabled)")
    parser.add_argument("--skip_json_output", action='store_true', help="Skip JSON serialization for speed (default: disabled - JSON output enabled)")
    parser.add_argument("--rt_correction_json", type=str, default=None, help="Path to JSON file with fitted RT correction models (optional)")
    parser.add_argument("--trafoxmls_dir", type=str, default=None, help="Directory containing trafoXML files for AlignTree method (OpenMS RT transformations)")
    parser.add_argument("--rt_alignment_method", type=str, default='polynomial', help="RT alignment method: 'polynomial', 'spline', 'loess', 'rbf', 'piecewise', 'aligntree' (determines whether to use trafoXML)")
    parser.add_argument("--rt_correction_mode", type=str, choices=['none', 'pre-kdtree', 'per-distance'], default='none', help="RT correction mode: 'none' (default): no correction, 'pre-kdtree': correct before KD-tree building, 'per-distance': correct during distance calculation")
    parser.add_argument("--rt_source", type=str, choices=['y_center', 'y_center_geo', 'rt_start', 'rt_end'], default='y_center', help="Which RT field to use for RT correction (default: y_center)")
    parser.add_argument("--mz_rt_weight_ratio", type=float, default=1.0, help="Weight ratio for RT relative to m/z in distance calculations (default: 1.0 for equal weighting). Values > 1.0 emphasize RT, < 1.0 emphasize m/z")
    parser.add_argument("--disable_multiprocessing", action='store_true', help="Disable multiprocessing (default: enabled)")
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print(f"Feature Pairing")
    print(f"{'='*70}")
    
    # Log RT correction configuration received from command line
    print(f"\nRT Correction Configuration (from command-line args):")
    print(f"  - rt_correction_mode: {repr(args.rt_correction_mode)}")
    print(f"  - rt_correction_json: {repr(args.rt_correction_json)}")
    print(f"  - rt_source: {repr(args.rt_source)}")
    json_exists = os.path.exists(args.rt_correction_json) if args.rt_correction_json else False
    print(f"  - JSON file exists: {json_exists}")
    
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
    
    # Capture the TRUE "before filtering" count - total features from TSV files BEFORE any cutoffs
    total_features_before_filtering = sum(len(feature_list) for feature_list in all_feature_data_lists)
    print(f"Total features loaded from TSV files: {total_features_before_filtering}")
    
    
    # [DIAGNOSTIC] Print memory state after loading all features
    get_memory_info("AFTER LOADING ALL FEATURES", print_details=True)
    print_structure_sizes("INITIAL STRUCTURE SIZES", {
        'all_feature_data_lists': all_feature_data_lists
    })
    
    # [DELETION POINT #1] feature_files no longer used after this point
    # Save the count before deleting since it's needed in output_data
    num_input_files = len(feature_files)
    del feature_files
    gc.collect()
    
    # DEBUG: Print command-line arguments related to trafoXML and RT correction
    print(f"\n[DEBUG] Command-line argument values:")
    print(f"  - args.trafoxmls_dir: {repr(args.trafoxmls_dir)}")
    print(f"  - args.rt_correction_json: {repr(args.rt_correction_json)}")
    print(f"  - args.rt_correction_mode: {repr(args.rt_correction_mode)}")
    
    # Validate: aligntree method requires trafoXML directory
    if args.rt_alignment_method == 'aligntree' and not args.trafoxmls_dir:
        print(f"\n✗ ERROR: RT alignment method 'aligntree' requires a trafoXML directory")
        print(f"  You selected: --rt_alignment_method aligntree")
        print(f"  But did not provide: --trafoxmls_dir <path>")
        print(f"\n  AlignTree uses pre-calculated RT transformations from trafoXML files.")
        print(f"  You must provide the directory containing trafoXML files from OpenMS RT alignment.")
        print(f"\n  Example usage:")
        print(f"    python pair_features.py \\")
        print(f"      --input_files ... \\")
        print(f"      --output_pkl ... \\")
        print(f"      --output_json ... \\")
        print(f"      --rt_alignment_method aligntree \\")
        print(f"      --trafoxmls_dir <path/to/trafoXML/files>")
        sys.exit(1)
    
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
    
    # Apply AlignTree RT transformations if trafoXML directory is provided AND aligntree method is selected
    # NOTE: Only use trafoXML if rt_alignment_method is 'aligntree'
    if args.trafoxmls_dir and args.rt_alignment_method == 'aligntree':
        print(f"\n[PROGRESS] Loading AlignTree RT transformations from: {args.trafoxmls_dir}")
        print(f"  (Using trafoXML-based RT corrections, not JSON-based fitting)")
        try:
            # Build feature ID -> transformed RT mapping
            fid_to_rt_dict = build_trafoxmls_fid_dict(args.trafoxmls_dir)
            
            # Apply transformations to all features
            print(f"[PROGRESS] Applying RT transformations to features...")
            all_feature_data_lists = apply_aligntree_rt_transformations(all_feature_data_lists, fid_to_rt_dict)
            print(f"  ✓ AlignTree transformations applied successfully")
            
            # [DELETION POINT #2] fid_to_rt_dict no longer used after this point
            del fid_to_rt_dict
            gc.collect()
            
            # Set rt_correction_mode to 'none' since corrections already applied via transformations
            args.rt_correction_mode = 'none'
        except Exception as e:
            print(f"  ✗ ERROR applying AlignTree transformations: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    elif args.trafoxmls_dir and args.rt_alignment_method != 'aligntree':
        # TrafoXML directory provided but not using aligntree method - warn and skip
        print(f"\n[WARNING] TrafoXML directory provided but rt_alignment_method is '{args.rt_alignment_method}' (not 'aligntree')")
        print(f"  Skipping trafoXML transformations - will use {args.rt_alignment_method} method instead")
    elif args.rt_correction_mode and args.rt_correction_mode != 'none':
        # JSON-based corrections will be applied by feature_pairing_optimized
        pass
    
    # Pair features using KD-tree based optimized method (always skip component building)
    paired_features, network_edges, file_index_map = feature_pairing_optimized(all_feature_data_lists, best_match_only=True, skip_match_building=True, normalize_coordinates=normalize_coordinates, distance_calc_before_scaling=args.distance_calc_before_scaling, rt_correction_json=args.rt_correction_json, rt_correction_mode=args.rt_correction_mode, rt_source=args.rt_source, mz_rt_weight_ratio=args.mz_rt_weight_ratio, mz_cutoff=args.mz_cutoff, rt_cutoff=args.rt_cutoff, edges_cutoff=args.edges_cutoff)
    
    # [DIAGNOSTIC] Print memory state after pairing completes
    get_memory_info("AFTER FEATURE PAIRING", print_details=True)
    print_structure_sizes("STRUCTURE SIZES AFTER PAIRING", {
        'network_edges': network_edges,
        'all_feature_data_lists': all_feature_data_lists
    })
    
    print(f"Built network graph with {sum(len(e) for e in network_edges.values())} edges")
    print(f"  File index mapping: {file_index_map}")
    print(f"  [Note: Edges were filtered during pairing based on cutoff parameters]")
    
    # Count unique vertices (features) in network_edges after cutoff filtering
    # Each feature is identified by (file_idx, feature_idx) tuple
    unique_vertices = set()
    for ref_file_idx, edge_list in network_edges.items():
        for edge in edge_list:
            # Add reference feature vertex
            ref_vertex = (edge['ref_file_idx'], edge['current_feature_idx'])
            unique_vertices.add(ref_vertex)
            # Add matched feature vertex
            match_vertex = (edge['match_file_idx'], edge['matched_feature_idx'])
            unique_vertices.add(match_vertex)
    
    total_features_after_filter = len(unique_vertices)
    print(f"Total unique features (vertices) in network edges after filtering: {total_features_after_filter}")
    # Filter edges by cutoff if specified
    # [OPTIMIZATION] Filtering now happens during Step 4 pairing, not as post-processing
    if args.mz_cutoff is not None and args.rt_cutoff is not None:
        print(f"✓ Edges were filtered during pairing with coordinate cutoffs: mz={args.mz_cutoff}, rt={args.rt_cutoff}")
    elif args.edges_cutoff is not None:
        print(f"✓ Edges were filtered during pairing with distance cutoff: {args.edges_cutoff}")
    else:
        print(f"ℹ No edge filtering (no cutoff parameters specified)")

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
    
    # Extract and save unpaired vertices (features with no edges in network_edges)
    unpaired_vertices = extract_and_save_unpaired_vertices(network_edges, all_feature_data_lists, args.output_unpaired_json)
    
    # Build combined pep_ident and prot_ident lookup from all feature data for export to build_network_graph.py
    # IMPORTANT: Build this before deleting all_feature_data_lists to avoid losing feature metadata
    print("\nBuilding combined pep_ident and prot_ident lookup from all features...")
    ident_lookup = {}
    for file_idx, file_features in enumerate(all_feature_data_lists):
        for feat_idx, feature in enumerate(file_features):
            if isinstance(feature, dict):
                filename = feature.get('filename', f'file_{file_idx}')
                pep_idents = feature.get('pep_ident', [])
                prot_idents = feature.get('prot_ident', [])
                x_center = feature.get('x_center', 0.0)
                y_center = feature.get('y_center', 0.0)
                intensity = feature.get('intensity', 0.0)
                openms_fid = feature.get('openms_fid', '')
                ms2_scans = feature.get('ms2_scans', '')
                charge = feature.get('charge', 0)
                
                # Create lookup key: "filename_feature_idx"
                lookup_key = f"{filename}_{feat_idx}"
                ident_lookup[lookup_key] = {
                    'pep_ident': pep_idents if pep_idents else [],
                    'prot_ident': prot_idents if prot_idents else [],
                    'x_center': x_center,
                    'y_center': y_center,
                    'intensity': intensity,
                    'openms_fid': openms_fid,
                    'ms2_scans': ms2_scans,
                    'charge': charge
                }
    
    print(f"Extracted identifications for {len(ident_lookup)} features")
    
    # [DELETION POINT] all_feature_data_lists no longer needed - free memory before saving
    print("Clearing feature data from memory...")
    del all_feature_data_lists
    del unpaired_vertices
    gc.collect()
    print("  ✓ Feature data freed (~500MB)")
    
    # [DIAGNOSTIC] Print memory state after clearing features
    get_memory_info("AFTER CLEARING FEATURES", print_details=True)
    print_structure_sizes("STRUCTURE SIZES AFTER FEATURE CLEANUP", {
        'network_edges': network_edges,
        'ident_lookup': ident_lookup
    })

    # Save results with edges and ident_lookup only (no paired_features)
    # Save results with edges and ident_lookup only (no paired_features)
    print("[DEBUG] Creating output_data dictionary...")
    mem_before_output_create = get_memory_info("BEFORE OUTPUT CREATION", print_details=True)
    
    output_data = {
        'network_edges': network_edges,  # Include network graph structure
        'ident_lookup': ident_lookup,  # Combined lookup: {filename_feat_idx: {pep_ident: [...], prot_ident: [...]}}
        'total_features_before_filtering': total_features_before_filtering,  # True count from TSV files before cutoffs
        'total_features_after_filtering': total_features_after_filter,  # Count of unique features in edges after cutoff filtering
        'cutoff_params': {
            'mz_cutoff': args.mz_cutoff,
            'rt_cutoff': args.rt_cutoff,
            'edges_cutoff': args.edges_cutoff
        },
        'pairing_parameters': {
            'pairing_method': 'kdtree',
            'best_match_only': True,
            'distance_calc_before_scaling': args.distance_calc_before_scaling,
            'normalize_coordinates': normalize_coordinates,
            'postpair_normalize_coordinates': args.postpair_normalize_coordinates,
            'postpair_note': 'If postpair_normalize_coordinates=true: distance is original euclidean distance; distance_rescaled is from rescaled coordinate space',
            'match_cutoff': args.match_cutoff,
            'input_files': num_input_files,
            'skip_matchfinder': False  # Always False now (always skip component building)
        }
    }
    
    mem_after_output_create = get_memory_info("AFTER OUTPUT CREATION", print_details=True)
    print(f"[DEBUG] Memory delta for output_data creation: {mem_after_output_create - mem_before_output_create:.1f} MB")
    
    # OPTIMIZED: Save pickle (fast binary format) - NO paired_features, only edges and ident_lookup
    print("[DEBUG] Starting pickle save...")
    mem_before_pickle = get_memory_info("BEFORE PICKLE SAVE", print_details=True)
    print("Saving output (pickle) - edges and ident_lookup only...")
    with open(args.output_pkl, 'wb') as f:
        pickle.dump(output_data, f)
    mem_after_pickle = get_memory_info("AFTER PICKLE SAVE", print_details=True)
    print(f"[DEBUG] Memory delta for pickle save: {mem_after_pickle - mem_before_pickle:.1f} MB")
    print(f"✓ Saved edges and ident_lookup (pickle): {os.path.basename(args.output_pkl)}")
    
    # OPTIMIZED: Skip JSON by default (optional) - JSON serialization of 45K+ features is slow
    if not args.skip_json_output:
        print("[DEBUG] Starting JSON save...")
        mem_before_json = get_memory_info("BEFORE JSON SAVE", print_details=True)
        print("Saving output (JSON)...  (this may take a while, use --skip_json_output to skip)")
        with open(args.output_json, 'w') as f:
            json.dump(output_data, f, indent=2)
        mem_after_json = get_memory_info("AFTER JSON SAVE", print_details=True)
        print(f"[DEBUG] Memory delta for JSON save: {mem_after_json - mem_before_json:.1f} MB")
        print(f"✓ Saved edges and ident_lookup (JSON): {os.path.basename(args.output_json)}")
    else:
        print(f"⊘ Skipped JSON output (--skip_json_output enabled)")
    
    # Save combined ident lookup (pep_ident + prot_ident) separately for fast loading in build_network_graph.py
    print("[DEBUG] Starting ident_lookup save...")
    mem_before_ident = get_memory_info("BEFORE IDENT_LOOKUP SAVE", print_details=True)
    ident_lookup_path = args.output_json.replace('.json', '_ident_lookup.json')
    with open(ident_lookup_path, 'w') as f:
        json.dump(ident_lookup, f, indent=2)
    mem_after_ident = get_memory_info("AFTER IDENT_LOOKUP SAVE", print_details=True)
    print(f"[DEBUG] Memory delta for ident_lookup save: {mem_after_ident - mem_before_ident:.1f} MB")
    print(f"✓ Saved combined ident lookup: {os.path.basename(ident_lookup_path)}")
    
    # Save edges dictionary separately for direct use with build_network_graph.py
    print("[DEBUG] Starting edges save...")
    mem_before_edges = get_memory_info("BEFORE EDGES SAVE", print_details=True)
    edges_output_path = args.output_edges_pkl
    if not edges_output_path:
        # Infer from output_pkl: replace 'paired_features' with 'edges' or add '_edges' suffix
        edges_output_path = str(args.output_pkl).replace('paired_features', 'edges')
        if edges_output_path == args.output_pkl:
            # If replacement didn't work, add suffix
            edges_output_path = args.output_pkl.replace('.pkl', '_edges.pkl')
    
    with open(edges_output_path, 'wb') as f:
        pickle.dump(network_edges, f)
    mem_after_edges_pkl = get_memory_info("AFTER EDGES PKL SAVE", print_details=True)
    print(f"[DEBUG] Memory delta for edges pkl save: {mem_after_edges_pkl - mem_before_edges:.1f} MB")
    print(f"✓ Saved edges dictionary (pickle): {os.path.basename(edges_output_path)}")
    print(f"  Optionally use with: python build_network_graph.py --input_pkl {os.path.basename(edges_output_path)}")
    
    # Also save edges JSON for reporting and visualization pipelines
    print("[DEBUG] Starting edges JSON save...")
    edges_json_path = edges_output_path.replace('.pkl', '.json')
    if not args.skip_json_output:
        try:
            # DIAGNOSTIC: Add a marker to the first edge so we know this code is running
            if network_edges:
                first_file_idx = list(network_edges.keys())[0]
                if network_edges[first_file_idx]:
                    # Add diagnostic field to FIRST EDGE ONLY
                    network_edges[first_file_idx][0]['__DIAGNOSTIC_CODE_EXECUTED__'] = True
            
            # Convert network_edges dict to JSON-serializable format
            edges_for_json = {}
            for edge_key, edge_data in network_edges.items():
                # Convert tuple keys to strings for JSON serialization
                str_key = str(edge_key)
                edges_for_json[str_key] = edge_data
            
            # Replace file indices with filenames FOR JSON OUTPUT ONLY
            # (keeping in-memory dict small with indices)
            print("[DEBUG] Converting file indices to filenames during JSON write...")
            edges_converted = 0
            for file_idx_str, edge_list in edges_for_json.items():
                for edge in edge_list:
                    # Build current_file nested structure from index
                    ref_file_idx = edge.get('ref_file_idx')
                    if ref_file_idx is not None and ref_file_idx in file_index_map:
                        edge['current_file'] = {
                            'filename': file_index_map[ref_file_idx],
                            'feature_idx': edge.get('current_feature_idx', -1)
                        }
                    
                    # Build matched_file nested structure from index
                    match_file_idx = edge.get('match_file_idx')
                    if match_file_idx is not None and match_file_idx in file_index_map:
                        edge['matched_file'] = {
                            'filename': file_index_map[match_file_idx],
                            'feature_idx': edge.get('matched_feature_idx', -1)
                        }
                    
                    # Remove index fields to clean up JSON output
                    for key in ['ref_file_idx', 'match_file_idx', 'current_feature_idx', 'matched_feature_idx']:
                        edge.pop(key, None)
                    
                    edges_converted += 1
            
            print(f"[DEBUG] Converted {edges_converted} edge records with filenames")
            
            mem_before_edges_json = get_memory_info("BEFORE EDGES JSON SAVE", print_details=True)
            with open(edges_json_path, 'w') as f:
                json.dump(edges_for_json, f, indent=2)
            mem_after_edges_json = get_memory_info("AFTER EDGES JSON SAVE", print_details=True)
            print(f"[DEBUG] Memory delta for edges JSON save: {mem_after_edges_json - mem_before_edges_json:.1f} MB")
            print(f"✓ Saved edges dictionary (JSON): {os.path.basename(edges_json_path)}")
            
            # [AGGRESSIVE CLEANUP] Delete the duplicate edges_for_json immediately after saving
            print("  Cleaning up temporary edges_for_json copy...")
            del edges_for_json
            gc.collect()
            mem_after_cleanup = get_memory_info("AFTER EDGES_FOR_JSON CLEANUP", print_details=False)
            print(f"    Memory freed: {mem_after_edges_json - mem_after_cleanup:.1f} MB")
        except Exception as e:
            print(f"⊘ Failed to save edges JSON: {e}")
    else:
        print(f"⊘ Skipped edges JSON output (--skip_json_output enabled)")
    
    # [DELETION POINT #7] Final cleanup after all output operations complete
    print("Final memory cleanup...")
    del output_data
    del network_edges
    del ident_lookup
    gc.collect()
    print("  ✓ Output data structures freed")
    
    # Save metadata JSON (pairing parameters, cutoff info, feature counts)
    # This contains all metadata from output_data EXCEPT network_edges and ident_lookup
    print("[DEBUG] Starting metadata JSON save...")
    metadata_json_path = args.output_json.replace('.json', '_metadata.json')
    metadata = {
        'total_features_before_filtering': total_features_before_filtering,
        'total_features_after_filtering': total_features_after_filter,
        'cutoff_params': {
            'mz_cutoff': args.mz_cutoff,
            'rt_cutoff': args.rt_cutoff,
            'edges_cutoff': args.edges_cutoff
        },
        'pairing_parameters': {
            'pairing_method': 'kdtree',
            'best_match_only': True,
            'distance_calc_before_scaling': args.distance_calc_before_scaling,
            'normalize_coordinates': normalize_coordinates,
            'postpair_normalize_coordinates': args.postpair_normalize_coordinates,
            'postpair_note': 'If postpair_normalize_coordinates=true: distance is original euclidean distance; distance_rescaled is from rescaled coordinate space',
            'match_cutoff': args.match_cutoff,
            'input_files': num_input_files,
            'skip_matchfinder': False
        }
    }
    
    try:
        with open(metadata_json_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"✓ Saved pairing metadata: {os.path.basename(metadata_json_path)}")
    except Exception as e:
        print(f"✗ Failed to save metadata JSON: {e}")
    
    # NOTE: Filename restoration is now done DURING JSON save (see edges JSON write section above)
    # This eliminates the memory spike from re-reading the entire JSON file
    print("[DEBUG] Skipping post-process filename restoration (already done during JSON write)")


if __name__ == "__main__":
    main()
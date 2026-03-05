#!/usr/bin/env python3
"""
Extract feature data from TSV file and optionally overlay on heatmap.
Refactored from map_mzml_batch_tsv.py
"""

import argparse
import pickle
import json
import numpy as np
import pandas as pd
import os
import ast
from typing import Dict, List


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


def parse_nested_list(s):
    """Parse nested lists from TSV - handles structure like [[1,2,3], [4,5,6]]."""
    if isinstance(s, str):
        try:
            parsed = ast.literal_eval(s)
            if parsed and isinstance(parsed[0], list):
                result = []
                for inner_list in parsed:
                    inner_result = []
                    for x in inner_list:
                        try:
                            inner_result.append(float(x))
                        except (ValueError, TypeError):
                            inner_result.append(float(str(x).strip()))
                    result.append(inner_result)
                return result
            else:
                result = []
                for x in parsed:
                    try:
                        result.append(float(x))
                    except (ValueError, TypeError):
                        result.append(float(str(x).strip()))
                return result
        except Exception as e:
            print(f"    Warning: Failed to parse '{s[:50]}...': {e}")
            return []
    return []


def parse_scalar(s):
    """Parse scalar values from TSV - handles single numeric values."""
    if isinstance(s, (int, float)):
        return int(s) if isinstance(s, int) else int(float(s))
    if isinstance(s, str):
        try:
            val = ast.literal_eval(s)
            return int(val) if isinstance(val, (int, float)) else int(float(val))
        except:
            return None
    return None


def extract_feature_data(features_tsv: str, filename_base: str, round_up_to: int = 2, rt_start_trim: float = 0, rt_end_trim: float = 0) -> List[Dict]:
    """
    Extract feature data from TSV file.
    
    Args:
        features_tsv: Path to TSV feature file
        filename_base: Base filename for output
        round_up_to: Decimal places for m/z rounding
        rt_start_trim: Seconds to trim from the start of the measurement (retention time)
        rt_end_trim: Maximum RT value to keep in seconds (keep only features with center <= this value, 0 = no limit)
    """
    print(f"  Loading features from TSV...")
    features_data = pd.read_csv(features_tsv, sep='\t')
    print(f"    Loaded {len(features_data)} rows")
    
    # Validate that this is a feature TSV, not a PSM TSV
    required_columns = {'l_mz_start', 'l_mz_end', 'l_rt_start', 'l_rt_end', 'charge'}
    actual_columns = set(features_data.columns)
    missing_columns = required_columns - actual_columns
    
    if missing_columns:
        print(f"\n✗ ERROR: Input TSV file is missing required columns:")
        print(f"  Missing: {', '.join(sorted(missing_columns))}")
        print(f"  Available columns: {', '.join(sorted(list(actual_columns)[:5]))}...")
        print(f"\n  This appears to be a PSM (identification) TSV file, not a feature TSV file.")
        print(f"  Expected input: Feature TSV from OpenMS feature extraction")
        print(f"  Expected columns: {', '.join(sorted(required_columns))}")
        
        # Check if this looks like a PSM file
        if 'psm_id' in actual_columns or 'peptide_seq' in actual_columns:
            print(f"\n  This looks like a PSM/identification result file.")
            print(f"  Please use the feature TSV files from the quantification step instead.")
        return []
    
    if len(features_data) == 0:
        print(f"  ✗ TSV file contains no data rows")
        return []
    
    print(f"    Loaded {len(features_data)} features")
    if rt_start_trim > 0 or rt_end_trim > 0:
        start_msg = f"start +{rt_start_trim:.1f} sec" if rt_start_trim > 0 else ""
        end_msg = f"end <= {rt_end_trim:.1f} sec" if rt_end_trim > 0 else ""
        trim_msg = ", ".join(filter(None, [start_msg, end_msg]))
        print(f"    RT trimming: {trim_msg}")
    
    feature_data_list = []
    box_count = 0
    skipped_count = 0
    
    print(f"  Processing {len(features_data)} features...")
    for idx, row in features_data.iterrows():
        if idx % 100 == 0:
            print(f"    Processing feature {idx}/{len(features_data)}", end='\r')
        
        # Parse list-based coordinates and intensities from TSV
        mz_starts = parse_list(row['l_mz_start'])
        mz_ends = parse_list(row['l_mz_end'])
        rt_starts = parse_list(row['l_rt_start'])
        rt_ends = parse_list(row['l_rt_end'])
        charge = parse_scalar(row['charge'])
        
        # Parse nested lists of actual coordinate data
        l_mass_to_charges_raw = parse_nested_list(row['l_mass_to_charges'])
        l_intensities_raw = parse_nested_list(row['l_intensities'])
        l_retention_times_raw = parse_nested_list(row['l_retention_times'])
        
        pep_ident = parse_list_flexible(row['l_pep_ident'])
        
        # Build feature_ls: list of (m/z, retention_time, intensity) tuples
        feature_ls = []
        
        if l_mass_to_charges_raw and isinstance(l_mass_to_charges_raw[0], list):
            for mz_list, rt_list, int_list in zip(l_mass_to_charges_raw, l_retention_times_raw, l_intensities_raw):
                for mz, rt, intensity in zip(mz_list, rt_list, int_list):
                    feature_ls.append((float(mz), float(rt*60), float(intensity)))
        elif l_mass_to_charges_raw:
            for mz, rt, intensity in zip(l_mass_to_charges_raw, l_retention_times_raw, l_intensities_raw):
                feature_ls.append((float(mz), float(rt*60), float(intensity)))
        
        # Ensure we have valid boundaries
        if mz_starts and mz_ends and rt_starts and rt_ends and feature_ls:
            mz_start = min(mz_starts + mz_ends)
            mz_end = max(mz_starts + mz_ends)
            rt_start = min(rt_starts + rt_ends)  # In seconds
            rt_end = max(rt_starts + rt_ends)    # In seconds
            
            # Calculate geometric center from feature_ls coordinates
            if feature_ls:
                mz_coords = [f[0] for f in feature_ls]
                rt_coords = [f[1] for f in feature_ls]
                x_center_geo = (min(mz_coords) + max(mz_coords)) / 2
                y_center_geo = (min(rt_coords) + max(rt_coords)) / 2
                x_center = x_center_geo
                y_center = y_center_geo
                #TODO: try Linear and Quadratic Discriminant Analysis with covariance ellipsoid for centroid calculation as alternative
                #TODO: try maximum from first isotop, try overall max intensity isotope as centroid
                #TODO: Try nereast real point from virtual centroid as centroid -> later for validation
                # Calculate intensity-weighted center of mass
                intensity_sum = sum(f[2] for f in feature_ls if f[2] > 0)
                
                if intensity_sum > 0:
                    for mz, rt, intensity in feature_ls:
                        if intensity > 0:
                            x_vec = (mz - x_center_geo) * (intensity / intensity_sum)
                            y_vec = (rt - y_center_geo) * (intensity / intensity_sum)
                            x_center += x_vec
                            y_center += y_vec
                
                # Apply RT trimming filter based on y_center (centered RT position in seconds)
                if rt_start_trim > 0 or rt_end_trim > 0:
                    rt_start_min = rt_start_trim if rt_start_trim > 0 else float('-inf')  # Absolute lower bound
                    rt_end_max = rt_end_trim if rt_end_trim > 0 else float('inf')  # Absolute upper bound
                    
                    # Skip feature if its center is outside the allowed RT range
                    if not (rt_start_min <= y_center <= rt_end_max):
                        skipped_count += 1
                        continue
                
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
                    'charge': charge,
                    'filename': filename_base,
                    'feature_count': len(feature_ls)
                })
                box_count += 1
        else:
            skipped_count += 1
    
    print(f"  Extracted {box_count} features ({skipped_count} skipped)")
    return feature_data_list


def main():
    parser = argparse.ArgumentParser(
        description="Extract feature data from TSV file"
    )
    parser.add_argument("--tsv", required=True, help="Path to TSV feature file")
    parser.add_argument("--mzml", help="Path to mzML file (optional, for future use)")
    parser.add_argument("--spectrum_pkl", help="Path to spectrum pickle file (optional, for future use)")
    parser.add_argument("--output_pkl", required=True, help="Output pickle file path")
    parser.add_argument("--output_json", required=True, help="Output JSON file path")
    parser.add_argument("--round_up_to", type=int, default=2, help="Decimal places for m/z rounding")
    parser.add_argument("--feature_mode", default="CoM", help="Feature center calculation mode")
    parser.add_argument("--generate_diagnostic", type=bool, default=False, help="Generate diagnostic plots")
    parser.add_argument("--rt_start_trim", type=float, default=0, help="Seconds to trim from the start of retention time")
    parser.add_argument("--rt_end_trim", type=float, default=0, help="Maximum RT value to keep (keep only features with RT <= this value, 0 = no limit)")
    
    args = parser.parse_args()
    
    basename = os.path.splitext(os.path.basename(args.tsv))[0]
    
    print(f"\n{'='*70}")
    print(f"Extracting Features: {basename}")
    print(f"{'='*70}")
    
    # Extract feature data
    feature_data_list = extract_feature_data(
        args.tsv,
        basename,
        round_up_to=args.round_up_to,
        rt_start_trim=args.rt_start_trim,
        rt_end_trim=args.rt_end_trim
    )
    
    # Save as pickle
    with open(args.output_pkl, 'wb') as f:
        pickle.dump(feature_data_list, f)
    print(f"✓ Saved pickle: {os.path.basename(args.output_pkl)}")
    
    # Save as JSON
    with open(args.output_json, 'w') as f:
        json.dump(feature_data_list, f, indent=2)
    print(f"✓ Saved JSON: {os.path.basename(args.output_json)}")


if __name__ == "__main__":
    main()

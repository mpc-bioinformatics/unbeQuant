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


def extract_feature_data(features_tsv: str, filename_base: str, round_up_to: int = 2) -> List[Dict]:
    """
    Extract feature data from TSV file.
    """
    print(f"  Loading features from TSV...")
    features_data = pd.read_csv(features_tsv, sep='\t')
    print(f"    Loaded {len(features_data)} features")
    
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
                feature_ls.append((float(mz), float((rt*60)), float(intensity)))
        
        # Ensure we have valid boundaries
        if mz_starts and mz_ends and rt_starts and rt_ends and feature_ls:
            mz_start = min(mz_starts + mz_ends)
            mz_end = max(mz_starts + mz_ends)
            rt_start = min(rt_starts + rt_ends)
            rt_end = max(rt_starts + rt_ends)
            
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
    
    args = parser.parse_args()
    
    basename = os.path.splitext(os.path.basename(args.tsv))[0]
    
    print(f"\n{'='*70}")
    print(f"Extracting Features: {basename}")
    print(f"{'='*70}")
    
    # Extract feature data
    feature_data_list = extract_feature_data(
        args.tsv,
        basename,
        round_up_to=args.round_up_to
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

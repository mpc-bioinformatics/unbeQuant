#!/usr/bin/env python3
"""
Process mzML file and extract spectral data.
Refactored from map_mzml_batch_tsv.py
"""

import argparse
import pickle
import numpy as np
import pyopenms as oms
import os
from typing import Tuple, List, Dict


def process_mzml_file(mzml_path: str, round_up_to: int = 2) -> Tuple[np.ndarray, List[float], List[float], Dict, Dict]:
    """
    Process a single mzML file.
    Returns: (spec_total, RTINSECONDS_arr, mz_total_arr, mz_dict, rt_dict)
    """
    print(f"  Loading mzML: {os.path.basename(mzml_path)}")
    
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


def main():
    parser = argparse.ArgumentParser(
        description="Process mzML file and extract spectral data"
    )
    parser.add_argument("--mzml", required=True, help="Path to mzML file")
    parser.add_argument("--output_pickle", required=True, help="Output pickle file path")
    parser.add_argument("--round_up_to", type=int, default=2, help="Decimal places for m/z rounding")
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print(f"Processing mzML: {os.path.basename(args.mzml)}")
    print(f"{'='*70}")
    
    spec_total, RTINSECONDS_arr, mz_total_arr, mz_dict, rt_dict = process_mzml_file(
        args.mzml, 
        round_up_to=args.round_up_to
    )
    
    # Save to pickle
    data = {
        'spec_total': spec_total,
        'RTINSECONDS_arr': RTINSECONDS_arr,
        'mz_total_arr': mz_total_arr,
        'mz_dict': mz_dict,
        'rt_dict': rt_dict
    }
    
    with open(args.output_pickle, 'wb') as f:
        pickle.dump(data, f)
    
    print(f"\n✓ Saved spectrum data to {os.path.basename(args.output_pickle)}")


if __name__ == "__main__":
    main()

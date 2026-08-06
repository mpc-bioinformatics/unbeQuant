#!/usr/bin/env python3
"""
Extract feature data from a per-run feature TSV.

Reads an OpenMS-derived feature TSV and emits one dict per feature (pickle + JSON) with
the coordinates the downstream clusterless workflow needs:

  - mz_start/mz_end/rt_start/rt_end : the raw bounding box, taken straight from the
                                      l_* boundary columns (NOT affected by feature_mode)
  - x_center_geo/y_center_geo       : geometric midrange of the selected peak set
  - x_center/y_center               : the feature "point" used for pairing — the geometric
                                      midrange plus (for weighted modes) an intensity
                                      centre-of-mass shift
  - identifications, charge, intensity, openms_fid, ms2_scans

`--feature_mode` selects WHICH peaks define the centre (see calculate_feature_center) and
whether the intensity weighting is applied. RT values coming from the nested peak lists are
in minutes and are converted to SECONDS here; the l_rt_* boundary columns are already in
seconds.

This script fails loudly: a missing column, an unparseable cell, or an unknown feature_mode
aborts the run with a message naming the offending row/column, rather than silently
emitting a degraded feature that only breaks further down the pipeline.

Consumed by: pair_features_1_clusterless.py (via the .pkl) and
compare_retention_time_alignment.py (via the .json).
"""

import argparse
import ast
import json
import math
import os
import pickle
import sys
from typing import Dict, List

import pandas as pd


# ============================================================================
# TSV field parsing
# ============================================================================

def _abort(message: str, column: str = None, row_idx=None, value=None) -> None:
    """Print a precise parse/validation failure and stop the run."""
    print(f"\n✗ ERROR: {message}")
    if column is not None:
        print(f"    column : {column}")
    if row_idx is not None:
        print(f"    row    : {row_idx}")
    if value is not None:
        print(f"    value  : {value!r}")
    sys.exit(1)


def parse_list(value, column: str, row_idx) -> List[float]:
    """Parse a flat list literal (e.g. "[1.0, 2.0]") into floats. Aborts if unparseable."""
    try:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
        return [float(x) for x in parsed]
    except Exception as e:
        _abort(f"could not parse a numeric list ({e})", column, row_idx, value)


def parse_list_flexible(value, column: str, row_idx) -> List:
    """
    Parse a flat list literal, keeping non-numeric entries as strings.

    Used for identification columns, whose entries are peptide/protein accessions.
    An empty list ("[]") is valid — most features carry no identification.
    """
    try:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
        iter(parsed)
    except Exception as e:
        _abort(f"could not parse an identification list ({e})", column, row_idx, value)

    result = []
    for x in parsed:
        try:
            result.append(float(x))
        except (ValueError, TypeError):
            result.append(str(x))
    return result


def parse_nested_list(value, column: str, row_idx) -> List:
    """
    Parse a nested list literal like "[[1,2],[3,4]]" (one inner list per isotope) into floats.

    Flat input (no isotope nesting) is returned as a flat float list. Aborts if unparseable.
    """
    def to_float(x):
        try:
            return float(x)
        except (ValueError, TypeError):
            return float(str(x).strip())

    try:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
        if parsed and isinstance(parsed[0], list):
            return [[to_float(x) for x in inner] for inner in parsed]
        return [to_float(x) for x in parsed]
    except Exception as e:
        _abort(f"could not parse a nested numeric list ({e})", column, row_idx, value)


def parse_scalar(value, column: str, row_idx) -> int:
    """Parse a single integer value. Aborts if the cell is empty/NaN or not numeric."""
    try:
        if isinstance(value, str):
            value = ast.literal_eval(value)
        if isinstance(value, float) and math.isnan(value):
            raise ValueError("cell is empty (NaN)")
        return int(value)
    except Exception as e:
        _abort(f"could not parse an integer ({e})", column, row_idx, value)


# ============================================================================
# Feature centre calculation
# ============================================================================

def calculate_feature_center(l_mass_to_charges_raw: List, l_intensities_raw: List,
                             l_retention_times_raw: List, feature_mode: str):
    """
    Reduce one feature's raw per-isotope peak lists to a single centre point.

    Everything the feature_mode controls lives in the one dispatch block below: which peaks
    the mode uses, and whether the centre is intensity-weighted. RT is converted from
    minutes to SECONDS here (rt * 60).

        mode                          peaks used                          weighted?
        ----------------------------  ----------------------------------  ---------
        "CoM"                         every peak of every isotope         yes
        "first_mean"                  every peak of the first isotope     yes
        "first_mean_unweighted_mean"  every peak of the first isotope     no
        "first_mean_unweighted_first" earliest-eluting peak, 1st isotope  no
        "single_peak"                 most intense peak, all isotopes     yes (no-op: 1 peak)

    "weighted" means the centre is the intensity-weighted mean of the peaks; unweighted
    leaves it at the plain geometric midrange (so x/y_center == x/y_center_geo).

    Returns (x_center, y_center, x_center_geo, y_center_geo, n_peaks), or None if the
    feature has no usable peaks. Raises ValueError on an unknown feature_mode.
    """
    if not l_mass_to_charges_raw:
        return None

    # Flat input (no isotope nesting) is treated as a single isotope, so the code below
    # only ever deals with the nested [isotope][peak] shape.
    if not isinstance(l_mass_to_charges_raw[0], list):
        l_mass_to_charges_raw = [l_mass_to_charges_raw]
        l_retention_times_raw = [l_retention_times_raw]
        l_intensities_raw = [l_intensities_raw]

    def peaks_of(isotope_count):
        """(m/z, RT_seconds, intensity) of the first N isotopes; None = all isotopes."""
        return [
            (float(mz), float(rt * 60), float(intensity))
            for mz_list, rt_list, int_list in zip(l_mass_to_charges_raw[:isotope_count],
                                                  l_retention_times_raw[:isotope_count],
                                                  l_intensities_raw[:isotope_count])
            for mz, rt, intensity in zip(mz_list, rt_list, int_list)
        ]

    # ---- What each mode does, all in one place -----------------------------
    if feature_mode == "CoM":
        peaks, weighted = peaks_of(None), True

    elif feature_mode == "first_mean":
        peaks, weighted = peaks_of(1), True

    elif feature_mode == "first_mean_unweighted_mean":
        peaks, weighted = peaks_of(1), False

    elif feature_mode == "first_mean_unweighted_first":
        first_isotope = peaks_of(1)
        peaks = [min(first_isotope, key=lambda p: p[1])] if first_isotope else []
        weighted = False

    elif feature_mode == "single_peak":
        # Strongest peak anywhere in the feature. The > -1 floor is historical: a feature
        # whose peaks are all at or below -1 intensity yields no centre at all.
        strongest = max(peaks_of(None), key=lambda p: p[2], default=None)
        peaks = [strongest] if strongest is not None and strongest[2] > -1 else []
        weighted = True

    else:
        raise ValueError(f"Unknown feature_mode: {feature_mode!r}")

    if not peaks:
        return None

    # ---- Centre of the selected peaks --------------------------------------
    mz_coords = [p[0] for p in peaks]
    rt_coords = [p[1] for p in peaks]
    x_center_geo = (min(mz_coords) + max(mz_coords)) / 2
    y_center_geo = (min(rt_coords) + max(rt_coords)) / 2
    x_center, y_center = x_center_geo, y_center_geo

    # Intensity-weighted centre of mass. Algebraically this equals
    # sum(coord * I) / sum(I) over the positive-intensity peaks: the weights sum to 1, so
    # the x/y_center_geo seed cancels out exactly.
    # NOTE: the loop variable is deliberately peak_intensity, not intensity — see the
    # feature-level `intensity` in extract_feature_data(), which this used to shadow.
    intensity_sum = sum(p[2] for p in peaks if p[2] > 0)
    if weighted and intensity_sum > 0:
        for mz, rt, peak_intensity in peaks:
            if peak_intensity > 0:
                x_center += (mz - x_center_geo) * (peak_intensity / intensity_sum)
                y_center += (rt - y_center_geo) * (peak_intensity / intensity_sum)

    return x_center, y_center, x_center_geo, y_center_geo, len(peaks)


# ============================================================================
# Extraction
# ============================================================================

# Every column this script reads. All are mandatory: silently defaulting a missing column
# would emit features with a null intensity/charge or empty identifications, which only
# fail (confusingly) much further downstream.
REQUIRED_COLUMNS = [
    'l_mz_start', 'l_mz_end', 'l_rt_start', 'l_rt_end',
    'charge', 'intensity',
    'l_mass_to_charges', 'l_intensities', 'l_retention_times',
    'l_pep_ident', 'l_raw_pep_ident', 'l_prot_ident',
]

# Pairs where either spelling is accepted (they differ between OpenMS export versions).
# The first name listed wins if both are present.
REQUIRED_EITHER = [
    ('id', 'openms_fid'),
    ('ms2_scans', 'l_ms2_scans'),
]


def _validate_feature_tsv(features_data: pd.DataFrame, features_tsv: str) -> None:
    """
    Abort with a precise message unless every column this script reads is present.

    Also detects the common mistake of passing a PSM/identification TSV instead of a
    feature TSV.
    """
    actual_columns = set(features_data.columns)
    missing = [c for c in REQUIRED_COLUMNS if c not in actual_columns]
    missing_either = [pair for pair in REQUIRED_EITHER if not set(pair) & actual_columns]

    if missing or missing_either:
        print(f"\n✗ ERROR: Input TSV is missing required column(s): {features_tsv}")
        for col in missing:
            print(f"    - {col}")
        for pair in missing_either:
            print(f"    - {' or '.join(pair)}")
        print(f"\n  Columns found ({len(actual_columns)}): {', '.join(sorted(actual_columns))}")
        if 'psm_id' in actual_columns or 'peptide_seq' in actual_columns:
            print("\n  This looks like a PSM/identification result file, not a feature TSV.")
            print("  Use the feature TSV files from the quantification step instead.")
        else:
            print("\n  Expected input: feature TSV from OpenMS feature extraction.")
        sys.exit(1)

    if len(features_data) == 0:
        print(f"\n✗ ERROR: TSV file contains no data rows: {features_tsv}")
        sys.exit(1)


def extract_feature_data(features_tsv: str, filename_base: str, feature_mode: str,
                         rt_start_trim: float = 0, rt_end_trim: float = 0) -> List[Dict]:
    """
    Extract the per-feature dicts from a feature TSV.

    Args:
        features_tsv: Path to the TSV feature file.
        filename_base: Base filename recorded on every feature (used as the file identity
            downstream).
        feature_mode: Which peaks define the centre and whether intensity weighting applies
            (see calculate_feature_center).
        rt_start_trim: Drop features whose y_center is below this RT (seconds); 0 = no bound.
        rt_end_trim: Drop features whose y_center is above this RT (seconds); 0 = no bound.

    Returns:
        List of feature dicts. Any cell that cannot be parsed aborts the run; features are
        only dropped for the structural reasons counted in the summary (no bounding box,
        no usable peaks, or trimmed out by the RT bounds).
    """
    print("  Loading features from TSV...")
    features_data = pd.read_csv(features_tsv, sep='\t')
    _validate_feature_tsv(features_data, features_tsv)
    print(f"    Loaded {len(features_data)} features")

    if rt_start_trim > 0 or rt_end_trim > 0:
        start_msg = f"start +{rt_start_trim:.1f} sec" if rt_start_trim > 0 else ""
        end_msg = f"end <= {rt_end_trim:.1f} sec" if rt_end_trim > 0 else ""
        print(f"    RT trimming: {', '.join(filter(None, [start_msg, end_msg]))}")

    # Resolve the either/or column spellings once (validation guarantees one of each).
    columns = set(features_data.columns)
    fid_col = 'id' if 'id' in columns else 'openms_fid'
    ms2_col = 'ms2_scans' if 'ms2_scans' in columns else 'l_ms2_scans'

    feature_data_list = []
    skipped_no_bounds = 0   # empty m/z or RT boundary list -> no bounding box
    skipped_no_peaks = 0    # feature_mode selected no usable peaks
    skipped_rt_trim = 0     # centre outside the requested RT window

    print(f"  Processing {len(features_data)} features...")
    for idx, row in features_data.iterrows():
        if idx % 100 == 0:
            print(f"    Processing feature {idx}/{len(features_data)}", end='\r')

        # Bounding box + charge (the l_rt_* columns are already in seconds).
        mz_starts = parse_list(row['l_mz_start'], 'l_mz_start', idx)
        mz_ends = parse_list(row['l_mz_end'], 'l_mz_end', idx)
        rt_starts = parse_list(row['l_rt_start'], 'l_rt_start', idx)
        rt_ends = parse_list(row['l_rt_end'], 'l_rt_end', idx)
        charge = parse_scalar(row['charge'], 'charge', idx)

        # Per-isotope peak lists (RT here is in minutes).
        l_mass_to_charges_raw = parse_nested_list(row['l_mass_to_charges'], 'l_mass_to_charges', idx)
        l_intensities_raw = parse_nested_list(row['l_intensities'], 'l_intensities', idx)
        l_retention_times_raw = parse_nested_list(row['l_retention_times'], 'l_retention_times', idx)

        pep_ident = parse_list_flexible(row['l_pep_ident'], 'l_pep_ident', idx)
        raw_pep_ident = parse_list_flexible(row['l_raw_pep_ident'], 'l_raw_pep_ident', idx)
        prot_ident = parse_list_flexible(row['l_prot_ident'], 'l_prot_ident', idx)

        intensity = float(row['intensity'])
        openms_fid = row[fid_col]
        ms2_scans = row[ms2_col]

        if not (mz_starts and mz_ends and rt_starts and rt_ends):
            skipped_no_bounds += 1
            continue

        # The peaks that define this feature's centre.
        center = calculate_feature_center(l_mass_to_charges_raw, l_intensities_raw,
                                          l_retention_times_raw, feature_mode)
        if not center:
            skipped_no_peaks += 1
            continue

        x_center, y_center, x_center_geo, y_center_geo, n_peaks = center

        mz_start = min(mz_starts + mz_ends)
        mz_end = max(mz_starts + mz_ends)
        rt_start = min(rt_starts + rt_ends)  # seconds
        rt_end = max(rt_starts + rt_ends)    # seconds

        # RT trimming, applied to the feature centre.
        if rt_start_trim > 0 or rt_end_trim > 0:
            rt_min = rt_start_trim if rt_start_trim > 0 else float('-inf')
            rt_max = rt_end_trim if rt_end_trim > 0 else float('inf')
            if not (rt_min <= y_center <= rt_max):
                skipped_rt_trim += 1
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
            'l_raw_pep_ident': raw_pep_ident,
            'prot_ident': prot_ident,
            'intensity': intensity,
            'openms_fid': openms_fid,
            'ms2_scans': ms2_scans,
            'charge': charge,
            'filename': filename_base,
            'feature_count': n_peaks,
        })

    total_skipped = skipped_no_bounds + skipped_no_peaks + skipped_rt_trim
    print(f"  Extracted {len(feature_data_list)} features ({total_skipped} skipped)")
    if total_skipped:
        print(f"    skipped - no bounding box  : {skipped_no_bounds}")
        print(f"    skipped - no usable peaks  : {skipped_no_peaks}")
        print(f"    skipped - outside RT window: {skipped_rt_trim}")
    return feature_data_list


# ============================================================================
# main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Extract feature data from TSV file")
    parser.add_argument("--tsv", required=True, help="Path to TSV feature file")
    parser.add_argument("--output_pkl", required=True, help="Output pickle file path")
    parser.add_argument("--output_json", required=True, help="Output JSON file path")
    parser.add_argument("--feature_mode", required=True,
                        choices=["CoM", "first_mean", "single_peak",
                                 "first_mean_unweighted_mean", "first_mean_unweighted_first"],
                        help="Which peaks define the feature centre: 'CoM' (all isotopes), "
                             "'first_mean' (first isotope), 'single_peak' (strongest peak), "
                             "'first_mean_unweighted_mean' (first isotope, geometric midrange, "
                             "no intensity weighting), 'first_mean_unweighted_first' "
                             "(first isotope, earliest-eluting peak)")
    parser.add_argument("--rt_start_trim", type=float, default=0,
                        help="Drop features whose centre RT is below this value (seconds, 0 = no bound)")
    parser.add_argument("--rt_end_trim", type=float, default=0,
                        help="Drop features whose centre RT is above this value (seconds, 0 = no bound)")

    args = parser.parse_args()

    basename = os.path.splitext(os.path.basename(args.tsv))[0]

    print(f"\n{'=' * 70}")
    print(f"Extracting Features: {basename}")
    print(f"{'=' * 70}")

    feature_data_list = extract_feature_data(
        args.tsv,
        basename,
        feature_mode=args.feature_mode,
        rt_start_trim=args.rt_start_trim,
        rt_end_trim=args.rt_end_trim
    )

    with open(args.output_pkl, 'wb') as f:
        pickle.dump(feature_data_list, f)
    print(f"✓ Saved pickle: {os.path.basename(args.output_pkl)}")

    with open(args.output_json, 'w') as f:
        json.dump(feature_data_list, f, indent=2)
    print(f"✓ Saved JSON: {os.path.basename(args.output_json)}")


if __name__ == "__main__":
    main()

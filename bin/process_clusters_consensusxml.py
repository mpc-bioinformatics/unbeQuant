#!/usr/bin/env python
"""
Process filtered/clusters ConsensusXML file and generate TSV reports with comparison to original.

This script processes a ConsensusXML file (e.g., filtered clusters) in the same way as the 
main workflow processes consensus files. It generates TSV reports and compares the filtered
consensus with the original consensus to show what was retained/filtered.

Usage (standalone):
    python process_clusters_consensusxml.py \
        -consensus results/feature_analysis/feature_data_lists/clusters_fraction1.0_filtered.consensusXML \
        -feature_tsvs_dir results/feature_analysis/features_with_annotated_identifications \
        -out_dir results/feature_analysis/clusters_analysis

Usage (nextflow):
    Provides the same parameters via nextflow process

"""

import argparse
import os
import sys
import json
import csv
import glob
from pathlib import Path
from ast import literal_eval
from datetime import datetime
from collections import defaultdict

import tqdm
import pandas as pd
import numpy as np
import pyopenms

# For visualization
try:
    import plotly.graph_objects as go
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

csv.field_size_limit(sys.maxsize)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Process filtered/clusters ConsensusXML and generate reports'
    )
    parser.add_argument('-consensus', required=True,
                       help='Input ConsensusXML file (e.g., filtered clusters)')
    parser.add_argument('-original_consensus', default='',
                       help='Original ConsensusXML file for comparison (optional)')
    parser.add_argument('-feature_tsvs_dir', required=True,
                       help='Directory containing feature TSV files')
    parser.add_argument('-out_dir', default='results/feature_analysis/clusters_analysis',
                       help='Output directory for results')
    parser.add_argument('-verbose', action='store_true',
                       help='Enable verbose output')
    
    return parser.parse_args()


def load_feature_tsvs(feature_dir, verbose=False):
    """Load all feature TSV files from directory."""
    dict_of_single_features = {}
    filenames = []
    
    pattern = os.path.join(feature_dir, '*__.tsv')  # Match feature TSV files
    tsv_files = glob.glob(pattern)
    
    if not tsv_files:
        # Try alternative pattern
        pattern = os.path.join(feature_dir, '*.tsv')
        tsv_files = [f for f in glob.glob(pattern) if not 'cutoff' in f and not 'plot' in f]
    
    if verbose:
        print(f"Found {len(tsv_files)} feature TSV files\n")
    
    for ftsv in sorted(tsv_files):
        filename_key = os.path.basename(ftsv).split('_____')[0]
        
        try:
            # Try loading with literal_eval for list columns
            df = pd.read_csv(ftsv, sep='\t', 
                            converters={
                                col: literal_eval for col in [
                                    'l_pep_ident', 'l_raw_pep_ident', 'l_prot_ident',
                                    'l_mz_start', 'l_mz_end', 'l_rt_start', 'l_rt_end',
                                    'l_retention_times', 'l_mass_to_charges', 'l_intensities',
                                    'l_ms2_scans'
                                ] if col in pd.read_csv(ftsv, sep='\t', nrows=1).columns
                            })
        except Exception as e:
            if verbose:
                print(f"  Warning loading {filename_key}: {e}, trying without converters")
            df = pd.read_csv(ftsv, sep='\t')
        
        dict_of_single_features[filename_key] = df
        filenames.append(filename_key)
        
        if verbose:
            print(f"  ✓ {filename_key}: {len(df)} features")
    
    return dict_of_single_features, sorted(filenames)


def process_consensusxml(consensus_file, feature_dict, filenames, verbose=False):
    """Process ConsensusXML file and extract feature information."""
    
    if verbose:
        print(f"\nLoading ConsensusXML: {consensus_file}")
    
    cmap = pyopenms.ConsensusMap()
    pyopenms.ConsensusXMLFile().load(consensus_file, cmap)
    
    # Get column headers
    map_list = cmap.getColumnHeaders()
    for key, val in map_list.items():
        map_list[key] = ".".join(val.filename.split(".")[:-1])
    
    if verbose:
        print(f"Loaded {len([None for _ in cmap])} consensus features from {len(map_list)} files")
    
    # Define output headers (same structure as original workflow)
    header_per_file = [
        "intensity",
        "l_pep_ident",
        "l_raw_pep_ident",
        "l_prot_ident",
        "openms_fid",
        "charge",
        "l_ms2_scans",
        "l_mz_start",
        "l_mz_end",
        "l_rt_start",
        "l_rt_end",
        "l_retention_times",
        "l_mass_to_charges",
        "l_intensities"
    ]
    
    header_per_file_reduced = [
        "intensity",
        "l_pep_ident",
        "l_raw_pep_ident",
        "l_prot_ident",
        "openms_fid",
        "charge",
        "l_ms2_scans",
    ]
    
    header_per_minimal = [
        "intensity",
        "openms_fid",
        "charge",
        "l_ms2_scans",
    ]
    
    global_headers = [
        "first_iso_global_min_mz",
        "first_iso_global_max_mz",
        "first_iso_global_min_rt",
        "first_iso_global_max_rt",
        "feature_global_min_mz",
        "feature_global_max_mz",
        "feature_global_min_rt",
        "feature_global_max_rt"
    ]
    
    # Extract features
    features_list = []
    sanity_check = {f: [False]*len(feature_dict[f]) for f in filenames}
    
    if verbose:
        print(f"Extracting consensus features...")
    
    for ce in tqdm.tqdm(cmap, total=len([None for _ in cmap]), unit="features", disable=not verbose):
        entry = ["e_" + str(ce.getUniqueId())]
        entry_reduced = ["e_" + str(ce.getUniqueId())]
        entry_minimal = ["e_" + str(ce.getUniqueId())]
        
        feature_list = ce.getFeatureList()
        files_involved = [map_list[fl.getMapIndex()] for fl in feature_list]
        feature_ids = [fl.getUniqueId() for fl in feature_list]
        
        # Preload rows
        file_rows = []
        for fi, fi_fid in zip(files_involved, feature_ids):
            if fi not in feature_dict:
                if verbose:
                    print(f"    Warning: File {fi} not found in feature dict, skipping")
                continue
            
            filter_row = feature_dict[fi]["openms_fid"] == "f_" + str(fi_fid)
            if not filter_row.any():
                if verbose:
                    print(f"    Warning: Feature f_{fi_fid} not found in {fi}")
                continue
            
            file_rows.append(feature_dict[fi][filter_row])
            sanity_check[fi][file_rows[-1].index[0]] = True
        
        if not file_rows:
            continue
        
        # Calculate global stats
        mzs, mze, rts, rte = [], [], [], []
        for f in filenames:
            if f in files_involved:
                fidx = files_involved.index(f)
                if fidx < len(file_rows):
                    try:
                        mz_start = file_rows[fidx]["l_mz_start"].values[0]
                        mz_end = file_rows[fidx]["l_mz_end"].values[0]
                        rt_start = file_rows[fidx]["l_rt_start"].values[0]
                        rt_end = file_rows[fidx]["l_rt_end"].values[0]
                        
                        if isinstance(mz_start, list) and len(mz_start) > 0:
                            mzs.append(float(mz_start[0]))
                        if isinstance(mz_end, list) and len(mz_end) > 0:
                            mze.append(float(mz_end[0]))
                        if isinstance(rt_start, list) and len(rt_start) > 0:
                            rts.append(float(rt_start[0]))
                        if isinstance(rt_end, list) and len(rt_end) > 0:
                            rte.append(float(rt_end[0]))
                    except:
                        pass
        
        g_mzs, g_mze, g_rts, g_rte = [], [], [], []
        for f in filenames:
            if f in files_involved:
                fidx = files_involved.index(f)
                if fidx < len(file_rows):
                    try:
                        for val in file_rows[fidx]["l_mz_start"].values[0]:
                            g_mzs.append(float(val))
                        for val in file_rows[fidx]["l_mz_end"].values[0]:
                            g_mze.append(float(val))
                        for val in file_rows[fidx]["l_rt_start"].values[0]:
                            g_rts.append(float(val))
                        for val in file_rows[fidx]["l_rt_end"].values[0]:
                            g_rte.append(float(val))
                    except:
                        pass
        
        global_info = [
            min(mzs) if mzs else 0,
            max(mze) if mze else 0,
            min(rts) if rts else 0,
            max(rte) if rte else 0,
            min(g_mzs) if g_mzs else 0,
            max(g_mze) if g_mze else 0,
            min(g_rts) if g_rts else 0,
            max(g_rte) if g_rte else 0
        ]
        
        # Build full entry
        entry_data = {"full": entry, "reduced": entry_reduced, "minimal": entry_minimal}
        
        for h in header_per_file:
            for f in filenames:
                if f in files_involved:
                    fidx = files_involved.index(f)
                    if fidx < len(file_rows):
                        entry_data["full"].append(file_rows[fidx][h].values[0])
                    else:
                        entry_data["full"].append(None)
                else:
                    entry_data["full"].append(None)
        
        for h in header_per_file_reduced:
            for f in filenames:
                if f in files_involved:
                    fidx = files_involved.index(f)
                    if fidx < len(file_rows):
                        entry_data["reduced"].append(file_rows[fidx][h].values[0])
                    else:
                        entry_data["reduced"].append(None)
                else:
                    entry_data["reduced"].append(None)
        
        for h in header_per_minimal:
            for f in filenames:
                if f in files_involved:
                    fidx = files_involved.index(f)
                    if fidx < len(file_rows):
                        entry_data["minimal"].append(file_rows[fidx][h].values[0])
                    else:
                        entry_data["minimal"].append(None)
                else:
                    entry_data["minimal"].append(None)
        
        entry_data["full"].extend(global_info)
        entry_data["reduced"].extend(global_info)
        entry_data["minimal"].extend(global_info)
        
        features_list.append(entry_data)
    
    return features_list, header_per_file, header_per_file_reduced, header_per_minimal, global_headers, filenames


def write_tsv_outputs(features_list, header_per_file, header_per_file_reduced, header_per_minimal,
                      global_headers, filenames, out_dir, verbose=False):
    """Write TSV output files."""
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Create headers
    row_header_full = ["openms_ceid"]
    for h in header_per_file:
        for f in filenames:
            row_header_full.append(f + "_____" + h)
    row_header_full.extend(global_headers)
    
    row_header_reduced = ["openms_ceid"]
    for h in header_per_file_reduced:
        for f in filenames:
            row_header_reduced.append(f + "_____" + h)
    row_header_reduced.extend(global_headers)
    
    row_header_minimal = ["openms_ceid"]
    for h in header_per_minimal:
        for f in filenames:
            row_header_minimal.append(f + "_____" + h)
    row_header_minimal.extend(global_headers)
    
    # Write files
    output_files = {}
    
    with open(os.path.join(out_dir, 'clusters_full.tsv'), 'w', newline='') as f_full, \
         open(os.path.join(out_dir, 'clusters_reduced.tsv'), 'w', newline='') as f_reduced, \
         open(os.path.join(out_dir, 'clusters_minimal.tsv'), 'w', newline='') as f_minimal:
        
        writer_full = csv.writer(f_full, delimiter='\t')
        writer_reduced = csv.writer(f_reduced, delimiter='\t')
        writer_minimal = csv.writer(f_minimal, delimiter='\t')
        
        writer_full.writerow(row_header_full)
        writer_reduced.writerow(row_header_reduced)
        writer_minimal.writerow(row_header_minimal)
        
        for feature_data in features_list:
            writer_full.writerow(feature_data["full"])
            writer_reduced.writerow(feature_data["reduced"])
            writer_minimal.writerow(feature_data["minimal"])
    
    output_files['full'] = os.path.join(out_dir, 'clusters_full.tsv')
    output_files['reduced'] = os.path.join(out_dir, 'clusters_reduced.tsv')
    output_files['minimal'] = os.path.join(out_dir, 'clusters_minimal.tsv')
    
    if verbose:
        print(f"\n✓ Generated TSV files:")
        for key, path in output_files.items():
            print(f"  - {path}")
    
    return output_files


def generate_xlsx_reports(tsv_files, out_dir, verbose=False):
    """Generate XLSX reports from TSV files using write_xlsx_report.py."""
    
    import subprocess
    
    xlsx_files = []
    
    for key, tsv_file in tsv_files.items():
        xlsx_name = f"clusters_{key}.report.xlsx"
        xlsx_path = os.path.join(out_dir, xlsx_name)
        
        try:
            cmd = [
                'python', 
                os.path.join(os.path.dirname(__file__), 'write_xlsx_report.py'),
                '-i', tsv_file,
                '-o', xlsx_path
            ]
            subprocess.run(cmd, check=True, capture_output=True)
            xlsx_files.append(xlsx_path)
            
            if verbose:
                print(f"  ✓ Generated {xlsx_name}")
        except Exception as e:
            print(f"  ✗ Failed to generate {xlsx_name}: {e}")
    
    return xlsx_files


def create_comparison_plots(consensus_file, original_file, tsv_files, out_dir, verbose=False):
    """Create comparison visualizations between filtered and original consensus."""
    
    if not PLOTLY_AVAILABLE:
        print("  ⚠ Plotly not available, skipping visualizations")
        return
    
    if verbose:
        print(f"\nGenerating comparison visualizations...")
    
    # Load both consensus files
    cmap_filtered = pyopenms.ConsensusMap()
    pyopenms.ConsensusXMLFile().load(consensus_file, cmap_filtered)
    
    cmap_original = pyopenms.ConsensusMap()
    pyopenms.ConsensusXMLFile().load(original_file, cmap_original)
    
    num_filtered = len([None for _ in cmap_filtered])
    num_original = len([None for _ in cmap_original])
    
    # 1. Feature count comparison
    fig1 = go.Figure(data=[
        go.Bar(name='Original Consensus', x=['Features'], y=[num_original], marker_color='#1f77b4'),
        go.Bar(name='Filtered Clusters', x=['Features'], y=[num_filtered], marker_color='#ff7f0e')
    ])
    fig1.update_layout(
        title='Consensus Features: Original vs Filtered',
        barmode='group',
        yaxis_title='Number of Features',
        template='plotly_white'
    )
    fig1.write_html(os.path.join(out_dir, 'comparison_01_feature_counts.html'))
    
    if verbose:
        print(f"  ✓ Feature count comparison: {num_original} → {num_filtered} ({100*num_filtered/num_original:.1f}%)")
    
    # 2. M/Z and RT distributions
    mzs_filtered = [ce.getMonoMZ() for ce in cmap_filtered]
    rts_filtered = [ce.getRT() for ce in cmap_filtered]
    intensities_filtered = [ce.getIntensity() for ce in cmap_filtered]
    
    mzs_original = [ce.getMonoMZ() for ce in cmap_original]
    rts_original = [ce.getRT() for ce in cmap_original]
    intensities_original = [ce.getIntensity() for ce in cmap_original]
    
    # M/Z distribution
    fig2 = go.Figure()
    fig2.add_trace(go.Histogram(x=mzs_original, name='Original', opacity=0.7, nbinsx=50))
    fig2.add_trace(go.Histogram(x=mzs_filtered, name='Filtered', opacity=0.7, nbinsx=50))
    fig2.update_layout(
        title='M/Z Distribution: Original vs Filtered',
        xaxis_title='M/Z',
        yaxis_title='Count',
        barmode='overlay',
        template='plotly_white'
    )
    fig2.write_html(os.path.join(out_dir, 'comparison_02_mz_distribution.html'))
    
    # RT distribution
    fig3 = go.Figure()
    fig3.add_trace(go.Histogram(x=rts_original, name='Original', opacity=0.7, nbinsx=50))
    fig3.add_trace(go.Histogram(x=rts_filtered, name='Filtered', opacity=0.7, nbinsx=50))
    fig3.update_layout(
        title='Retention Time Distribution: Original vs Filtered',
        xaxis_title='RT (seconds)',
        yaxis_title='Count',
        barmode='overlay',
        template='plotly_white'
    )
    fig3.write_html(os.path.join(out_dir, 'comparison_03_rt_distribution.html'))
    
    # Intensity distribution (log scale)
    fig4 = go.Figure()
    fig4.add_trace(go.Histogram(x=np.log10(np.array(intensities_original)+1), name='Original', opacity=0.7, nbinsx=50))
    fig4.add_trace(go.Histogram(x=np.log10(np.array(intensities_filtered)+1), name='Filtered', opacity=0.7, nbinsx=50))
    fig4.update_layout(
        title='Intensity Distribution (log10): Original vs Filtered',
        xaxis_title='log10(Intensity)',
        yaxis_title='Count',
        barmode='overlay',
        template='plotly_white'
    )
    fig4.write_html(os.path.join(out_dir, 'comparison_04_intensity_distribution.html'))
    
    # Scatter plot: M/Z vs RT colored by origin
    fig5 = go.Figure()
    fig5.add_trace(go.Scatter(
        x=mzs_original, y=rts_original, mode='markers',
        name='Original',
        marker=dict(size=5, opacity=0.5, color='#1f77b4')
    ))
    fig5.add_trace(go.Scatter(
        x=mzs_filtered, y=rts_filtered, mode='markers',
        name='Filtered',
        marker=dict(size=5, opacity=0.7, color='#ff7f0e')
    ))
    fig5.update_layout(
        title='M/Z vs RT: Original vs Filtered',
        xaxis_title='M/Z',
        yaxis_title='RT (seconds)',
        template='plotly_white'
    )
    fig5.write_html(os.path.join(out_dir, 'comparison_05_mz_vs_rt.html'))
    
    if verbose:
        print(f"  ✓ Generated 5 comparison visualizations")


def generate_summary(consensus_file, features_list, out_dir, num_original=None, verbose=False):
    """Generate summary statistics."""
    
    summary = {
        'consensus_file': consensus_file,
        'num_features': len(features_list),
        'output_dir': out_dir,
        'generated_at': datetime.now().isoformat(),
        'statistics': {}
    }
    
    if num_original:
        retention_rate = 100 * len(features_list) / num_original
        summary['num_original_features'] = num_original
        summary['retention_rate'] = retention_rate
        summary['filtered_features'] = num_original - len(features_list)
    
    # Write summary JSON
    with open(os.path.join(out_dir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Print summary
    print(f"\n{'='*70}")
    print(f"Processing Complete")
    print(f"{'='*70}\n")
    print(f"Consensus File: {os.path.basename(consensus_file)}")
    print(f"Consensus Features: {len(features_list)}")
    
    if num_original:
        print(f"Original Features: {num_original}")
        print(f"Retention Rate: {retention_rate:.1f}%")
        print(f"Filtered Out: {num_original - len(features_list)}")
    
    print(f"\nOutputs saved to: {out_dir}")
    print(f"\n✓ Processing complete\n")


if __name__ == "__main__":
    args = parse_args()
    
    # Check files exist
    if not os.path.exists(args.consensus):
        print(f"✗ ConsensusXML not found: {args.consensus}")
        sys.exit(1)
    
    if not os.path.exists(args.feature_tsvs_dir):
        print(f"✗ Feature TSV directory not found: {args.feature_tsvs_dir}")
        sys.exit(1)
    
    try:
        # Load feature TSV files
        print(f"Loading feature data from: {args.feature_tsvs_dir}")
        feature_dict, filenames = load_feature_tsvs(args.feature_tsvs_dir, args.verbose)
        print()
        
        # Process consensus
        features_list, header_full, header_reduced, header_minimal, global_headers, filenames = \
            process_consensusxml(args.consensus, feature_dict, filenames, args.verbose)
        
        # Create output directory
        os.makedirs(args.out_dir, exist_ok=True)
        
        # Write TSV outputs
        print(f"\nWriting TSV outputs to: {args.out_dir}")
        tsv_files = write_tsv_outputs(
            features_list, header_full, header_reduced, header_minimal, 
            global_headers, filenames, args.out_dir, args.verbose
        )
        
        # Generate XLSX reports
        print(f"\nGenerating XLSX reports...")
        xlsx_files = generate_xlsx_reports(tsv_files, args.out_dir, args.verbose)
        
        # Generate comparison if original provided
        if args.original_consensus and os.path.exists(args.original_consensus):
            cmap_original = pyopenms.ConsensusMap()
            pyopenms.ConsensusXMLFile().load(args.original_consensus, cmap_original)
            num_original = len([None for _ in cmap_original])
            
            create_comparison_plots(args.consensus, args.original_consensus, 
                                   tsv_files, args.out_dir, args.verbose)
        else:
            num_original = None
        
        # Generate summary
        generate_summary(args.consensus, features_list, args.out_dir, num_original, args.verbose)
        
    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

#!/usr/bin/env python
"""
Process ConsensusXML file and extract feature information into TSV format.

This script reads an OpenMS ConsensusXML file and extracts consensus feature 
information into structured TSV output files.

Usage:
    python process_consensusxml_standalone.py \
        -consensus <consensusXML_file> \
        -out_tsv <output_file>
"""

import argparse
import os
import sys
import json
import csv

import pandas as pd
import pyopenms
from datetime import datetime

csv.field_size_limit(sys.maxsize)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Process ConsensusXML file and extract feature information'
    )
    parser.add_argument('-consensus', required=True,
                       help='Input ConsensusXML file')
    parser.add_argument('-featurexmls_tsvs', default='',
                       help='Input feature TSV files (comma-separated). If provided, will merge data.')
    parser.add_argument('-out_tsv', default='consensus_features.tsv',
                       help='Output TSV file with consensus feature data')
    parser.add_argument('-out_json', default='',
                       help='Optional: Output JSON file with consensus feature data')
    parser.add_argument('-verbose', action='store_true',
                       help='Enable verbose output')
    
    return parser.parse_args()


def extract_feature_info(consensus_file, feature_tsvs_str='', verbose=False):
    """
    Extract consensus feature information from ConsensusXML file.
    
    Args:
        consensus_file: Path to consensusXML file
        feature_tsvs_str: Comma-separated list of feature TSV files (optional)
        verbose: Enable verbose output
        
    Returns:
        tuple: (features_list, column_headers, metadata)
    """
    if verbose:
        print(f"\n{'='*70}")
        print(f"Processing ConsensusXML: {consensus_file}")
        print(f"{'='*70}\n")
    
    # Load ConsensusXML
    print("Loading ConsensusXML file...")
    cmap = pyopenms.ConsensusMap()
    pyopenms.ConsensusXMLFile().load(consensus_file, cmap)
    
    if verbose:
        print(f"  Loaded ConsensusMap with {len([None for _ in cmap])} consensus features\n")
    
    # Get column headers (maps file indices to filenames)
    map_list = cmap.getColumnHeaders()
    if verbose:
        print(f"Input files in consensus ({len(map_list)} files):")
        for key, val in sorted(map_list.items()):
            print(f"  [{key}] {val.filename}")
        print()
    
    # Load feature TSV files if provided
    dict_of_single_features = {}
    filenames = []
    
    if feature_tsvs_str and feature_tsvs_str.strip():
        print(f"Loading feature TSV files...")
        for ftsv in feature_tsvs_str.split(','):
            ftsv = ftsv.strip()
            if not ftsv or not os.path.exists(ftsv):
                if verbose:
                    print(f"  ⚠ Skipping missing file: {ftsv}")
                continue
            
            filename_key = os.path.basename(ftsv).replace('_feature_data.tsv', '').replace('.tsv', '')
            try:
                # Try loading with literal_eval for list columns
                from ast import literal_eval
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
                # Fallback: load without converters
                df = pd.read_csv(ftsv, sep='\t')
            
            dict_of_single_features[filename_key] = df
            filenames.append(filename_key)
            if verbose:
                print(f"  ✓ Loaded {filename_key}: {len(df)} features")
        
        print()
    
    # Extract consensus features
    features = []
    
    print(f"Extracting consensus features...")
    for ce_idx, ce in enumerate(cmap):
        feature_data = {
            'openms_ceid': f"e_{ce.getUniqueId()}",
            'num_features': len(ce.getFeatureList()),
            'mz': ce.getMonoMZ(),
            'rt': ce.getRT(),
            'intensity': ce.getIntensity(),
            'charge': ce.getCharge(),
        }
        
        # Get feature list with file mapping
        feature_list = ce.getFeatureList()
        feature_data['files'] = []
        
        for fl in feature_list:
            file_idx = fl.getMapIndex()
            feature_id = fl.getUniqueId()
            
            # Map to filename from column headers
            file_info = map_list.get(file_idx)
            if file_info:
                filename = os.path.basename(file_info.filename).replace('_feature_data.tsv', '').replace('.featureXML', '')
                feature_data['files'].append({
                    'file_index': file_idx,
                    'filename': filename,
                    'feature_id': feature_id,
                    'mz': fl.getMZ(),
                    'rt': fl.getRT(),
                    'intensity': fl.getIntensity(),
                    'charge': fl.getCharge()
                })
        
        features.append(feature_data)
    
    if verbose:
        print(f"  ✓ Extracted {len(features)} consensus features\n")
    
    # Define output columns
    header_cols = [
        'openms_ceid',
        'num_features',
        'consensus_mz',
        'consensus_rt',
        'consensus_intensity',
        'consensus_charge',
        'files_involved'
    ]
    
    return features, header_cols, {
        'total_features': len(features),
        'num_input_files': len(map_list),
        'input_filenames': [os.path.basename(v.filename) for v in map_list.values()],
        'processed_at': datetime.now().isoformat()
    }


def write_tsv_output(features, headers, output_file, verbose=False):
    """Write features to TSV file."""
    print(f"Writing TSV output: {output_file}")
    
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        
        for feature in features:
            row = {
                'openms_ceid': feature['openms_ceid'],
                'num_features': feature['num_features'],
                'consensus_mz': f"{feature['mz']:.6f}",
                'consensus_rt': f"{feature['rt']:.2f}",
                'consensus_intensity': f"{feature['intensity']:.2f}",
                'consensus_charge': feature['charge'],
                'files_involved': json.dumps([f['filename'] for f in feature['files']])
            }
            writer.writerow(row)
    
    print(f"  ✓ Wrote {len(features)} consensus features to {output_file}\n")


def write_json_output(features, metadata, output_file, verbose=False):
    """Write features to JSON file."""
    print(f"Writing JSON output: {output_file}")
    
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    output_data = {
        'metadata': metadata,
        'features': []
    }
    
    for feature in features:
        output_data['features'].append({
            'openms_ceid': feature['openms_ceid'],
            'num_features': feature['num_features'],
            'mz': feature['mz'],
            'rt': feature['rt'],
            'intensity': feature['intensity'],
            'charge': feature['charge'],
            'files': feature['files']
        })
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"  ✓ Wrote {len(features)} consensus features to {output_file}\n")


def generate_summary(features, metadata, verbose=False):
    """Generate and print summary statistics."""
    print(f"{'='*70}")
    print(f"Summary Statistics")
    print(f"{'='*70}")
    print(f"\nInput Files: {metadata['num_input_files']}")
    for fname in metadata['input_filenames']:
        print(f"  - {fname}")
    
    print(f"\nConsensus Features: {metadata['total_features']}")
    
    # Calculate statistics
    mz_values = [f['mz'] for f in features]
    rt_values = [f['rt'] for f in features]
    intensity_values = [f['intensity'] for f in features]
    charge_values = [f['charge'] for f in features]
    features_per_consensus = [f['num_features'] for f in features]
    
    print(f"\nM/Z Statistics:")
    print(f"  Range: {min(mz_values):.4f} - {max(mz_values):.4f}")
    print(f"  Mean: {sum(mz_values)/len(mz_values):.4f}")
    
    print(f"\nRT Statistics (seconds):")
    print(f"  Range: {min(rt_values):.2f} - {max(rt_values):.2f}")
    print(f"  Mean: {sum(rt_values)/len(rt_values):.2f}")
    
    print(f"\nIntensity Statistics:")
    intensity_avg = sum(intensity_values) / len(intensity_values) if intensity_values else 0
    print(f"  Range: {min(intensity_values):.2e} - {max(intensity_values):.2e}")
    print(f"  Mean: {intensity_avg:.2e}")
    
    print(f"\nFeatures per Consensus Element:")
    print(f"  Min: {min(features_per_consensus)}")
    print(f"  Max: {max(features_per_consensus)}")
    print(f"  Mean: {sum(features_per_consensus)/len(features_per_consensus):.1f}")
    
    charge_dist = {}
    for c in charge_values:
        charge_dist[c] = charge_dist.get(c, 0) + 1
    print(f"\nCharge Distribution:")
    for charge in sorted(charge_dist.keys()):
        print(f"  +{charge}: {charge_dist[charge]}")
    
    print(f"\nProcessed: {metadata['processed_at']}\n")


if __name__ == "__main__":
    args = parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.consensus):
        print(f"✗ ConsensusXML file not found: {args.consensus}")
        sys.exit(1)
    
    try:
        # Extract features
        features, headers, metadata = extract_feature_info(
            args.consensus,
            args.featurexmls_tsvs,
            args.verbose
        )
        
        # Write outputs
        write_tsv_output(features, headers, args.out_tsv, args.verbose)
        
        if args.out_json:
            write_json_output(features, metadata, args.out_json, args.verbose)
        
        # Generate summary
        generate_summary(features, metadata, args.verbose)
        
        print(f"✓ Processing complete")
        
    except Exception as e:
        print(f"✗ Error processing ConsensusXML: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

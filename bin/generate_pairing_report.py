#!/usr/bin/env python3
"""
Generate summary report of feature pairing results.
"""

import argparse
import json
import os
from typing import Dict, List


def generate_pairing_report(paired_json: str, output_summary: str, output_stats: str):
    """Generate summary report and statistics from paired features."""
    
    with open(paired_json, 'r') as f:
        data = json.load(f)
    
    # Handle both old format (list) and new format (dict with statistics)
    if isinstance(data, dict) and 'paired_features' in data:
        paired_features = data['paired_features']
        pep_ident_stats = data.get('pep_ident_statistics', {})
        total_features_before_filter = data.get('total_features_before_filter', len(paired_features))
        total_features_after_filter = data.get('total_features_after_filter', sum(len(f.get('individual_features', [])) for f in paired_features))
    else:
        paired_features = data
        pep_ident_stats = {}
        total_features_before_filter = sum(len(f.get('individual_features', [])) if isinstance(f, dict) else 1 for f in paired_features)
        total_features_after_filter = total_features_before_filter
    
    print(f"  Processing {len(paired_features)} paired feature groups...")
    
    # Calculate statistics
    total_features = len(paired_features)
    match_scores = [f.get('match_score', 0) for f in paired_features]
    files_matched = [f.get('num_files_matched', 0) for f in paired_features]
    
    avg_match_score = sum(match_scores) / len(match_scores) if match_scores else 0
    max_match_score = max(match_scores) if match_scores else 0
    min_match_score = min(match_scores) if match_scores else 0
    
    # Count features matched across different numbers of files
    file_match_distribution = {}
    for f in paired_features:
        num_files = f.get('num_files_matched', 0)
        file_match_distribution[num_files] = file_match_distribution.get(num_files, 0) + 1
    
    # Count features with peptide identifications
    features_with_ident = sum(1 for f in paired_features if f.get('pep_ident'))
    
    # Generate summary text file
    with open(output_summary, 'w') as f:
        f.write("="*70 + "\n")
        f.write("FEATURE PAIRING SUMMARY REPORT\n")
        f.write("="*70 + "\n\n")
        
        f.write("OVERVIEW\n")
        f.write("-" * 70 + "\n")
        f.write(f"Total paired feature groups: {total_features}\n")
        f.write(f"Total features before filtering: {total_features_before_filter}\n")
        f.write(f"Total features after filtering: {total_features_after_filter}\n")
        f.write(f"Features with peptide identifications: {features_with_ident} ({100*features_with_ident/total_features:.1f}%)\n\n")
        
        f.write("MATCH SCORE STATISTICS\n")
        f.write("-" * 70 + "\n")
        f.write(f"Average match score: {avg_match_score:.4f}\n")
        f.write(f"Maximum match score: {max_match_score:.4f}\n")
        f.write(f"Minimum match score: {min_match_score:.4f}\n\n")
        
        f.write("FILE MATCHING DISTRIBUTION\n")
        f.write("-" * 70 + "\n")
        for num_files in sorted(file_match_distribution.keys()):
            count = file_match_distribution[num_files]
            percentage = 100 * count / total_features
            f.write(f"Matched across {num_files} file(s): {count} ({percentage:.1f}%)\n")
        
        # Add peptide identification statistics if available
        if pep_ident_stats:
            f.write("\nPEPTIDE IDENTIFICATION STATISTICS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Groups with all features sharing at least one common pep_ident: {pep_ident_stats.get('groups_with_common_pep_idents', 0)}\n")
            f.write(f"Groups with all features having mismatching pep_idents (no common ones): {pep_ident_stats.get('groups_with_mismatching_pep_idents', 0)}\n")
            f.write(f"Groups with no pep_idents: {pep_ident_stats.get('groups_with_no_pep_idents', 0)}\n")
            f.write(f"Groups with partial pep_idents (some features have none): {pep_ident_stats.get('groups_with_partial_matching', 0)}\n")
            if pep_ident_stats.get('pep_ident_match_distribution'):
                f.write("\nMatching distribution:\n")
                for key, count in pep_ident_stats['pep_ident_match_distribution'].items():
                    f.write(f"  {key}: {count}\n")
        
        f.write("\n" + "="*70 + "\n")
    
    # Generate statistics CSV file
    with open(output_stats, 'w') as f:
        f.write("metric,value\n")
        f.write(f"total_features_before_filter,{total_features_before_filter}\n")
        f.write(f"total_features_after_filter,{total_features_after_filter}\n")
        f.write(f"total_features,{total_features}\n")
        f.write(f"features_with_identifications,{features_with_ident}\n")
        f.write(f"average_match_score,{avg_match_score:.6f}\n")
        f.write(f"max_match_score,{max_match_score:.6f}\n")
        f.write(f"min_match_score,{min_match_score:.6f}\n")
        
        for num_files in sorted(file_match_distribution.keys()):
            count = file_match_distribution[num_files]
            f.write(f"matched_across_{num_files}_files,{count}\n")
        
        # Write peptide identification statistics
        if pep_ident_stats:
            f.write(f"groups_with_common_pep_idents,{pep_ident_stats.get('groups_with_common_pep_idents', 0)}\n")
            f.write(f"groups_with_mismatching_pep_idents,{pep_ident_stats.get('groups_with_mismatching_pep_idents', 0)}\n")
            f.write(f"groups_with_no_pep_idents,{pep_ident_stats.get('groups_with_no_pep_idents', 0)}\n")
            f.write(f"groups_with_partial_matching,{pep_ident_stats.get('groups_with_partial_matching', 0)}\n")
            
            for key, count in pep_ident_stats.get('pep_ident_match_distribution', {}).items():
                f.write(f"pep_ident_distribution_{key},{count}\n")
    
    print(f"✓ Summary report generated")
    print(f"✓ Statistics CSV generated")


def main():
    parser = argparse.ArgumentParser(
        description="Generate summary report of feature pairing results"
    )
    parser.add_argument("--paired_json", required=True, help="Path to paired features JSON file")
    parser.add_argument("--output_summary", required=True, help="Output summary text file path")
    parser.add_argument("--output_stats", required=True, help="Output statistics CSV file path")
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print(f"Generating Pairing Report")
    print(f"{'='*70}")
    
    generate_pairing_report(args.paired_json, args.output_summary, args.output_stats)
    
    print(f"\n✓ Report saved to {os.path.basename(args.output_summary)}")
    print(f"✓ Statistics saved to {os.path.basename(args.output_stats)}")


if __name__ == "__main__":
    main()

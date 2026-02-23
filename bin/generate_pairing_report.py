#!/usr/bin/env python3
"""
Generate summary report of feature pairing results.
"""

import argparse
import json
import os
from typing import Dict, List
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import base64
from io import BytesIO
import pandas as pd


def generate_pairing_report(paired_json: str, output_summary: str, output_stats: str, edges_json: str = None, output_html: str = None, graph_composition_json: str = None):
    """Generate summary report and statistics from paired features."""
    
    with open(paired_json, 'r') as f:
        data = json.load(f)
    
    # Load graph composition if provided
    graph_composition = None
    if graph_composition_json and os.path.exists(graph_composition_json):
        try:
            with open(graph_composition_json, 'r') as f:
                graph_composition = json.load(f)
            print(f"  Loaded graph composition from {os.path.basename(graph_composition_json)}")
        except Exception as e:
            print(f"  Warning: Could not load graph composition: {e}")
    
    # Check if this is an edges file format (has numeric keys with lists of edges)
    is_edges_format = isinstance(data, dict) and all(isinstance(v, list) for v in data.values())
    
    # Extract pairing parameters
    pairing_parameters = {}
    cutoff_params = {}
    
    if is_edges_format:
        # Handle edges format directly
        paired_features = []
        for charge_key, edges_list in data.items():
            paired_features.append({'num_files_matched': 0, 'pep_ident': None})
        pep_ident_stats = {}
        total_features_before_filter = len(paired_features)
        total_features_after_filter = len(paired_features)
        # Use the paired_json as edges_json if edges_json not provided
        if not edges_json:
            edges_json = paired_json
    # Handle both old format (list) and new format (dict with statistics)
    elif isinstance(data, dict) and 'paired_features' in data:
        paired_features = data['paired_features']
        pep_ident_stats = data.get('pep_ident_statistics', {})
        total_features_before_filter = data.get('total_features_before_filter', len(paired_features))
        total_features_after_filter = data.get('total_features_after_filter', sum(len(f.get('individual_features', [])) for f in paired_features))
        pairing_parameters = data.get('pairing_parameters', {})
        cutoff_params = data.get('cutoff_params', {})
    else:
        paired_features = data
        pep_ident_stats = {}
        total_features_before_filter = sum(len(f.get('individual_features', [])) if isinstance(f, dict) else 1 for f in paired_features)
        total_features_after_filter = total_features_before_filter
    
    print(f"  Processing {len(paired_features)} paired feature groups...")
    
    # Extract distances from edges JSON if provided
    distances = []
    if edges_json and os.path.exists(edges_json):
        try:
            with open(edges_json, 'r') as f:
                edges_data = json.load(f)
            # Flatten the nested structure to collect all distances
            for charge_key, edges_list in edges_data.items():
                for edge in edges_list:
                    if 'distance' in edge:
                        distances.append(edge['distance'])
            print(f"  Extracted {len(distances)} distances from edges file")
        except Exception as e:
            print(f"  Warning: Could not extract distances from edges file: {e}")
    
    # Calculate statistics
    total_features = len(paired_features)
    files_matched = [f.get('num_files_matched', 0) for f in paired_features if isinstance(f, dict)]
    
    # Extract intra-group distances from new match structure (distances_to_group_features)
    intra_group_distances = []
    for match in paired_features:
        if isinstance(match, dict) and 'individual_features' in match:
            for feature in match['individual_features']:
                distances_dict = feature.get('distances_to_group_features', {})
                intra_group_distances.extend(distances_dict.values())
    
    avg_intra_distance = sum(intra_group_distances) / len(intra_group_distances) if intra_group_distances else 0
    max_intra_distance = max(intra_group_distances) if intra_group_distances else 0
    min_intra_distance = min(intra_group_distances) if intra_group_distances else 0
    
    # Count features matched across different numbers of files
    file_match_distribution = {}
    for f in paired_features:
        if isinstance(f, dict):
            num_files = f.get('num_files_matched', 0)
            file_match_distribution[num_files] = file_match_distribution.get(num_files, 0) + 1
    
    # Count features with peptide identifications
    features_with_ident = sum(1 for f in paired_features if isinstance(f, dict) and f.get('pep_ident'))
    
    # Generate HTML report if requested
    if output_html:
        generate_html_report(
            output_html,
            total_features,
            total_features_before_filter,
            total_features_after_filter,
            features_with_ident,
            avg_intra_distance,
            max_intra_distance,
            min_intra_distance,
            file_match_distribution,
            pep_ident_stats,
            distances,
            graph_composition,
            intra_group_distances,
            pairing_parameters,
            cutoff_params
        )
    
    # Generate summary text file
    with open(output_summary, 'w') as f:
        f.write("="*70 + "\n")
        f.write("FEATURE PAIRING SUMMARY REPORT\n")
        f.write("="*70 + "\n\n")
        
        # Write pairing parameters at the beginning
        if pairing_parameters:
            f.write("PAIRING PARAMETERS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Pairing Method: {pairing_parameters.get('pairing_method', 'unknown')}\n")
            f.write(f"Best Match Only: {pairing_parameters.get('best_match_only', False)}\n")
            f.write(f"Normalize Coordinates: {pairing_parameters.get('normalize_coordinates', True)}\n")
            f.write(f"Distance Calc Before Scaling: {pairing_parameters.get('distance_calc_before_scaling', False)}\n")
            f.write(f"Match Cutoff: {pairing_parameters.get('match_cutoff', 0.0)}\n")
            f.write(f"Input Files: {pairing_parameters.get('input_files', 'unknown')}\n")
            f.write(f"Skip Matchfinder: {pairing_parameters.get('skip_matchfinder', False)}\n\n")
        
        if cutoff_params:
            f.write("FILTERING PARAMETERS\n")
            f.write("-" * 70 + "\n")
            if cutoff_params.get('mz_cutoff') is not None:
                f.write(f"M/Z Cutoff: {cutoff_params.get('mz_cutoff')}\n")
            if cutoff_params.get('rt_cutoff') is not None:
                f.write(f"RT Cutoff: {cutoff_params.get('rt_cutoff')}\n")
            if cutoff_params.get('edges_cutoff') is not None:
                f.write(f"Edges Distance Cutoff: {cutoff_params.get('edges_cutoff')}\n")
            if not any(v is not None for v in cutoff_params.values()):
                f.write("No filtering applied\n")
            f.write("\n")
        
        f.write("OVERVIEW\n")
        f.write("-" * 70 + "\n")
        f.write(f"Total paired feature groups: {total_features}\n")
        f.write(f"Total features before filtering: {total_features_before_filter}\n")
        f.write(f"Total features after filtering: {total_features_after_filter}\n")
        f.write(f"Features with peptide identifications: {features_with_ident} ({100*features_with_ident/total_features:.1f}%)\n\n")
        
        f.write("INTRA-GROUP DISTANCE STATISTICS\n")
        f.write("-" * 70 + "\n")
        f.write(f"Average intra-group distance: {avg_intra_distance:.4f}\n")
        f.write(f"Maximum intra-group distance: {max_intra_distance:.4f}\n")
        f.write(f"Minimum intra-group distance: {min_intra_distance:.4f}\n")
        f.write(f"Total intra-group pairwise distances: {len(intra_group_distances)}\n\n")
        
        if distances:
            f.write("DISTANCE STATISTICS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Total paired edges: {len(distances)}\n")
            f.write(f"Average distance: {sum(distances)/len(distances):.4f}\n")
            f.write(f"Median distance: {sorted(distances)[len(distances)//2]:.4f}\n")
            f.write(f"Min distance: {min(distances):.4f}\n")
            f.write(f"Max distance: {max(distances):.4f}\n\n")
        
        # Add graph composition statistics if available
        if graph_composition:
            comp_summary = graph_composition.get('composition_summary', {})
            comp_dist = graph_composition.get('component_composition', {})
            if comp_summary:
                f.write("GRAPH COMPOSITION STATISTICS\n")
                f.write("-" * 70 + "\n")
                f.write(f"Total connected components: {comp_summary.get('total_components')}\n")
                f.write(f"Total features in components: {comp_summary.get('total_features')}\n")
                f.write(f"Average features per component: {comp_summary.get('avg_features_per_component', 0):.2f}\n")
                f.write(f"Average files per component: {comp_summary.get('avg_files_per_component', 0):.2f}\n")
                f.write(f"Max features in any component: {comp_summary.get('max_features_in_component')}\n")
                f.write(f"Max files in any component: {comp_summary.get('max_files_in_component')}\n")
                
                # Show distribution by composition
                if comp_dist:
                    f.write("\nComponent distribution by (features, files):\n")
                    for comp_key in sorted(comp_dist.keys(), key=lambda x: tuple(map(int, x.split('_')))):
                        count = comp_dist[comp_key]
                        parts = comp_key.split('_')
                        f.write(f"  {parts[0]} features from {parts[1]} file(s): {count} component{'s' if count != 1 else ''}\n")
                f.write("\n")
        
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
        f.write(f"average_intra_group_distance,{avg_intra_distance:.6f}\n")
        f.write(f"max_intra_group_distance,{max_intra_distance:.6f}\n")
        f.write(f"min_intra_group_distance,{min_intra_distance:.6f}\n")
        f.write(f"total_intra_group_distances,{len(intra_group_distances)}\n")
        
        # Write distance statistics
        if distances:
            f.write(f"total_paired_edges,{len(distances)}\n")
            f.write(f"average_distance,{sum(distances)/len(distances):.6f}\n")
            f.write(f"median_distance,{sorted(distances)[len(distances)//2]:.6f}\n")
            f.write(f"min_distance,{min(distances):.6f}\n")
            f.write(f"max_distance,{max(distances):.6f}\n")
            f.write(f"q1_distance,{sorted(distances)[len(distances)//4]:.6f}\n")
            f.write(f"q3_distance,{sorted(distances)[3*len(distances)//4]:.6f}\n")
        
        for num_files in sorted(file_match_distribution.keys()):
            count = file_match_distribution[num_files]
            f.write(f"matched_across_{num_files}_files,{count}\n")
        
        # Write graph composition statistics
        if graph_composition:
            comp_summary = graph_composition.get('composition_summary', {})
            if comp_summary:
                f.write(f"total_components,{comp_summary.get('total_components')}\n")
                f.write(f"total_features_in_components,{comp_summary.get('total_features')}\n")
                f.write(f"avg_features_per_component,{comp_summary.get('avg_features_per_component', 0):.6f}\n")
                f.write(f"avg_files_per_component,{comp_summary.get('avg_files_per_component', 0):.6f}\n")
                f.write(f"max_features_in_component,{comp_summary.get('max_features_in_component')}\n")
                f.write(f"max_files_in_component,{comp_summary.get('max_files_in_component')}\n")
                
                # Write individual composition distribution
                comp_dist = graph_composition.get('component_composition', {})
                for comp_key, count in comp_dist.items():
                    f.write(f"component_composition_{comp_key},{count}\n")
        
        # Write peptide identification statistics
        if pep_ident_stats:
            f.write(f"groups_with_common_pep_idents,{pep_ident_stats.get('groups_with_common_pep_idents', 0)}\n")
            f.write(f"groups_with_mismatching_pep_idents,{pep_ident_stats.get('groups_with_mismatching_pep_idents', 0)}\n")
            f.write(f"groups_with_no_pep_idents,{pep_ident_stats.get('groups_with_no_pep_idents', 0)}\n")
            f.write(f"groups_with_partial_matching,{pep_ident_stats.get('groups_with_partial_matching', 0)}\n")
            
            for key, count in pep_ident_stats.get('pep_ident_match_distribution', {}).items():
                f.write(f"pep_ident_distribution_{key},{count}\n")
    
    print(f"✓ Summary report generated")
    if output_html:
        print(f"✓ HTML report generated")
    print(f"✓ Statistics CSV generated")


def create_distance_histogram_base64(distances: List[float]) -> str:
    """Create a distance distribution histogram and return as base64 encoded image."""
    if not distances or len(distances) == 0:
        return None
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create histogram
    ax.hist(distances, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax.set_xlabel('Distance', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title('Distance Distribution Among Paired Features', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    # Add statistics as text
    mean_dist = sum(distances) / len(distances)
    median_dist = sorted(distances)[len(distances)//2]
    ax.axvline(mean_dist, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_dist:.3f}')
    ax.axvline(median_dist, color='green', linestyle='--', linewidth=2, label=f'Median: {median_dist:.3f}')
    ax.legend(fontsize=10)
    
    # Convert to base64
    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=100, bbox_inches='tight')
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode()
    plt.close(fig)
    
    return image_base64


def generate_html_report(
    output_html: str,
    total_features: int,
    total_features_before_filter: int,
    total_features_after_filter: int,
    features_with_ident: int,
    avg_intra_distance: float,
    max_intra_distance: float,
    min_intra_distance: float,
    file_match_distribution: Dict,
    pep_ident_stats: Dict,
    distances: List[float],
    graph_composition: Dict = None,
    intra_group_distances: List[float] = None,
    pairing_parameters: Dict = None,
    cutoff_params: Dict = None
):
    """Generate a beautiful HTML report with visualizations."""
    
    # Create histogram
    histogram_base64 = create_distance_histogram_base64(distances) if distances else None
    
    # Calculate distance statistics
    dist_stats = {}
    if distances:
        dist_stats = {
            'count': len(distances),
            'mean': sum(distances) / len(distances),
            'median': sorted(distances)[len(distances)//2],
            'min': min(distances),
            'max': max(distances),
            'q1': sorted(distances)[len(distances)//4],
            'q3': sorted(distances)[3*len(distances)//4],
        }
    
    # Generate pairing parameters HTML
    params_html = ""
    if pairing_parameters:
        params_html = f"""
        <div class="stats-section">
            <h3>Pairing Parameters</h3>
            <div class="stats-grid">
                <div class="stat-box">
                    <div class="stat-value">{pairing_parameters.get('pairing_method', 'unknown').upper()}</div>
                    <div class="stat-label">Pairing Method</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{'YES' if pairing_parameters.get('best_match_only', False) else 'NO'}</div>
                    <div class="stat-label">Best Match Only</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{'YES' if pairing_parameters.get('normalize_coordinates', True) else 'NO'}</div>
                    <div class="stat-label">Normalize Coordinates</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{'YES' if pairing_parameters.get('distance_calc_before_scaling', False) else 'NO'}</div>
                    <div class="stat-label">Distance Before Scaling</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{pairing_parameters.get('match_cutoff', 0.0)}</div>
                    <div class="stat-label">Match Cutoff</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{pairing_parameters.get('input_files', 'unknown')}</div>
                    <div class="stat-label">Input Files</div>
                </div>
            </div>
        </div>
        """
    
    # Generate filtering parameters HTML
    filter_html = ""
    if cutoff_params and any(v is not None for v in cutoff_params.values()):
        filter_items_html = ""
        if cutoff_params.get('mz_cutoff') is not None:
            filter_items_html += f"<div>M/Z Cutoff: <strong>{cutoff_params.get('mz_cutoff')}</strong></div>"
        if cutoff_params.get('rt_cutoff') is not None:
            filter_items_html += f"<div>RT Cutoff: <strong>{cutoff_params.get('rt_cutoff')}</strong></div>"
        if cutoff_params.get('edges_cutoff') is not None:
            filter_items_html += f"<div>Edges Distance Cutoff: <strong>{cutoff_params.get('edges_cutoff')}</strong></div>"
        
        if filter_items_html:
            filter_html = f"""
            <div class="stats-section">
                <h3>Filtering Parameters</h3>
                <div style="background: #f8f9fa; padding: 15px; border-radius: 8px; font-size: 1.05em;">
                    {filter_items_html}
                </div>
            </div>
            """
    
    # Generate file distribution HTML
    file_dist_html = ""
    for num_files in sorted(file_match_distribution.keys()):
        count = file_match_distribution[num_files]
        percentage = 100 * count / total_features
        bar_width = min(95, percentage)
        file_dist_html += f"""
        <div class="distribution-row">
            <div class="dist-label">{num_files} file(s)</div>
            <div class="dist-bar-container">
                <div class="dist-bar" style="width: {bar_width}%"></div>
            </div>
            <div class="dist-value">{count} ({percentage:.1f}%)</div>
        </div>
        """
    
    # Generate identification statistics HTML
    ident_html = ""
    if pep_ident_stats:
        ident_html = f"""
        <div class="stats-section">
            <h3>Peptide Identification Statistics</h3>
            <div class="stats-grid">
                <div class="stat-box">
                    <div class="stat-value">{pep_ident_stats.get('groups_with_common_pep_idents', 0)}</div>
                    <div class="stat-label">Common Identifications</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{pep_ident_stats.get('groups_with_mismatching_pep_idents', 0)}</div>
                    <div class="stat-label">Mismatching Identifications</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{pep_ident_stats.get('groups_with_no_pep_idents', 0)}</div>
                    <div class="stat-label">No Identifications</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{pep_ident_stats.get('groups_with_partial_matching', 0)}</div>
                    <div class="stat-label">Partial Matching</div>
                </div>
            </div>
        </div>
        """
    
    # Generate graph composition HTML
    comp_html = ""
    if graph_composition:
        comp_summary = graph_composition.get('composition_summary', {})
        comp_dist = graph_composition.get('component_composition', {})
        
        if comp_summary:
            # Create composition distribution table
            comp_rows_html = ""
            for comp_key in sorted(comp_dist.keys(), key=lambda x: tuple(map(int, x.split('_')))):
                count = comp_dist[comp_key]
                parts = comp_key.split('_')
                num_features = parts[0]
                num_files = parts[1]
                comp_rows_html += f"""
                <div class="comp-row">
                    <span class="comp-label">{num_features} features from {num_files} file(s)</span>
                    <span class="comp-count">{count} component{'s' if count != 1 else ''}</span>
                </div>
                """
            
            comp_html = f"""
            <div class="stats-section">
                <h3>Graph Composition</h3>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{comp_summary.get('total_components', 0)}</div>
                        <div class="stat-label">Total Components</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{comp_summary.get('avg_features_per_component', 0):.2f}</div>
                        <div class="stat-label">Avg Features/Component</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{comp_summary.get('avg_files_per_component', 0):.2f}</div>
                        <div class="stat-label">Avg Files/Component</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{comp_summary.get('max_features_in_component', 0)}</div>
                        <div class="stat-label">Max Features/Component</div>
                    </div>
                </div>
                <h4 style="margin-top: 20px; color: #333;">Component Distribution</h4>
                <div class="composition-list">
                    {comp_rows_html}
                </div>
            </div>
            """
    
    # Generate distance statistics section
    dist_section_html = ""
    if dist_stats:
        dist_section_html = f"""
        <div class="stats-section">
            <h3>Distance Statistics</h3>
            <div class="stats-grid">
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['count']}</div>
                    <div class="stat-label">Total Paired Edges</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['mean']:.4f}</div>
                    <div class="stat-label">Mean Distance</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['median']:.4f}</div>
                    <div class="stat-label">Median Distance</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['min']:.4f}</div>
                    <div class="stat-label">Min Distance</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['max']:.4f}</div>
                    <div class="stat-label">Max Distance</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{dist_stats['q3']-dist_stats['q1']:.4f}</div>
                    <div class="stat-label">IQR (Q3-Q1)</div>
                </div>
            </div>
        </div>
        """
    
    ident_percentage = 100 * features_with_ident / total_features if total_features > 0 else 0
    
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Feature Pairing Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 12px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}
        
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        
        .header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        
        .content {{
            padding: 40px;
        }}
        
        .stats-section {{
            margin-bottom: 40px;
        }}
        
        .stats-section h3 {{
            color: #333;
            margin-bottom: 20px;
            font-size: 1.5em;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .stat-box {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            transition: transform 0.3s ease, box-shadow 0.3s ease;
        }}
        
        .stat-box:hover {{
            transform: translateY(-5px);
            box-shadow: 0 8px 25px rgba(0,0,0,0.15);
        }}
        
        .stat-value {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
            margin-bottom: 10px;
        }}
        
        .stat-label {{
            font-size: 0.9em;
            color: #555;
            font-weight: 500;
        }}
        
        .overview-box {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 40px;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
        }}
        
        .overview-item {{
            text-align: center;
        }}
        
        .overview-value {{
            font-size: 2.5em;
            font-weight: bold;
            margin-bottom: 10px;
        }}
        
        .overview-label {{
            font-size: 0.95em;
            opacity: 0.9;
        }}
        
        .distribution-section {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        
        .distribution-row {{
            display: flex;
            align-items: center;
            margin-bottom: 15px;
            gap: 15px;
        }}
        
        .dist-label {{
            min-width: 100px;
            font-weight: 500;
            color: #333;
        }}
        
        .dist-bar-container {{
            flex: 1;
            height: 30px;
            background: #ddd;
            border-radius: 5px;
            overflow: hidden;
        }}
        
        .dist-bar {{
            height: 100%;
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            transition: width 0.3s ease;
        }}
        
        .dist-value {{
            min-width: 120px;
            text-align: right;
            font-weight: 500;
            color: #666;
        }}
        
        .composition-list {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin-top: 15px;
        }}
        
        .comp-row {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 10px 0;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        .comp-row:last-child {{
            border-bottom: none;
        }}
        
        .comp-label {{
            font-weight: 500;
            color: #333;
        }}
        
        .comp-count {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.9em;
            font-weight: 500;
        }}
        
        .histogram-container {{
            text-align: center;
            margin: 40px 0;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 10px;
        }}
        
        .histogram-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }}
        
        .histogram-container h3 {{
            margin-bottom: 20px;
            color: #333;
        }}
        
        .footer {{
            background: #f8f9fa;
            padding: 20px 40px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
            border-top: 1px solid #ddd;
        }}
        
        @media (max-width: 768px) {{
            .stats-grid {{
                grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            }}
            
            .header h1 {{
                font-size: 1.8em;
            }}
            
            .distribution-row {{
                flex-direction: column;
                align-items: flex-start;
            }}
            
            .dist-bar-container {{
                width: 100%;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>📊 Feature Pairing Report</h1>
            <p>Comprehensive Analysis of Paired Features</p>
        </div>
        
        <div class="content">
            <!-- Overview Section -->
            <div class="overview-box">
                <div class="overview-item">
                    <div class="overview-value">{total_features:,}</div>
                    <div class="overview-label">Paired Groups</div>
                </div>
                <div class="overview-item">
                    <div class="overview-value">{features_with_ident:,}</div>
                    <div class="overview-label">With Identifications</div>
                </div>
                <div class="overview-item">
                    <div class="overview-value">{ident_percentage:.1f}%</div>
                    <div class="overview-label">Identification Rate</div>
                </div>
                <div class="overview-item">
                    <div class="overview-value">{avg_intra_distance:.4f}</div>
                    <div class="overview-label">Avg Intra-group Distance</div>
                </div>
            </div>
            
            <!-- Pairing Parameters Section -->
            {params_html}
            
            <!-- Filtering Parameters Section -->
            {filter_html}
            
            <!-- Feature Summary Section -->
            <div class="stats-section">
                <h3>Feature Summary</h3>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{total_features_before_filter:,}</div>
                        <div class="stat-label">Features Before Filtering</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{total_features_after_filter:,}</div>
                        <div class="stat-label">Features After Filtering</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{max_intra_distance:.4f}</div>
                        <div class="stat-label">Max Intra-group Distance</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{min_intra_distance:.4f}</div>
                        <div class="stat-label">Min Intra-group Distance</div>
                    </div>
                </div>
            </div>
            
            <!-- Distance Statistics -->
            {dist_section_html}
            
            <!-- Distance Histogram -->
            {f'<div class="histogram-container"><h3>Distance Distribution</h3><img src="data:image/png;base64,{histogram_base64}" alt="Distance Distribution Histogram"></div>' if histogram_base64 else ''}
            
            <!-- File Matching Distribution -->
            <div class="stats-section">
                <h3>File Matching Distribution</h3>
                <div class="distribution-section">
                    {file_dist_html}
                </div>
            </div>
            
            <!-- Graph Composition Statistics -->
            {comp_html}
            
            <!-- Peptide Identification Statistics -->
            {ident_html}
        </div>
        
        <div class="footer">
            Generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')} | Feature Pairing Analysis Report
        </div>
    </div>
</body>
</html>
"""
    
    with open(output_html, 'w') as f:
        f.write(html_content)



def main():
    parser = argparse.ArgumentParser(
        description="Generate summary report of feature pairing results"
    )
    parser.add_argument("--paired_json", required=True, help="Path to paired features JSON file")
    parser.add_argument("--output_summary", required=True, help="Output summary text file path")
    parser.add_argument("--output_stats", required=True, help="Output statistics CSV file path")
    parser.add_argument("--edges_json", default=None, help="Path to edges JSON file (for distance statistics and histogram)")
    parser.add_argument("--output_html", default=None, help="Output HTML report file path")
    parser.add_argument("--graph_composition_json", default=None, help="Path to graph composition JSON file (from build_network_graph --output_composition)")
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print(f"Generating Pairing Report")
    print(f"{'='*70}")
    
    generate_pairing_report(args.paired_json, args.output_summary, args.output_stats, args.edges_json, args.output_html, args.graph_composition_json)
    
    print(f"\n✓ Report saved to {os.path.basename(args.output_summary)}")
    print(f"✓ Statistics saved to {os.path.basename(args.output_stats)}")
    if args.output_html:
        print(f"✓ HTML report saved to {os.path.basename(args.output_html)}")


if __name__ == "__main__":
    main()

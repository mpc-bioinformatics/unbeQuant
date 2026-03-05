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


def generate_pairing_report(paired_json: str, output_summary: str, output_stats: str, edges_json: str = None, output_html: str = None, graph_composition_json: str = None, clusters_json: str = None, parameters: Dict = None):
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
    
    # Load cluster data if provided
    clusters_data = None
    if clusters_json and os.path.exists(clusters_json):
        try:
            with open(clusters_json, 'r') as f:
                clusters_data = json.load(f)
            print(f"  Loaded cluster data from {os.path.basename(clusters_json)}")
        except Exception as e:
            print(f"  Warning: Could not load cluster data: {e}")
    
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
    # and analyze file uniqueness and multiplicity within each group
    file_match_distribution = {}
    file_match_uniqueness = {}  # Track file uniqueness and multiplicity
    
    for f in paired_features:
        if isinstance(f, dict):
            num_files = f.get('num_files_matched', 0)
            file_match_distribution[num_files] = file_match_distribution.get(num_files, 0) + 1
            
            # Analyze file uniqueness and multiplicity for this group
            if num_files not in file_match_uniqueness:
                file_match_uniqueness[num_files] = {
                    'all_unique': 0,  # All files unique (one feature per file)
                    'with_duplicates': 0,  # At least one file has multiple features
                    'total': 0,
                    'file_composition': {},  # file_list -> count
                    'multiplicity_distribution': {},  # e.g., "1x2,1x1" -> count (1 file with 2, 1 with 1)
                    'file_appearance': {}  # filename -> count (how many times each file appears in matches)
                }
            
            file_match_uniqueness[num_files]['total'] += 1
            
            # Count file occurrences in this group
            individual_features = f.get('individual_features', [])
            if individual_features:
                file_counts = {}  # filename -> count
                for feat in individual_features:
                    fname = feat.get('filename')
                    file_counts[fname] = file_counts.get(fname, 0) + 1
                
                # Track file appearance frequency
                for fname in file_counts.keys():
                    if fname not in file_match_uniqueness[num_files]['file_appearance']:
                        file_match_uniqueness[num_files]['file_appearance'][fname] = 0
                    file_match_uniqueness[num_files]['file_appearance'][fname] += 1
                
                # Check if all files are unique (each appears exactly once)
                multiplicity_list = sorted(file_counts.values(), reverse=True)
                is_all_unique = all(count == 1 for count in multiplicity_list)
                
                if is_all_unique:
                    file_match_uniqueness[num_files]['all_unique'] += 1
                else:
                    file_match_uniqueness[num_files]['with_duplicates'] += 1
                
                # Track file composition
                files_in_group = sorted(set(fname for fname in file_counts.keys()))
                file_comp_str = ', '.join(files_in_group)
                if file_comp_str not in file_match_uniqueness[num_files]['file_composition']:
                    file_match_uniqueness[num_files]['file_composition'][file_comp_str] = 0
                file_match_uniqueness[num_files]['file_composition'][file_comp_str] += 1
                
                # Track multiplicity pattern (e.g., "2x2,1x1" for 2 files with 2 features, 1 with 1)
                multiplicity_pattern = ','.join([f"{count}x{multiplicity_list.count(count)}" for count in sorted(set(multiplicity_list), reverse=True)])
                if multiplicity_pattern not in file_match_uniqueness[num_files]['multiplicity_distribution']:
                    file_match_uniqueness[num_files]['multiplicity_distribution'][multiplicity_pattern] = 0
                file_match_uniqueness[num_files]['multiplicity_distribution'][multiplicity_pattern] += 1
    
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
            file_match_uniqueness,
            pep_ident_stats,
            distances,
            graph_composition,
            clusters_data,
            intra_group_distances,
            pairing_parameters,
            cutoff_params,
            parameters
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
    file_match_uniqueness: Dict,
    pep_ident_stats: Dict,
    distances: List[float],
    graph_composition: Dict = None,
    clusters_data: Dict = None,
    intra_group_distances: List[float] = None,
    pairing_parameters: Dict = None,
    cutoff_params: Dict = None,
    parameters: Dict = None
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
    
    # Generate file distribution HTML with clickable bars
    file_dist_html = """
    <div id="fileDistChart" class="distribution-section">
    """
    
    # Generate data for JavaScript
    file_dist_data = {}
    for num_files in sorted(file_match_distribution.keys()):
        count = file_match_distribution[num_files]
        percentage = 100 * count / total_features
        
        uniqueness_info = file_match_uniqueness.get(num_files, {})
        all_unique_count = uniqueness_info.get('all_unique', 0)
        with_duplicates_count = uniqueness_info.get('with_duplicates', 0)
        file_compositions = uniqueness_info.get('file_composition', {})
        multiplicity_dist = uniqueness_info.get('multiplicity_distribution', {})
        
        # Create composition breakdown
        composition_breakdown = []
        for file_comp, comp_count in sorted(file_compositions.items(), key=lambda x: -x[1]):
            composition_breakdown.append({
                'composition': file_comp,
                'count': comp_count,
                'percentage': 100 * comp_count / count if count > 0 else 0
            })
        
        # Create multiplicity breakdown
        multiplicity_breakdown = []
        for mult_pattern, mult_count in sorted(multiplicity_dist.items(), key=lambda x: -x[1]):
            multiplicity_breakdown.append({
                'pattern': mult_pattern,
                'count': mult_count,
                'percentage': 100 * mult_count / count if count > 0 else 0
            })
        
        # Create file appearance breakdown
        file_appearance_dist = uniqueness_info.get('file_appearance', {})
        file_appearance_breakdown = []
        for fname, appearance_count in sorted(file_appearance_dist.items(), key=lambda x: -x[1]):
            file_appearance_breakdown.append({
                'filename': fname,
                'count': appearance_count,
                'percentage': 100 * appearance_count / count if count > 0 else 0
            })
        
        file_dist_data[num_files] = {
            'total': count,
            'percentage': percentage,
            'all_unique': all_unique_count,
            'with_duplicates': with_duplicates_count,
            'composition': composition_breakdown,
            'multiplicity': multiplicity_breakdown,
            'file_appearance': file_appearance_breakdown
        }
        
        bar_width = min(95, percentage)
        file_dist_html += f"""
        <div class="distribution-row clickable-bar" onclick="showFileDistributionDetails({num_files})" data-num-files="{num_files}" style="cursor: pointer;">
            <div class="dist-label">{num_files} file(s)</div>
            <div class="dist-bar-container">
                <div class="dist-bar" style="width: {bar_width}%"></div>
            </div>
            <div class="dist-value">{count} ({percentage:.1f}%)</div>
            <div class="dist-click-hint">📊 Click for details</div>
        </div>
        """
    
    file_dist_html += """
    </div>
    <script id="fileDistData" type="application/json">
    """ + json.dumps(file_dist_data) + """
    </script>
    """
    
    # Generate cluster size distribution HTML (if clustering data available)
    cluster_dist_html = ""
    cluster_dist_data = {}
    cluster_summary = {}
    
    if clusters_data:
        # Process clusters to build size distribution
        cluster_sizes = []
        component_clusters = {}
        
        for comp_key, comp_data in clusters_data.items():
            if isinstance(comp_data, dict) and 'clusters' in comp_data:
                component_id = comp_data.get('component_id', comp_key)
                method = comp_data.get('method', 'unknown')
                clusters = comp_data.get('clusters', {})
                
                component_clusters[component_id] = {
                    'method': method,
                    'clusters': {},
                    'total_clusters': len(clusters)
                }
                
                for cluster_key, cluster_data in clusters.items():
                    num_vertices = cluster_data.get('num_vertices', 0)
                    vertices = cluster_data.get('vertices', [])
                    
                    # Check for duplicates (same file appearing multiple times in cluster)
                    file_counts = {}
                    for vertex in vertices:
                        fname = vertex.get('filename', 'unknown')
                        file_counts[fname] = file_counts.get(fname, 0) + 1
                    
                    has_duplicates = any(count > 1 for count in file_counts.values())
                    
                    cluster_sizes.append(num_vertices)
                    component_clusters[component_id]['clusters'][cluster_key] = {
                        'size': num_vertices,
                        'has_duplicates': has_duplicates,
                        'file_counts': file_counts,
                        'vertices': vertices
                    }
        
        # Build size distribution
        size_dist = {}
        for size in cluster_sizes:
            size_dist[size] = size_dist.get(size, 0) + 1
        
        # Calculate summary statistics
        if cluster_sizes:
            cluster_summary = {
                'total_clusters': len(cluster_sizes),
                'avg_cluster_size': sum(cluster_sizes) / len(cluster_sizes),
                'min_cluster_size': min(cluster_sizes),
                'max_cluster_size': max(cluster_sizes),
                'median_cluster_size': sorted(cluster_sizes)[len(cluster_sizes)//2],
                'components_with_clustering': len(component_clusters)
            }
        
        # Generate cluster distribution HTML
        if size_dist:
            cluster_dist_html = """
    <div id="clusterDistChart" class="distribution-section">
    """
            
            total_clusters = sum(size_dist.values())
            
            for size in sorted(size_dist.keys()):
                count = size_dist[size]
                percentage = 100 * count / total_clusters
                bar_width = min(95, percentage)
                
                # Calculate statistics for clusters of this size
                clusters_with_duplicates = 0
                clusters_without_duplicates = 0
                total_files_in_clusters = 0
                components_with_this_size = set()
                
                for comp_id, comp_info in component_clusters.items():
                    for cluster_key, cluster_info in comp_info['clusters'].items():
                        if cluster_info['size'] == size:
                            components_with_this_size.add(comp_id)
                            if cluster_info['has_duplicates']:
                                clusters_with_duplicates += 1
                            else:
                                clusters_without_duplicates += 1
                            total_files_in_clusters += len(cluster_info['file_counts'])
                
                avg_files_per_cluster = total_files_in_clusters / count if count > 0 else 0
                
                cluster_dist_data[size] = {
                    'count': count,
                    'percentage': percentage,
                    'clusters_with_duplicates': clusters_with_duplicates,
                    'clusters_without_duplicates': clusters_without_duplicates,
                    'avg_files_per_cluster': avg_files_per_cluster,
                    'components_represented': len(components_with_this_size)
                }
                
                cluster_dist_html += f"""
        <div class="distribution-row clickable-bar" onclick="showClusterDistributionDetails({size})" data-cluster-size="{size}" style="cursor: pointer;">
            <div class="dist-label">{size} vertex/vertices</div>
            <div class="dist-bar-container">
                <div class="dist-bar" style="width: {bar_width}%"></div>
            </div>
            <div class="dist-value">{count} cluster{'s' if count != 1 else ''} ({percentage:.1f}%)</div>
            <div class="dist-click-hint">📊 Click for details</div>
        </div>
        """
            
            cluster_dist_html += """
    </div>
    <script id="clusterDistData" type="application/json">
    """ + json.dumps(cluster_dist_data) + """
    </script>
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
    
    # Generate parameters section HTML
    params_html = ""
    if parameters:
        # Parameter tooltips
        tooltips = {
            'mz_cutoff': 'Maximum m/z value difference allowed for pairing two features',
            'rt_cutoff': 'Maximum retention time difference (seconds) allowed for pairing',
            'edges_cutoff': 'Maximum Euclidean distance in scaled space for feature pairing',
            'best_match_only': 'Keep only the best scoring match for each feature pair',
            'optimize_pairing': 'Use optimized KD-tree indexing instead of basic matching',
            'rt_start_trim': 'Exclude features with RT values below this threshold (seconds)',
            'rt_end_trim': 'Exclude features with RT values above this threshold (seconds)',
            'postpair_normalize_coordinates': 'Rescale coordinate axes to balance ranges after pairing for better distance calculation',
        }
        
        params_html = """
        <div class="stats-section">
            <h3>Analysis Configuration</h3>
        """
        
        # Pairing parameters
        pairing_params = [
            ('mz_cutoff', 'M/Z Cutoff'),
            ('rt_cutoff', 'RT Cutoff'),
            ('edges_cutoff', 'Euclidean Cutoff'),
            ('best_match_only', 'Best Match Only'),
            ('optimize_pairing', 'Optimized Pairing'),
            ('rt_start_trim', 'RT Start Trim'),
            ('rt_end_trim', 'RT End Trim'),
            ('postpair_normalize_coordinates', 'Post-Pair Normalization'),
        ]
        
        if any(parameters.get(p[0]) for p in pairing_params):
            params_html += """
            <div class="params-subsection">
                <h4 class="params-subtitle">Pairing Configuration</h4>
                <div class="stats-grid">
            """
            for param_key, param_label in pairing_params:
                value = parameters.get(param_key)
                if value is not None:
                    display_value = str(value) if not isinstance(value, bool) else ('Yes' if value else 'No')
                    tooltip = tooltips.get(param_key, '')
                    tooltip_attr = f' title="{tooltip}"' if tooltip else ''
                    params_html += f"""
                    <div class="stat-box"{tooltip_attr}>
                        <div class="stat-label">{param_label}</div>
                        <div class="stat-value" style="font-size: 1.5em; margin-top: 8px;">{display_value}</div>
                    </div>
                    """
            params_html += """
                </div>
            </div>
            """
        
        # Graph parameter tooltips
        graph_tooltips = {
            'graph_edge_cutoff': 'Maximum edge distance for network graph edges',
            'graph_mz_cutoff': 'Filter graph edges by maximum m/z difference',
            'graph_rt_cutoff': 'Filter graph edges by maximum retention time difference',
            'graph_test_fraction': 'Use only a fraction of edges for testing/analysis (0-1)',
            'graph_layout_engine': 'Graphviz layout algorithm: sfdp (force-directed), dot (hierarchical), neato (spring), etc.',
        }
        
        # Graph parameters
        graph_params = [
            ('graph_edge_cutoff', 'Edge Cutoff'),
            ('graph_mz_cutoff', 'M/Z Cutoff'),
            ('graph_rt_cutoff', 'RT Cutoff'),
            ('graph_test_fraction', 'Test Fraction'),
            ('graph_layout_engine', 'Layout Engine'),
        ]
        
        if any(parameters.get(p[0]) for p in graph_params):
            params_html += """
            <div class="params-subsection">
                <h4 class="params-subtitle">Network Graph Configuration</h4>
                <div class="stats-grid">
            """
            for param_key, param_label in graph_params:
                value = parameters.get(param_key)
                if value is not None and (not isinstance(value, (int, float)) or value != 1):
                    display_value = str(value)
                    tooltip = graph_tooltips.get(param_key, '')
                    tooltip_attr = f' title="{tooltip}"' if tooltip else ''
                    params_html += f"""
                    <div class="stat-box"{tooltip_attr}>
                        <div class="stat-label">{param_label}</div>
                        <div class="stat-value" style="font-size: 1.5em; margin-top: 8px;">{display_value}</div>
                    </div>
                    """
            params_html += """
                </div>
            </div>
            """
        
        # Clustering parameter tooltips
        clustering_tooltips = {
            'enable_clustering': 'Enable community detection to identify clusters in the feature network',
            'clustering_method': 'Algorithm for community detection: edge_betweenness (slow, accurate), louvain, walktrap, label_propagation',
            'clustering_use_weights': 'Use edge distance weights in clustering algorithm for better community structure',
            'clustering_weight_mode': 'Weight calculation: inverse (1/(1+distance)) gives higher weight to closer features, distance uses raw distances',
        }
        
        # Clustering parameters
        clustering_params = [
            ('enable_clustering', 'Clustering Enabled'),
            ('clustering_method', 'Clustering Method'),
            ('clustering_use_weights', 'Use Weights'),
            ('clustering_weight_mode', 'Weight Mode'),
        ]
        
        if parameters.get('enable_clustering') or parameters.get('clustering_method'):
            params_html += """
            <div class="params-subsection">
                <h4 class="params-subtitle">Clustering Configuration</h4>
                <div class="stats-grid">
            """
            for param_key, param_label in clustering_params:
                value = parameters.get(param_key)
                if value is not None:
                    display_value = str(value) if not isinstance(value, bool) else ('Yes' if value else 'No')
                    tooltip = clustering_tooltips.get(param_key, '')
                    tooltip_attr = f' title="{tooltip}"' if tooltip else ''
                    params_html += f"""
                    <div class="stat-box"{tooltip_attr}>
                        <div class="stat-label">{param_label}</div>
                        <div class="stat-value" style="font-size: 1.5em; margin-top: 8px;">{display_value}</div>
                    </div>
                    """
            params_html += """
                </div>
            </div>
            """
        
        params_html += """
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
        
        .clickable-bar {{
            transition: all 0.3s ease;
            border-radius: 8px;
            padding: 10px !important;
            margin-bottom: 10px !important;
        }}
        
        .clickable-bar:hover {{
            background: rgba(102, 126, 234, 0.1);
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.2);
            transform: translateX(5px);
        }}
        
        .dist-click-hint {{
            font-size: 0.8em;
            color: #999;
            margin-left: auto;
            margin-right: 10px;
            opacity: 0;
            transition: opacity 0.3s ease;
        }}
        
        .clickable-bar:hover .dist-click-hint {{
            opacity: 1;
        }}
        
        /* Modal Styles */
        .modal {{
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0, 0, 0, 0.4);
            animation: fadeIn 0.3s ease;
        }}
        
        @keyframes fadeIn {{
            from {{ opacity: 0; }}
            to {{ opacity: 1; }}
        }}
        
        .modal-content {{
            background: white;
            margin: 5% auto;
            padding: 30px;
            border-radius: 12px;
            width: 90%;
            max-width: 600px;
            box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
            animation: slideDown 0.3s ease;
            max-height: 80vh;
            overflow-y: auto;
        }}
        
        @keyframes slideDown {{
            from {{
                transform: translateY(-50px);
                opacity: 0;
            }}
            to {{
                transform: translateY(0);
                opacity: 1;
            }}
        }}
        
        .modal-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
            padding-bottom: 15px;
            border-bottom: 2px solid #667eea;
        }}
        
        .modal-title {{
            font-size: 1.5em;
            font-weight: bold;
            color: #333;
        }}
        
        .close-btn {{
            font-size: 1.5em;
            cursor: pointer;
            color: #999;
            transition: color 0.3s ease;
        }}
        
        .close-btn:hover {{
            color: #333;
        }}
        
        .modal-detail-section {{
            margin-bottom: 25px;
        }}
        
        .modal-detail-title {{
            font-size: 1.1em;
            font-weight: 600;
            color: #667eea;
            margin-bottom: 12px;
            border-left: 4px solid #667eea;
            padding-left: 10px;
        }}
        
        .detail-stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
            margin-bottom: 15px;
        }}
        
        .detail-stat-box {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        
        .detail-stat-num {{
            font-size: 1.8em;
            font-weight: bold;
            color: #667eea;
        }}
        
        .detail-stat-label {{
            font-size: 0.85em;
            color: #555;
            margin-top: 5px;
        }}
        
        .composition-table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border: 1px solid #ddd;
            border-radius: 8px;
            overflow: hidden;
        }}
        
        .composition-table thead {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
        }}
        
        .composition-table th {{
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}
        
        .composition-table td {{
            padding: 12px;
            border-bottom: 1px solid #eee;
        }}
        
        .composition-table tbody tr:hover {{
            background: #f5f5f5;
        }}
        
        .composition-table tbody tr:last-child td {{
            border-bottom: none;
        }}
        
        .has-duplicates {{
            color: #e74c3c !important;
            background-color: #fadbd8 !important;
            font-weight: bold;
        }}
        
        .no-duplicates {{
            color: #27ae60 !important;
            background-color: #d5f4e6 !important;
            font-weight: bold;
        }}
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
        
        .file-appearance-row {{
            display: flex;
            align-items: center;
            margin-bottom: 15px;
            gap: 15px;
            padding: 10px;
            border-radius: 6px;
            background: #f9f9f9;
        }}
        
        .file-app-label {{
            min-width: 150px;
            font-weight: 500;
            color: #333;
            font-size: 0.95em;
            word-break: break-word;
        }}
        
        .file-app-bar-container {{
            flex: 1;
            height: 28px;
            background: #ddd;
            border-radius: 5px;
            overflow: hidden;
        }}
        
        .file-app-bar {{
            height: 100%;
            background: linear-gradient(90deg, #45B7D1 0%, #667eea 100%);
            transition: width 0.3s ease;
        }}
        
        .file-app-value {{
            min-width: 120px;
            text-align: right;
            font-weight: 500;
            color: #666;
            font-size: 0.9em;
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
            
            .file-appearance-row {{
                flex-direction: column;
                align-items: flex-start;
            }}
            
            .file-app-bar-container {{
                width: 100%;
            }}
            
            .file-app-label {{
                min-width: auto;
            }}
        }}
        
        /* Parameters section styling */
        .params-subsection {{
            margin-bottom: 30px;
        }}
        
        .params-subtitle {{
            color: #667eea;
            margin-bottom: 20px;
            font-size: 1.2em;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
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
            
            <!-- Cluster Size Distribution -->
            {f'''<div class="stats-section">
                <h3>Cluster Size Distribution</h3>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("total_clusters", 0)}</div>
                        <div class="stat-label">Total Clusters</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("avg_cluster_size", 0):.2f}</div>
                        <div class="stat-label">Avg Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("min_cluster_size", 0)}</div>
                        <div class="stat-label">Min Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("max_cluster_size", 0)}</div>
                        <div class="stat-label">Max Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("median_cluster_size", 0)}</div>
                        <div class="stat-label">Median Cluster Size</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{cluster_summary.get("components_with_clustering", 0)}</div>
                        <div class="stat-label">Components with Clustering</div>
                    </div>
                </div>
                <div class="distribution-section">
                    {cluster_dist_html}
                </div>
            </div>''' if cluster_summary else ''}
            
            <!-- Graph Composition Statistics -->
            {comp_html}
            
            <!-- Peptide Identification Statistics -->
            {ident_html}
            
            <!-- Analysis Configuration -->
            {params_html}
        </div>
        
        <div class="footer">
            Generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')} | Feature Pairing Analysis Report
        </div>
    </div>
    
    <!-- Modal for file distribution details -->
    <div id="fileDistModal" class="modal">
        <div class="modal-content">
            <div class="modal-header">
                <div class="modal-title" id="modalTitle">File Matching Details</div>
                <span class="close-btn" onclick="closeFileDistModal()">&times;</span>
            </div>
            <div id="modalBody"></div>
        </div>
    </div>
    
    <script>
        // Load file distribution data
        const fileDistData = JSON.parse(document.getElementById('fileDistData').textContent);
        
        function showFileDistributionDetails(numFiles) {{
            const data = fileDistData[numFiles];
            if (!data) return;
            
            // Build modal content
            let modalBody = document.getElementById('modalBody');
            
            let html = `
                <div class="modal-detail-section">
                    <div class="modal-detail-title">Overview</div>
                    <div class="detail-stats">
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.total}}</div>
                            <div class="detail-stat-label">Total Groups</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.percentage.toFixed(1)}}%</div>
                            <div class="detail-stat-label">Of All Groups</div>
                        </div>
                    </div>
                </div>
                
                <div class="modal-detail-section">
                    <div class="modal-detail-title">File Uniqueness & Multiplicity</div>
                    <div class="detail-stats">
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.all_unique}}</div>
                            <div class="detail-stat-label">All Files Unique</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.with_duplicates}}</div>
                            <div class="detail-stat-label">With Duplicates</div>
                        </div>
                    </div>
                    <p style="color: #666; font-size: 0.9em; margin-top: 10px;">
                        <strong>All Files Unique:</strong> Each of the ${{numFiles}} files appears exactly once<br>
                        <strong>With Duplicates:</strong> At least one file appears multiple times in the group
                    </p>
                </div>
            `;
            
            if (data.composition && data.composition.length > 0) {{
                html += `
                    <div class="modal-detail-section">
                        <div class="modal-detail-title">File Composition Breakdown</div>
                        <table class="composition-table">
                            <thead>
                                <tr>
                                    <th>Files Involved</th>
                                    <th style="text-align: right;">Count</th>
                                    <th style="text-align: right;">Percentage</th>
                                </tr>
                            </thead>
                            <tbody>
                `;
                
                data.composition.forEach(comp => {{
                    html += `
                        <tr>
                            <td>${{comp.composition}}</td>
                            <td style="text-align: right;"><strong>${{comp.count}}</strong></td>
                            <td style="text-align: right;">${{comp.percentage.toFixed(1)}}%</td>
                        </tr>
                    `;
                }});
                
                html += `
                            </tbody>
                        </table>
                    </div>
                `;
            }}
            
            if (data.multiplicity && data.multiplicity.length > 0) {{
                html += `
                    <div class="modal-detail-section">
                        <div class="modal-detail-title">Feature Multiplicity per File</div>
                        <p style="color: #666; font-size: 0.85em; margin-bottom: 12px;">
                            Shows how many features appear from each file: e.g., "2x2,1x1" means 2 files with 2 features each, 1 file with 1 feature
                        </p>
                        <table class="composition-table">
                            <thead>
                                <tr>
                                    <th>Multiplicity Pattern</th>
                                    <th style="text-align: right;">Count</th>
                                    <th style="text-align: right;">Percentage</th>
                                </tr>
                            </thead>
                            <tbody>
                `;
                
                data.multiplicity.forEach(mult => {{
                    html += `
                        <tr>
                            <td><strong>${{mult.pattern}}</strong></td>
                            <td style="text-align: right;"><strong>${{mult.count}}</strong></td>
                            <td style="text-align: right;">${{mult.percentage.toFixed(1)}}%</td>
                        </tr>
                    `;
                }});
                
                html += `
                            </tbody>
                        </table>
                    </div>
                `;
            }}
            
            if (data.file_appearance && data.file_appearance.length > 0) {{
                html += `
                    <div class="modal-detail-section">
                        <div class="modal-detail-title">File Appearance in Matches</div>
                        <p style="color: #666; font-size: 0.85em; margin-bottom: 15px;">
                            Shows which files appear in what percentage of the ${{data.total}} matching groups
                        </p>
                `;
                
                data.file_appearance.forEach(file => {{
                    const barWidth = Math.min(95, file.percentage);
                    html += `
                        <div class="file-appearance-row">
                            <div class="file-app-label">${{file.filename}}</div>
                            <div class="file-app-bar-container">
                                <div class="file-app-bar" style="width: ${{barWidth}}%"></div>
                            </div>
                            <div class="file-app-value">${{file.count}} (${{file.percentage.toFixed(1)}}%)</div>
                        </div>
                    `;
                }});
                
                html += `
                    </div>
                `;
            }}
            
            document.getElementById('modalTitle').textContent = `File Distribution Details: ${{numFiles}} File(s)`;
            modalBody.innerHTML = html;
            
            // Show modal
            document.getElementById('fileDistModal').style.display = 'block';
        }}
        
        function closeFileDistModal() {{
            document.getElementById('fileDistModal').style.display = 'none';
        }}
        
        // Close modal when clicking outside
        window.onclick = function(event) {{
            const modal = document.getElementById('fileDistModal');
            if (event.target === modal) {{
                modal.style.display = 'none';
            }}
        }}
    </script>
    
    <!-- Modal for cluster distribution details -->
    <div id="clusterDistModal" class="modal">
        <div class="modal-content">
            <div class="modal-header">
                <div class="modal-title" id="clusterModalTitle">Cluster Size Distribution Details</div>
                <span class="close-btn" onclick="closeClusterDistModal()">&times;</span>
            </div>
            <div id="clusterModalBody"></div>
        </div>
    </div>
    
    <script>
        // Load cluster distribution data (if available)
        let clusterDistData = {{}};
        const clusterDataElement = document.getElementById('clusterDistData');
        if (clusterDataElement) {{
            clusterDistData = JSON.parse(clusterDataElement.textContent);
        }}
        
        function showClusterDistributionDetails(clusterSize) {{
            const data = clusterDistData[clusterSize];
            if (!data) return;
            
            let modalBody = document.getElementById('clusterModalBody');
            
            const duplicatePercentage = data.count > 0 ? (100 * data.clusters_with_duplicates / data.count) : 0;
            const cleanPercentage = data.count > 0 ? (100 * data.clusters_without_duplicates / data.count) : 0;
            
            let html = `
                <div class="modal-detail-section">
                    <div class="modal-detail-title">Overview</div>
                    <div class="detail-stats">
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.count}}</div>
                            <div class="detail-stat-label">Clusters of This Size</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.percentage.toFixed(1)}}%</div>
                            <div class="detail-stat-label">Of All Clusters</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{clusterSize}}</div>
                            <div class="detail-stat-label">Feature(s) per Cluster</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.avg_files_per_cluster.toFixed(2)}}</div>
                            <div class="detail-stat-label">Avg Files per Cluster</div>
                        </div>
                    </div>
                </div>
                
                <div class="modal-detail-section">
                    <div class="modal-detail-title">Duplicate Detection Statistics</div>
                    <div class="detail-stats">
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.clusters_with_duplicates}}</div>
                            <div class="detail-stat-label">Clusters with Duplicates</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{duplicatePercentage.toFixed(1)}}%</div>
                            <div class="detail-stat-label">Percentage</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.clusters_without_duplicates}}</div>
                            <div class="detail-stat-label">Clean Clusters</div>
                        </div>
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{cleanPercentage.toFixed(1)}}%</div>
                            <div class="detail-stat-label">Percentage</div>
                        </div>
                    </div>
                    <p style="color: #666; font-size: 0.85em; margin-top: 10px;">
                        <strong>Clusters with Duplicates:</strong> Contain features from the same file appearing multiple times<br>
                        <strong>Clean Clusters:</strong> Each file appears at most once
                    </p>
                </div>
                
                <div class="modal-detail-section">
                    <div class="modal-detail-title">Component Distribution</div>
                    <div class="detail-stats">
                        <div class="detail-stat-box">
                            <div class="detail-stat-num">${{data.components_represented}}</div>
                            <div class="detail-stat-label">Components with This Cluster Size</div>
                        </div>
                    </div>
                    <p style="color: #666; font-size: 0.85em; margin-top: 10px;">
                        This cluster size is represented across ${{data.components_represented}} different components in the graph.
                    </p>
                </div>
            `;
            
            const pluralS = data.count !== 1 ? 's' : '';
            document.getElementById('clusterModalTitle').textContent = `Cluster Size ${{clusterSize}} (${{data.count}} cluster${{pluralS}})`;
            modalBody.innerHTML = html;
            
            // Show modal
            document.getElementById('clusterDistModal').style.display = 'block';
        }}
        
        function closeClusterDistModal() {{
            document.getElementById('clusterDistModal').style.display = 'none';
        }}
        
        // Close cluster modal when clicking outside
        window.onclick = function(event) {{
            const fileModal = document.getElementById('fileDistModal');
            const clusterModal = document.getElementById('clusterDistModal');
            if (event.target === fileModal) {{
                fileModal.style.display = 'none';
            }}
            if (event.target === clusterModal) {{
                clusterModal.style.display = 'none';
            }}
        }}
    </script>
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
    parser.add_argument("--clusters_json", default=None, help="Path to clusters JSON file (from build_network_graph --output-clusters-json)")
    
    # Pairing parameters
    parser.add_argument("--mz_cutoff", type=float, default=None, help="Maximum m/z difference for pairing")
    parser.add_argument("--rt_cutoff", type=float, default=None, help="Maximum RT difference for pairing")
    parser.add_argument("--edges_cutoff", type=float, default=None, help="Maximum euclidean distance cutoff for pairing")
    parser.add_argument("--best_match_only", action="store_true", help="Keep only best match for each feature pair")
    parser.add_argument("--optimize_pairing", action="store_true", help="Use optimized KD-tree based pairing")
    
    # Trimming parameters
    parser.add_argument("--rt_start_trim", type=float, default=None, help="Seconds to trim from start of retention time")
    parser.add_argument("--rt_end_trim", type=float, default=None, help="Maximum RT value to keep")
    
    # Network graph parameters
    parser.add_argument("--graph_edge_cutoff", type=float, default=None, help="Maximum euclidean distance for graph edges")
    parser.add_argument("--graph_mz_cutoff", type=float, default=None, help="Maximum m/z difference for graph")
    parser.add_argument("--graph_rt_cutoff", type=float, default=None, help="Maximum RT difference for graph")
    parser.add_argument("--graph_test_fraction", type=float, default=1, help="Fraction of edges to use for testing")
    parser.add_argument("--graph_layout_engine", default=None, help="Graphviz layout engine (dot, sfdp, fdp, neato, twopi, circo)")
    
    # Clustering parameters
    parser.add_argument("--enable_clustering", action="store_true", help="Enable community detection clustering")
    parser.add_argument("--clustering_method", default="edge_betweenness", help="Clustering method (louvain, walktrap, label_propagation, edge_betweenness)")
    parser.add_argument("--clustering_use_weights", action="store_true", help="Use edge distances as weights in clustering")
    parser.add_argument("--clustering_weight_mode", default="inverse", help="Weight mode (inverse or distance)")
    
    # Post-pairing parameters
    parser.add_argument("--postpair_normalize_coordinates", action="store_true", help="Rescale axes to balance ranges after pairing")
    
    # Input files tracking
    parser.add_argument("--input_files", nargs='*', default=[], help="List of input feature files used in pairing")
    
    args = parser.parse_args()
    
    # Build parameters dict
    parameters = {
        'mz_cutoff': args.mz_cutoff,
        'rt_cutoff': args.rt_cutoff,
        'edges_cutoff': args.edges_cutoff,
        'best_match_only': args.best_match_only,
        'optimize_pairing': args.optimize_pairing,
        'rt_start_trim': args.rt_start_trim,
        'rt_end_trim': args.rt_end_trim,
        'postpair_normalize_coordinates': args.postpair_normalize_coordinates,
        'graph_edge_cutoff': args.graph_edge_cutoff,
        'graph_mz_cutoff': args.graph_mz_cutoff,
        'graph_rt_cutoff': args.graph_rt_cutoff,
        'graph_test_fraction': args.graph_test_fraction,
        'graph_layout_engine': args.graph_layout_engine,
        'enable_clustering': args.enable_clustering,
        'clustering_method': args.clustering_method,
        'clustering_use_weights': args.clustering_use_weights,
        'clustering_weight_mode': args.clustering_weight_mode,
        'input_files': args.input_files,
    }
    
    print(f"\n{'='*70}")
    print(f"Generating Pairing Report")
    print(f"{'='*70}")
    
    generate_pairing_report(args.paired_json, args.output_summary, args.output_stats, args.edges_json, args.output_html, args.graph_composition_json, args.clusters_json, parameters)
    
    print(f"\n✓ Report saved to {os.path.basename(args.output_summary)}")
    print(f"✓ Statistics saved to {os.path.basename(args.output_stats)}")
    if args.output_html:
        print(f"✓ HTML report saved to {os.path.basename(args.output_html)}")


if __name__ == "__main__":
    main()

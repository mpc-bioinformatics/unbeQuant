#!/usr/bin/env python3
"""
2D Alignment Analysis - Scatter and Vector Field Visualization

Extensions to compare_retention_time_alignment.py that provide 2D visualization
of both m/z and RT corrections using scatter plots with encoding and vector fields.

Each matched pep_ident pair serves as an anchor point showing the displacement
in (m/z, RT) space.
"""

import argparse
import json
import glob
import os
from typing import Dict, List, Tuple
from collections import defaultdict
from datetime import datetime
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Import helper functions from main script
import sys
sys.path.insert(0, os.path.dirname(__file__))
from compare_retention_time_alignment import (
    load_feature_json, get_coordinate, resolve_mz_collisions,
    match_features_by_closest_mz_and_rt, round_with_precision
)


def compare_files_2d(file1_data: List[Dict], file2_data: List[Dict], 
                     file1_name: str, file2_name: str, one_to_one_only: bool = True,
                     mz_source: str = 'geo_center', rt_source: str = 'geo_center') -> Tuple[List[Dict], Dict]:
    """
    Compare two feature files and return BOTH m/z and RT differences (2D alignment data).
    
    Returns matched pairs with:
    - mz_file1, rt_file1: coordinates from file 1
    - mz_file2, rt_file2: coordinates from file 2
    - mz_correction: mz_file2 - mz_file1
    - rt_correction: rt_file2 - rt_file1
    - pep_ident: matched peptide identifier
    
    Args:
        file1_data, file2_data: Feature lists
        one_to_one_only: Only include 1:1 matches (pep_idents appearing exactly once in both files)
        mz_source, rt_source: Coordinate sources ('geo_center', 'center', 'start')
    
    Returns:
        Tuple of:
        - List of match dicts with 2D correction data
        - Statistics dict
    """
    # Group features by pep_ident
    file1_by_pep = defaultdict(list)
    file2_by_pep = defaultdict(list)
    
    for feature in file1_data:
        pep_idents = feature.get('l_raw_pep_ident', [])
        if pep_idents:
            if isinstance(pep_idents, list):
                for pep in pep_idents:
                    file1_by_pep[pep].append(feature)
            else:
                file1_by_pep[pep_idents].append(feature)
    
    for feature in file2_data:
        pep_idents = feature.get('l_raw_pep_ident', [])
        if pep_idents:
            if isinstance(pep_idents, list):
                for pep in pep_idents:
                    file2_by_pep[pep].append(feature)
            else:
                file2_by_pep[pep_idents].append(feature)
    
    # Find common pep_idents
    common_peps = set(file1_by_pep.keys()) & set(file2_by_pep.keys())
    
    # Filter to 1:1 matches if requested
    if one_to_one_only:
        common_peps = {pep for pep in common_peps 
                      if len(file1_by_pep[pep]) == 1 and len(file2_by_pep[pep]) == 1}
    
    if not common_peps:
        return [], {'num_matches': 0, 'mean_mz_diff': 0.0, 'mean_rt_diff': 0.0}
    
    # Resolve m/z collisions
    file1_first_per_pep = {pep: file1_by_pep[pep][0] for pep in common_peps}
    pep_to_mz = resolve_mz_collisions(file1_first_per_pep, mz_source=mz_source)
    
    # Collect all matched pairs
    all_matches = []
    mz_corrections = []
    rt_corrections = []
    
    for pep_ident in common_peps:
        f1_features = file1_by_pep[pep_ident]
        f2_features = file2_by_pep[pep_ident]
        
        # Match by closest point in (m/z, RT) space
        matches = match_features_by_closest_mz_and_rt(f1_features, f2_features,
                                                       mz_source=mz_source, rt_source=rt_source)
        
        for f1, f2 in matches:
            # Extract coordinates
            f1_mz = get_coordinate(f1, 'mz', mz_source)
            f1_rt = get_coordinate(f1, 'rt', rt_source)
            f2_mz = get_coordinate(f2, 'mz', mz_source)
            f2_rt = get_coordinate(f2, 'rt', rt_source)
            
            # Calculate corrections
            mz_diff = f2_mz - f1_mz
            rt_diff = f2_rt - f1_rt
            
            match_dict = {
                'mz_file1': f1_mz,
                'rt_file1': f1_rt,
                'mz_file2': f2_mz,
                'rt_file2': f2_rt,
                'mz_correction': mz_diff,
                'rt_correction': rt_diff,
                'pep_ident': pep_ident,
                'displacement_magnitude': np.sqrt(mz_diff**2 + rt_diff**2)
            }
            
            all_matches.append(match_dict)
            mz_corrections.append(mz_diff)
            rt_corrections.append(rt_diff)
    
    # Calculate statistics
    stats = {}
    if all_matches:
        stats['num_matches'] = len(all_matches)
        stats['mean_mz_diff'] = float(np.mean(mz_corrections))
        stats['std_mz_diff'] = float(np.std(mz_corrections))
        stats['mean_rt_diff'] = float(np.mean(rt_corrections))
        stats['std_rt_diff'] = float(np.std(rt_corrections))
        
        # Correlation between m/z and RT corrections
        if len(all_matches) > 2:
            corr = np.corrcoef(mz_corrections, rt_corrections)[0, 1]
            stats['mz_rt_correlation'] = float(corr) if not np.isnan(corr) else 0.0
        else:
            stats['mz_rt_correlation'] = 0.0
    else:
        stats = {
            'num_matches': 0,
            'mean_mz_diff': 0.0,
            'std_mz_diff': 0.0,
            'mean_rt_diff': 0.0,
            'std_rt_diff': 0.0,
            'mz_rt_correlation': 0.0
        }
    
    return all_matches, stats


def calculate_stats_from_matches(matches: List[Dict]) -> Dict:
    """
    Calculate statistics from a list of matches.
    
    Args:
        matches: List of match dictionaries with mz_correction, rt_correction
    
    Returns:
        Statistics dictionary
    """
    if not matches:
        return {
            'num_matches': 0,
            'mean_mz_diff': 0.0,
            'std_mz_diff': 0.0,
            'mean_rt_diff': 0.0,
            'std_rt_diff': 0.0,
            'mz_rt_correlation': 0.0
        }
    
    mz_corrections = [m['mz_correction'] for m in matches]
    rt_corrections = [m['rt_correction'] for m in matches]
    
    stats = {
        'num_matches': len(matches),
        'mean_mz_diff': float(np.mean(mz_corrections)),
        'std_mz_diff': float(np.std(mz_corrections)),
        'mean_rt_diff': float(np.mean(rt_corrections)),
        'std_rt_diff': float(np.std(rt_corrections)),
    }
    
    # Correlation between m/z and RT corrections
    if len(matches) > 2:
        corr = np.corrcoef(mz_corrections, rt_corrections)[0, 1]
        stats['mz_rt_correlation'] = float(corr) if not np.isnan(corr) else 0.0
    else:
        stats['mz_rt_correlation'] = 0.0
    
    return stats


def filter_matches_by_percentile(matches: List[Dict], percentile: float = 99.0, 
                                 metric: str = 'displacement_magnitude') -> List[Dict]:
    """
    Filter matches to keep only those within the specified percentile.
    
    Args:
        matches: List of match dictionaries
        percentile: Percentile threshold (0-100). Default 99 keeps all but top 1% outliers
        metric: Which metric to use for filtering ('displacement_magnitude', 'mz_correction', 'rt_correction')
    
    Returns:
        Filtered list of matches within the percentile threshold
    """
    if not matches or len(matches) < 2:
        return matches
    
    # Extract metric values
    metric_values = [abs(m.get(metric, m.get('displacement_magnitude', 0))) for m in matches]
    
    # Calculate percentile threshold
    threshold = np.percentile(metric_values, percentile)
    
    # Filter matches
    filtered = [m for m, val in zip(matches, metric_values) if val <= threshold]
    
    return filtered


def create_scatter_plot_2d(file1_name: str, comparisons: Dict[str, List[Dict]]) -> Tuple[go.Figure, Dict[str, str]]:
    """
    Create scatter plot showing 2D alignment with encoding.
    
    X-axis: M/Z from file1
    Y-axis: RT from file1
    Point size: Magnitude of m/z correction |Δm/z|
    Point color: RT correction amount (Δrt), with colorscale
    
    Args:
        file1_name: Reference file name
        comparisons: {file2_name: [match_dicts]}
    
    Returns:
        Tuple of (Plotly Figure, {file2_name: description_text})
    """
    fig = go.Figure()
    descriptions = {}
    
    for file2_name, matches in sorted(comparisons.items()):
        if not matches:
            continue
        
        # Extract data
        mz_values = [m['mz_file1'] for m in matches]
        rt_values = [m['rt_file1'] for m in matches]
        mz_corr_mag = np.abs([m['mz_correction'] for m in matches])
        rt_corr = [m['rt_correction'] for m in matches]
        pep_idents = [m['pep_ident'] for m in matches]
        
        # Create hover text
        hover_text = [
            f"<b>{file2_name}</b><br>" +
            f"M/Z: {mz:.4f}<br>" +
            f"RT: {rt:.2f}s<br>" +
            f"ΔM/Z: {dm:.4f}<br>" +
            f"ΔRT: {dr:.2f}s<br>" +
            f"Magnitude: {mag:.4f}<br>" +
            f"Pep ID: {pep}"
            for mz, rt, dm, dr, mag, pep in zip(
                mz_values, rt_values,
                [m['mz_correction'] for m in matches],
                rt_corr,
                [m['displacement_magnitude'] for m in matches],
                pep_idents
            )
        ]
        
        # Add scatter trace
        fig.add_trace(go.Scatter(
            x=mz_values,
            y=rt_values,
            mode='markers',
            name=f'vs {file2_name}',
            marker=dict(
                size=mz_corr_mag * 20 + 3,  # Scale for visibility, min size 3
                color=rt_corr,
                colorscale='RdBu',
                showscale=(file2_name == list(comparisons.keys())[0]),  # Show colorbar for first trace
                colorbar=dict(
                    title='ΔRT (sec)',
                    thickness=15,
                    len=0.7,
                    x=1.02
                ),
                line=dict(width=0.5, color='white'),
                opacity=0.7
            ),
            hovertext=hover_text,
            hoverinfo='text'
        ))
        
        # Add description
        num_matches = len(matches)
        mean_mz = np.mean(mz_corr_mag)
        mean_rt = np.mean(np.abs(rt_corr))
        descriptions[file2_name] = (
            f"{num_matches} matches | "
            f"Mean |ΔM/Z|: {mean_mz:.4f} | "
            f"Mean |ΔRT|: {mean_rt:.2f}s"
        )
    
    # Update layout
    fig.update_layout(
        title={
            'text': f'2D Alignment Scatter – {file1_name}<br><sub>Point size = |ΔM/Z|, Color = ΔRT</sub>',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18}
        },
        xaxis_title='M/Z (m/z)',
        yaxis_title='Retention Time (seconds)',
        hovermode='closest',
        height=700,
        width=1000,
        template='plotly_white',
        font=dict(size=12),
        legend=dict(
            title='File Comparisons',
            yanchor='top',
            y=0.99,
            xanchor='left',
            x=0.01
        )
    )
    
    return fig, descriptions


def draw_arrow_with_head(fig, x_start, y_start, x_end, y_end, color, hovertext, 
                         file2_name, arrow_size=0.04):
    """
    Draw an arrow line with a properly rotated arrowhead.
    
    The arrowhead points in the direction of (x_end - x_start, y_end - y_start).
    
    Args:
        fig: Plotly figure
        x_start, y_start: Starting point
        x_end, y_end: Ending point
        color: Arrow color
        hovertext: Hover text
        file2_name: File name for legend
        arrow_size: Size of arrowhead relative to arrow length
    """
    # Draw arrow line
    fig.add_trace(go.Scatter(
        x=[x_start, x_end],
        y=[y_start, y_end],
        mode='lines',
        line=dict(width=3.5, color=color),
        hoverinfo='skip',
        showlegend=False,
        name=file2_name
    ))
    
    # Calculate angle of arrow
    dx = x_end - x_start
    dy = y_end - y_start
    angle = np.arctan2(dy, dx)
    
    # Arrowhead as a small triangle at the endpoint
    # Two lines forming a V shape pointing in the direction of the arrow
    arrow_angle = 20 * np.pi / 180  # 20 degree angle for arrowhead sides
    
    # Left side of arrowhead
    x_left = x_end - arrow_size * np.cos(angle - arrow_angle)
    y_left = y_end - arrow_size * np.sin(angle - arrow_angle)
    
    # Right side of arrowhead
    x_right = x_end - arrow_size * np.cos(angle + arrow_angle)
    y_right = y_end - arrow_size * np.sin(angle + arrow_angle)
    
    # Draw left side of arrowhead
    fig.add_trace(go.Scatter(
        x=[x_left, x_end],
        y=[y_left, y_end],
        mode='lines',
        line=dict(width=3.5, color=color),
        hoverinfo='skip',
        showlegend=False,
        name=file2_name
    ))
    
    # Draw right side of arrowhead
    fig.add_trace(go.Scatter(
        x=[x_right, x_end],
        y=[y_right, y_end],
        mode='lines',
        line=dict(width=3.5, color=color),
        hovertext=hovertext,
        hoverinfo='text',
        showlegend=False,
        name=file2_name
    ))


def create_vector_field_plot(file1_name: str, comparisons: Dict[str, List[Dict]], 
                            cell_size_mz: float = 100.0, cell_size_rt: float = 100.0) -> Tuple[go.Figure, Dict[str, str]]:
    """
    Create vector field plot showing aggregated displacements in (m/z, RT) space.
    
    Divides the space into cells of fixed size (default 100×100) and shows one averaged arrow per cell.
    - Each cell is 100 m/z units × 100 RT seconds
    - X direction: Average m/z displacement (Δm/z)
    - Y direction: Average RT displacement (Δrt)
    - Arrow direction: Points in direction of average displacement
    - Arrow color: File comparison color
    
    Args:
        file1_name: Reference file name
        comparisons: {file2_name: [match_dicts]}
        cell_size_mz: M/Z range per cell (default 100.0)
        cell_size_rt: RT range per cell in seconds (default 100.0)
    
    Returns:
        Tuple of (Plotly Figure, {file2_name: description_text})
    """
    fig = go.Figure()
    descriptions = {}
    
    # Color palette for different file comparisons
    colors_palette = [
        '#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8',
        '#F7DC6F', '#BB8FCE', '#85C1E2', '#F8B195', '#C7CEEA'
    ]
    
    # Collect all m/z and RT values to determine grid bounds
    all_mz = []
    all_rt = []
    for matches in comparisons.values():
        all_mz.extend([m['mz_file1'] for m in matches])
        all_rt.extend([m['rt_file1'] for m in matches])
    
    if not all_mz or not all_rt:
        return fig, descriptions
    
    # Determine grid bounds and align to cell boundaries
    mz_min, mz_max = np.min(all_mz), np.max(all_mz)
    rt_min, rt_max = np.min(all_rt), np.max(all_rt)
    
    # Align to cell boundaries (round down/up)
    mz_min_aligned = np.floor(mz_min / cell_size_mz) * cell_size_mz
    mz_max_aligned = np.ceil(mz_max / cell_size_mz) * cell_size_mz
    rt_min_aligned = np.floor(rt_min / cell_size_rt) * cell_size_rt
    rt_max_aligned = np.ceil(rt_max / cell_size_rt) * cell_size_rt
    
    mz_range = mz_max_aligned - mz_min_aligned
    rt_range = rt_max_aligned - rt_min_aligned
    
    # Calculate per-grid-cell statistics
    for file_idx, (file2_name, matches) in enumerate(sorted(comparisons.items())):
        if not matches:
            continue
        
        # Get color for this file
        color = colors_palette[file_idx % len(colors_palette)]
        
        # Create grid cells and aggregate data
        grid_data = {}  # (cell_i, cell_j) -> {'mz_center': ..., 'rt_center': ..., 'displacements': [...]}
        
        for match in matches:
            mz = match['mz_file1']
            rt = match['rt_file1']
            dmz = match['mz_correction']
            drt = match['rt_correction']
            
            # Determine grid cell (which 100x100 block does this point belong to)
            cell_i = int((mz - mz_min_aligned) / cell_size_mz)
            cell_j = int((rt - rt_min_aligned) / cell_size_rt)
            
            # Clamp to grid bounds
            cell_i = np.clip(cell_i, 0, int(mz_range / cell_size_mz))
            cell_j = np.clip(cell_j, 0, int(rt_range / cell_size_rt))
            
            cell_key = (cell_i, cell_j)
            
            if cell_key not in grid_data:
                # Calculate grid cell center
                mz_center = mz_min_aligned + (cell_i + 0.5) * cell_size_mz
                rt_center = rt_min_aligned + (cell_j + 0.5) * cell_size_rt
                grid_data[cell_key] = {
                    'mz_center': mz_center,
                    'rt_center': rt_center,
                    'displacements': [],
                    'pep_idents': []
                }
            
            grid_data[cell_key]['displacements'].append((dmz, drt))
            grid_data[cell_key]['pep_idents'].append(match['pep_ident'])
        
        # Find max magnitude for auto-scaling
        max_mag = 0
        for cell_data in grid_data.values():
            if cell_data['displacements']:
                avg_dmz = np.mean([d[0] for d in cell_data['displacements']])
                avg_drt = np.mean([d[1] for d in cell_data['displacements']])
                mag = np.sqrt(avg_dmz**2 + avg_drt**2)
                max_mag = max(max_mag, mag)
        
        # Scale factor for arrows (make them span ~50% of cell for better visibility)
        scale_mz = 0.5 * cell_size_mz / (max_mag + 1e-6)
        scale_rt = 0.5 * cell_size_rt / (max_mag + 1e-6)
        
        # Draw arrows for each grid cell
        for (cell_i, cell_j), cell_data in sorted(grid_data.items()):
            mz_center = cell_data['mz_center']
            rt_center = cell_data['rt_center']
            displacements = cell_data['displacements']
            pep_idents = cell_data['pep_idents']
            
            # Calculate average displacement
            avg_dmz = np.mean([d[0] for d in displacements])
            avg_drt = np.mean([d[1] for d in displacements])
            avg_mag = np.sqrt(avg_dmz**2 + avg_drt**2)
            num_points = len(displacements)
            
            # Scaled endpoint
            mz_end = mz_center + avg_dmz * scale_mz
            rt_end = rt_center + avg_drt * scale_rt
            
            hovertext = (
                f"<b>{file2_name}</b><br>" +
                f"Cell center: ({mz_center:.1f}, {rt_center:.1f})<br>" +
                f"Avg ΔM/Z: {avg_dmz:.4f}<br>" +
                f"Avg ΔRT: {avg_drt:.2f}s<br>" +
                f"Avg magnitude: {avg_mag:.4f}<br>" +
                f"Points in cell: {num_points}"
            )
            
            # Draw arrow with properly rotated arrowhead
            draw_arrow_with_head(fig, mz_center, rt_center, mz_end, rt_end, 
                               color, hovertext, file2_name)
        
        # Add one legend entry for this file comparison
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='lines',
            line=dict(width=2, color=color),
            name=f'{file2_name}',
            hoverinfo='skip'
        ))
        
        # Statistics
        num_cells = len(grid_data)
        total_points = len(matches)
        avg_points_per_cell = total_points / num_cells if num_cells > 0 else 0
        grid_dims = (int(mz_range / cell_size_mz), int(rt_range / cell_size_rt))
        descriptions[file2_name] = (
            f"{num_cells} grid cells ({grid_dims[0]}×{grid_dims[1]}) of {cell_size_mz:.0f}×{cell_size_rt:.0f} | "
            f"{total_points} points | "
            f"Avg {avg_points_per_cell:.1f} points/cell"
        )
    
    
    # Update layout
    fig.update_layout(
        title={
            'text': f'2D Alignment Vector Field – {file1_name}<br><sub>Arrows show feature displacement in (m/z, RT) space</sub>',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18}
        },
        xaxis_title='M/Z (m/z)',
        yaxis_title='Retention Time (seconds)',
        hovermode='closest',
        height=700,
        width=1000,
        template='plotly_white',
        font=dict(size=12),
        legend=dict(
            title='File Comparisons',
            yanchor='top',
            y=0.99,
            xanchor='left',
            x=0.01
        )
    )
    
    return fig, descriptions


def generate_2d_alignment_report(file_data: Dict[str, List[Dict]], 
                                 file_mapping: Dict[str, str],
                                 output_html: str,
                                 mz_source: str = 'geo_center',
                                 rt_source: str = 'geo_center',
                                 one_to_one_only: bool = False,
                                 percentile_threshold: float = 99.0):
    """
    Generate standalone HTML report with 2D alignment visualizations.
    
    Includes:
    - Scatter plot with encoding (size=|ΔM/Z|, color=ΔRT)
    - Vector field (arrows showing displacements)
    - Statistics table
    
    Args:
        file_data: {filename: [features]}
        file_mapping: {basename: real_filename}
        output_html: Output HTML file path
        mz_source, rt_source: Coordinate sources
        one_to_one_only: Use only 1:1 matches
        percentile_threshold: Filter to this percentile to exclude outliers (default 99.0)
    """
    print(f"\n{'='*70}")
    print(f"Generating 2D Alignment Report...")
    print(f"{'='*70}\n")
    
    # Collect 2D alignment data
    all_comparisons_2d = {}
    all_stats_2d = {}
    
    for file1_name in sorted(file_data.keys()):
        file1_data = file_data[file1_name]
        comparisons = {}
        stats_dict = {}
        
        for file2_name in sorted(file_data.keys()):
            if file1_name == file2_name:
                continue
            
            file2_data = file_data[file2_name]
            print(f"Analyzing {file1_name} vs {file2_name}...", end='', flush=True)
            
            matches, stats = compare_files_2d(
                file1_data, file2_data, file1_name, file2_name,
                one_to_one_only=one_to_one_only,
                mz_source=mz_source,
                rt_source=rt_source
            )
            
            # Filter to percentile threshold to remove outliers
            original_count = len(matches) if matches else 0
            if matches:
                matches = filter_matches_by_percentile(matches, percentile=percentile_threshold)
                # Recalculate stats after filtering
                if matches:
                    stats = calculate_stats_from_matches(matches)
            filtered_count = len(matches) if matches else 0
            
            if matches:
                comparisons[file2_name] = matches
                stats_dict[file2_name] = stats
                if filtered_count < original_count:
                    print(f" ✓ ({filtered_count} of {original_count} matches after {percentile_threshold}% filtering)")
                else:
                    print(f" ✓ ({len(matches)} matches)")
            else:
                print(f" ⚠ (no matches)")
        
        if comparisons:
            all_comparisons_2d[file1_name] = comparisons
            all_stats_2d[file1_name] = stats_dict
    
    if not all_comparisons_2d:
        print("✗ No valid comparisons for 2D analysis")
        return
    
    # Create plots
    print(f"\nGenerating visualizations...\n")
    figures = []
    
    for file1_name in sorted(all_comparisons_2d.keys()):
        comparisons = all_comparisons_2d[file1_name]
        
        # Scatter plot
        scatter_fig, scatter_desc = create_scatter_plot_2d(file1_name, comparisons)
        figures.append(('scatter', file1_name, scatter_fig, scatter_desc))
        
        # Vector field plot
        vector_fig, vector_desc = create_vector_field_plot(file1_name, comparisons, 
                                                          cell_size_mz=100.0, cell_size_rt=100.0)
        figures.append(('vector', file1_name, vector_fig, vector_desc))
    
    # Generate HTML
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>2D Alignment Report</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            * {{
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                background-color: #f5f5f5;
                color: #333;
                line-height: 1.6;
            }}
            .header {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 40px 20px;
                text-align: center;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            .header h1 {{
                font-size: 2.5em;
                margin-bottom: 10px;
            }}
            .header p {{
                font-size: 1.1em;
                opacity: 0.9;
            }}
            .container {{
                max-width: 1400px;
                margin: 0 auto;
                padding: 30px 20px;
            }}
            .plot-section {{
                background: white;
                border-radius: 8px;
                padding: 20px;
                margin-bottom: 30px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            }}
            .plot-container {{
                width: 100%;
                height: 750px;
                margin-bottom: 15px;
            }}
            .plot-description {{
                background: #F0F8FF;
                border-left: 4px solid #667eea;
                padding: 12px;
                border-radius: 4px;
                font-size: 0.95em;
                color: #333;
            }}
            .stats {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 15px;
                margin-bottom: 20px;
            }}
            .stat-box {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 15px;
                border-radius: 6px;
                text-align: center;
            }}
            .stat-box .value {{
                font-size: 2em;
                font-weight: bold;
            }}
            .stat-box .label {{
                font-size: 0.9em;
                opacity: 0.9;
                margin-top: 5px;
            }}
            .info-box {{
                background: #FFFACD;
                border-left: 4px solid #FFD700;
                padding: 15px;
                margin-bottom: 20px;
                border-radius: 4px;
            }}
            .info-box p {{
                margin: 5px 0;
                font-size: 0.95em;
            }}
            .stats-table {{
                width: 100%;
                border-collapse: collapse;
                background: white;
                border-radius: 6px;
                overflow: hidden;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin-bottom: 30px;
            }}
            .stats-table th {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 15px;
                text-align: left;
                font-weight: 600;
            }}
            .stats-table td {{
                padding: 12px 15px;
                border-bottom: 1px solid #eee;
            }}
            .stats-table tr:hover {{
                background-color: #f9f9f9;
            }}
            .stats-table tr:last-child td {{
                border-bottom: none;
            }}
            .footer {{
                text-align: center;
                padding: 20px;
                color: #666;
                font-size: 0.9em;
            }}
            .section-title {{
                font-size: 1.5em;
                color: #667eea;
                margin-top: 30px;
                margin-bottom: 15px;
                padding-bottom: 10px;
                border-bottom: 2px solid #667eea;
            }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>2D Alignment Analysis Report</h1>
            <p>Scatter and Vector Field Visualization of Feature Displacements</p>
        </div>
        
        <div class="container">
            <div class="info-box">
                <p><strong>Analysis Description:</strong></p>
                <p><strong>Scatter Plot:</strong> Each point represents a matched feature pair. 
                Point size indicates the magnitude of m/z correction (|ΔM/Z|), and color shows the RT correction (ΔRT). 
                This reveals if m/z and RT corrections are correlated.</p>
                <p style="margin-top: 10px;"><strong>Vector Field:</strong> Each arrow shows the displacement 
                of a feature in (m/z, RT) space from file1 to file2. Arrow direction shows both m/z and RT 
                displacement simultaneously. This reveals the 2D "flow" pattern of instrument drift.</p>
            </div>
            
            <div class="info-box" style="background: #F0F8FF; border-left-color: #667eea;">
                <p><strong>Configuration:</strong></p>
                <p>Coordinate sources: M/Z = <code>{mz_source}</code>, RT = <code>{rt_source}</code></p>
                {f'<p>Filtering: Only 1:1 matches (pep_idents appearing exactly once) used</p>' if one_to_one_only else '<p>Filtering: All matches included</p>'}
            </div>
    """
    
    # Add summary stats
    total_files = len(file_data)
    total_comparisons = sum(len(v) for v in all_comparisons_2d.values())
    total_matches = sum(
        sum(len(matches) for matches in comps.values()) 
        for comps in all_comparisons_2d.values()
    )
    
    html_content += f"""
            <div class="stats">
                <div class="stat-box">
                    <div class="value">{total_files}</div>
                    <div class="label">Input Files</div>
                </div>
                <div class="stat-box">
                    <div class="value">{total_comparisons}</div>
                    <div class="label">Comparisons</div>
                </div>
                <div class="stat-box">
                    <div class="value">{total_matches}</div>
                    <div class="label">Total Matches</div>
                </div>
            </div>
    """
    
    # Stats table
    html_content += """
            <div style="background: white; border-radius: 8px; padding: 20px; margin-bottom: 30px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);">
                <h2 style="margin-bottom: 15px; color: #667eea;">2D Alignment Statistics</h2>
                <table class="stats-table">
                    <thead>
                        <tr>
                            <th>File 1</th>
                            <th>File 2</th>
                            <th>Matches</th>
                            <th>Mean ΔM/Z</th>
                            <th>Mean ΔRT (sec)</th>
                            <th>M/Z-RT Corr</th>
                        </tr>
                    </thead>
                    <tbody>
    """
    
    for file1_name in sorted(all_stats_2d.keys()):
        stats_dict = all_stats_2d[file1_name]
        mapped_file1 = file_mapping.get(file1_name, file1_name)
        
        for file2_name in sorted(stats_dict.keys()):
            stats = stats_dict[file2_name]
            mapped_file2 = file_mapping.get(file2_name, file2_name)
            
            num_matches = stats['num_matches']
            mean_mz = stats['mean_mz_diff']
            mean_rt = stats['mean_rt_diff']
            corr = stats['mz_rt_correlation']
            
            html_content += f"""
                        <tr>
                            <td><strong>{mapped_file1}</strong></td>
                            <td><strong>{mapped_file2}</strong></td>
                            <td style="text-align: center;">{num_matches}</td>
                            <td style="text-align: right;">{mean_mz:.6f}</td>
                            <td style="text-align: right;">{mean_rt:.2f}</td>
                            <td style="text-align: right;">{corr:.3f}</td>
                        </tr>
            """
    
    html_content += """
                    </tbody>
                </table>
            </div>
    """
    
    # Add plots
    current_file = None
    for plot_type, file1_name, fig, descriptions in figures:
        mapped_file = file_mapping.get(file1_name, file1_name)
        
        # Add section title for new file
        if current_file != file1_name:
            html_content += f'<h2 class="section-title">Reference File: {mapped_file}</h2>'
            current_file = file1_name
        
        # Add plot
        plot_html = fig.to_html(include_plotlyjs=False, div_id=f"plot_{plot_type}_{file1_name}")
        
        # Build descriptions
        desc_html = ""
        if descriptions:
            desc_html += '<div style="margin-top: 15px;">'
            for file2_name, desc_text in sorted(descriptions.items()):
                mapped_file2 = file_mapping.get(file2_name, file2_name)
                desc_html += f'<div class="plot-description"><strong>vs {mapped_file2}:</strong> {desc_text}</div>'
            desc_html += '</div>'
        
        plot_title = 'Scatter Plot (size=|ΔM/Z|, color=ΔRT)' if plot_type == 'scatter' else 'Vector Field'
        
        html_content += f"""
            <div class="plot-section">
                <h3 style="color: #667eea; margin-bottom: 10px;">{plot_title} – {mapped_file}</h3>
                <div class="plot-container" id="plot_{plot_type}_{file1_name}">
                    {plot_html}
                </div>
                {desc_html}
            </div>
        """
    
    html_content += """
        </div>
        
        <div class="footer">
            <p>Generated by 2D Alignment Analysis Tool</p>
        </div>
    </body>
    </html>
    """
    
    # Write output
    output_dir = os.path.dirname(output_html)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    with open(output_html, 'w') as f:
        f.write(html_content)
    
    print(f"\n✓ 2D Alignment Report generated: {output_html}")
    print(f"\nSummary:")
    print(f"  Files analyzed: {total_files}")
    print(f"  File comparisons: {total_comparisons}")
    print(f"  Total matches: {total_matches}")
    print(f"  Plots generated: {len(figures)} ({len(figures)//2} per reference file)")


def main():
    parser = argparse.ArgumentParser(
        description='2D Alignment Analysis - Scatter and Vector Field Visualization'
    )
    parser.add_argument('--input_dir', default='results/feature_analysis/feature_data_lists',
                       help='Directory containing feature_data.json files')
    parser.add_argument('--input_pattern', default='*_feature_data.json',
                       help='Glob pattern for feature data files')
    parser.add_argument('--output_html', default='results/2d_alignment_report.html',
                       help='Output HTML report file (default: results/2d_alignment_report.html)')
    parser.add_argument('--mz_source', default='geo_center', choices=['geo_center', 'center', 'start'],
                       help='Which m/z value to use: geo_center (default), center, or start')
    parser.add_argument('--rt_source', default='geo_center', choices=['geo_center', 'center', 'start'],
                       help='Which RT value to use: geo_center (default), center, or start')
    parser.add_argument('--allow_multi', action='store_true',
                       help='Allow pep_idents that appear multiple times (default: uses only 1:1 matches)')
    parser.add_argument('--percentile', type=float, default=99.0,
                       help='Percentile threshold to filter outliers (0-100, default 99.0). '
                            'Use 99 to exclude top 1%% of displacements, 100 to keep all.')
    
    args = parser.parse_args()
    
    # Determine filtering mode: default is to use only 1:1 pep_ident matches
    use_one_to_one = not args.allow_multi
    
    # Find and load files
    pattern = os.path.join(args.input_dir, args.input_pattern)
    json_files = sorted(glob.glob(pattern))
    
    if not json_files:
        print(f"✗ No files matching pattern: {pattern}")
        return
    
    print(f"\n{'='*70}")
    print(f"2D Alignment Analysis")
    print(f"{'='*70}")
    print(f"Found {len(json_files)} feature data files\n")
    
    # Load files
    file_data = {}
    file_mapping = {}
    
    for filepath in json_files:
        basename = os.path.basename(filepath).replace('_feature_data.json', '')
        print(f"Loading {basename}...", end='', flush=True)
        
        data = load_feature_json(filepath)
        if data:
            # Extract real filename
            if data and isinstance(data, list) and len(data) > 0 and 'filename' in data[0]:
                real_filename = data[0]['filename'].replace('_feature_data.json', '')
            else:
                real_filename = basename
            
            file_data[real_filename] = data
            file_mapping[basename] = real_filename
            print(f" ✓ ({len(data)} features)")
        else:
            print(f" ✗")
    
    if len(file_data) < 2:
        print(f"\n✗ Need at least 2 files for comparison")
        return
    
    print(f"\nFile mapping: {file_mapping}\n")
    
    # Generate report (using one_to_one setting from args)
    generate_2d_alignment_report(
        file_data, file_mapping, args.output_html,
        mz_source=args.mz_source,
        rt_source=args.rt_source,
        one_to_one_only=use_one_to_one,
        percentile_threshold=args.percentile
    )


if __name__ == '__main__':
    main()

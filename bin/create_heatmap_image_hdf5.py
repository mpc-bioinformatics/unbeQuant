#!/usr/bin/env python3
"""
Create heatmap image from processed mzML spectral data.
Refactored from map_mzml_batch_tsv.py
Optimized for fast processing with vectorized operations and batched HDF5 writes.
"""

import argparse
import pickle
import json
import numpy as np
import os
import gc
import h5py
from PIL import Image
from typing import Dict, List, Tuple

#TODO: Create an 3D-Heatmap
def initialize_rgb_dict(scale_colors: bool = True):
    """Precompute RGB dictionary once for all files."""
    if scale_colors:
        print("Precomputing RGB intensity dictionary...")
        
        def intensity_to_rgb(i):
            r = g = b = 0
            for bit_pos in range(24):
                if i & (1 << bit_pos):
                    channel = bit_pos % 3
                    bit_index = bit_pos // 3
                    
                    if channel == 0:
                        b |= (1 << bit_index)
                    elif channel == 1:
                        g |= (1 << bit_index)
                    else:
                        r |= (1 << bit_index)
            return (r, g, b)
        
        # Vectorized RGB conversion
        intensity_array = np.arange(16777216, dtype=np.uint32)
        rgb_array = np.zeros((16777216, 3), dtype=np.uint8)
        
        for bit_pos in range(24):
            mask = (intensity_array >> bit_pos) & 1
            channel = bit_pos % 3
            bit_index = bit_pos // 3
            rgb_array[:, channel] |= (mask << bit_index).astype(np.uint8)
        
        rgb_dict = {i: tuple(rgb_array[i]) for i in range(0, 16777216, 1)}
        print("RGB dictionary precomputed and ready")
        return rgb_array  # Return array for faster lookup
    return None


def intensity_to_rgb_lookup(intensity_int: int, rgb_array: np.ndarray, scale_colors: bool = True):
    """Look up RGB values using precomputed array."""
    if scale_colors and rgb_array is not None:
        return tuple(rgb_array[intensity_int])
    return rgb_int2tuple(intensity_int)


def intensity_to_rgb_vectorized(intensities: np.ndarray, rgb_array: np.ndarray, scale_colors: bool = True):
    """Vectorized RGB lookup for arrays of intensities."""
    if scale_colors and rgb_array is not None:
        return rgb_array[intensities.astype(np.uint32)]
    # Fallback for non-scaled colors
    return np.stack([
        (intensities // 256 // 256) % 256,
        (intensities // 256) % 256,
        intensities % 256
    ], axis=-1).astype(np.uint8)


def rgb_int2tuple(rgbint: int):
    """Convert RGB integer to tuple for non-scaled colors."""
    return (rgbint // 256 // 256 % 256, rgbint // 256 % 256, rgbint % 256)


def save_heatmap_to_hdf5(img_frame: np.ndarray, output_path: str):
    """Save heatmap image data to HDF5 format for efficient storage and retrieval."""
    print(f"  Saving heatmap to HDF5: {os.path.basename(output_path)}")
    with h5py.File(output_path, 'w') as hf:
        # Store RGB image data with compression
        hf.create_dataset('img_frame', data=img_frame, compression='gzip', compression_opts=4, dtype=np.uint8)
        # Store metadata
        hf.attrs['format'] = 'heatmap_image_v1'
        hf.attrs['img_height'] = img_frame.shape[0]
        hf.attrs['img_width'] = img_frame.shape[1]


def load_heatmap_from_hdf5(hdf5_path: str):
    """Load heatmap image data from HDF5 format."""
    with h5py.File(hdf5_path, 'r') as hf:
        img_frame = hf['img_frame'][:]
    return img_frame


def export_heatmap_to_png(hdf5_path: str, png_output_path: str):
    """Export heatmap from HDF5 to PNG format."""
    img_frame = load_heatmap_from_hdf5(hdf5_path)
    print(f"  Exporting PNG from HDF5: {os.path.basename(png_output_path)}")
    img_pil = Image.fromarray(img_frame, mode='RGB')
    img_pil.save(png_output_path, format='PNG')
    img_pil = None
    gc.collect()


def create_feature_visualization(spec_total: np.ndarray, RTINSECONDS_arr: List[float],
                                mz_total_arr: List[float], mz_dict: Dict, rt_dict: Dict,
                                features_data: List[Dict], output_png_path: str):
    """
    Create a feature visualization image with boxes and center markers on a blank background.
    This is optimized for features-only mode to show feature locations without intensity data.
    """
    print(f"  Creating feature visualization image...")
    
    img_height = int(len(RTINSECONDS_arr) + 10)
    mz_total_num = (max(mz_total_arr) - min(mz_total_arr)) * (10 ** 2)  # Assuming round_up_to=2
    img_width = int(mz_total_num) + 10
    
    # Create blank white image
    img_frame = Image.new('RGB', (img_width, img_height), color='white')
    draw = __import__('PIL.ImageDraw', fromlist=['ImageDraw']).ImageDraw(img_frame)
    
    print(f"  Drawing {len(features_data)} feature boxes on {img_width}x{img_height} canvas...")
    
    # Draw each feature as a box with center marker
    for feat_idx, feature in enumerate(features_data):
        try:
            # Get feature bounds
            mz_start = feature.get('mz_start')
            mz_end = feature.get('mz_end')
            rt_start = feature.get('rt_start')
            rt_end = feature.get('rt_end')
            x_center = feature.get('x_center')
            y_center = feature.get('y_center')
            charge = feature.get('charge', 0)
            
            if None in [mz_start, mz_end, rt_start, rt_end, x_center, y_center]:
                continue
            
            # Map to pixel coordinates using dictionaries
            # Round to match the dictionary keys
            mz_start_rounded = round(mz_start, 2)
            mz_end_rounded = round(mz_end, 2)
            rt_start_rounded = round(rt_start, 8)
            rt_end_rounded = round(rt_end, 8)
            
            # Get pixel indices from dictionaries
            mz_start_idx = mz_dict.get(mz_start_rounded)
            mz_end_idx = mz_dict.get(mz_end_rounded)
            rt_start_idx = rt_dict.get(rt_start_rounded)
            rt_end_idx = rt_dict.get(rt_end_rounded)
            
            # If not in dict, estimate from nearest values
            if mz_start_idx is None:
                mz_start_idx = int((mz_start - min(mz_total_arr)) * 100)
            if mz_end_idx is None:
                mz_end_idx = int((mz_end - min(mz_total_arr)) * 100)
            if rt_start_idx is None:
                rt_start_idx = int(rt_start)
            if rt_end_idx is None:
                rt_end_idx = int(rt_end)
            
            # Ensure indices are within bounds
            mz_start_idx = max(0, min(img_width - 1, int(mz_start_idx)))
            mz_end_idx = max(0, min(img_width - 1, int(mz_end_idx)))
            rt_start_idx = max(0, min(img_height - 1, int(rt_start_idx)))
            rt_end_idx = max(0, min(img_height - 1, int(rt_end_idx)))
            
            # Ensure start < end
            if mz_start_idx > mz_end_idx:
                mz_start_idx, mz_end_idx = mz_end_idx, mz_start_idx
            if rt_start_idx > rt_end_idx:
                rt_start_idx, rt_end_idx = rt_end_idx, rt_start_idx
            
            # Choose color based on charge state
            charge_colors = {
                1: (255, 0, 0),      # Red
                2: (0, 0, 255),      # Blue
                3: (0, 128, 0),      # Green
                4: (255, 165, 0),    # Orange
            }
            box_color = charge_colors.get(charge, (128, 128, 128))  # Gray for unknown
            
            # Draw box
            draw.rectangle(
                [(mz_start_idx, rt_start_idx), (mz_end_idx, rt_end_idx)],
                outline=box_color,
                width=1
            )
            
            # Draw center marker (cross and circle)
            center_x = int((mz_start_idx + mz_end_idx) / 2)
            center_y = int((rt_start_idx + rt_end_idx) / 2)
            marker_size = 4
            
            # Draw cross
            draw.line([(center_x - marker_size, center_y), (center_x + marker_size, center_y)], 
                     fill=box_color, width=1)
            draw.line([(center_x, center_y - marker_size), (center_x, center_y + marker_size)], 
                     fill=box_color, width=1)
            
            # Draw small circle around center
            draw.ellipse(
                [(center_x - marker_size - 1, center_y - marker_size - 1),
                 (center_x + marker_size + 1, center_y + marker_size + 1)],
                outline=box_color,
                width=1
            )
            
            if (feat_idx + 1) % 1000 == 0:
                print(f"    Drawn {feat_idx + 1}/{len(features_data)} features", end='\r')
        
        except Exception as e:
            print(f"  Warning: Failed to draw feature {feat_idx}: {e}")
            continue
    
    print(f"\n  Successfully drawn all {len(features_data)} features")
    
    # Save to PNG
    print(f"  Saving feature visualization to PNG: {os.path.basename(output_png_path)}")
    img_frame.save(output_png_path, format='PNG')
    img_frame = None
    gc.collect()
    
    return output_png_path


def generate_heatmap_report(png_path: str, hdf5_path: str, basename: str, 
                           spec_total: np.ndarray, log_scale: bool = True,
                           scale_colors: bool = True, invert_colors: bool = False,
                           features_only: bool = False, feature_list: List[Dict] = None):
    """
    Generate an HTML report with an interactive Plotly heatmap at the center point.
    """
    report_path = png_path.replace('.png', '_report.html')
    
    # Calculate statistics
    num_features = len(spec_total)
    rt_min = spec_total[:, 0].min()
    rt_max = spec_total[:, 0].max()
    mz_min = spec_total[:, 1].min()
    mz_max = spec_total[:, 1].max()
    intensity_min = spec_total[:, 2].min()
    intensity_max = spec_total[:, 2].max()
    intensity_mean = spec_total[:, 2].mean()
    
    # Get file sizes
    png_size = os.path.getsize(png_path) / (1024 * 1024) if os.path.exists(png_path) else 0
    hdf5_size = os.path.getsize(hdf5_path) / (1024 * 1024) if os.path.exists(hdf5_path) else 0
    
    # Create 2D heatmap grid from spec_total data
    # Binning RT and m/z into a 2D grid
    rt_bins = int(np.ceil((rt_max - rt_min) / 5)) + 1  # 5-second bins for RT
    mz_bins = int(np.ceil((mz_max - mz_min) * 10)) + 1  # 0.1 m/z bins
    
    heatmap_2d = np.zeros((rt_bins, mz_bins))
    
    # Bin the data into 2D grid
    for rt, mz, intensity in spec_total:
        rt_idx = min(int((rt - rt_min) / 5), rt_bins - 1)
        mz_idx = min(int((mz - mz_min) * 10), mz_bins - 1)
        heatmap_2d[rt_idx, mz_idx] += intensity
    
    # Apply log scale if requested
    if log_scale and np.any(heatmap_2d > 0):
        heatmap_2d = np.log10(heatmap_2d + 1)
    
    # Generate Plotly data in JSON format
    import json
    heatmap_json = json.dumps(heatmap_2d.tolist())
    
    rt_range = np.arange(rt_min, rt_max + 5, 5)[:rt_bins]
    mz_range = np.arange(mz_min, mz_max + 0.1, 0.1)[:mz_bins]
    
    # Initialize feature visualization variables
    shapes_json = None
    feature_centers_json = None
    
    # Create feature boxes visualization if feature_list is provided
    if features_only and feature_list:
        # Prepare feature boxes as rectangles only (no markers to avoid overlapping/lag)
        shapes = []
        hover_x = []
        hover_y = []
        hover_text = []
        
        # Charge to color mapping
        charge_colors = {
            1: 'rgb(255, 0, 0)',      # Red
            2: 'rgb(0, 0, 255)',      # Blue
            3: 'rgb(0, 128, 0)',      # Green
            4: 'rgb(255, 165, 0)',    # Orange
        }
        
        for idx, feature in enumerate(feature_list):
            # Get feature boundaries
            mz_start = feature.get('mz_start', 0)
            mz_end = feature.get('mz_end', 0)
            rt_start = feature.get('rt_start', 0)
            rt_end = feature.get('rt_end', 0)
            
            # Skip if boundaries are invalid
            if mz_start >= mz_end or rt_start >= rt_end:
                continue
            
            # Get center coordinates for hover point
            if 'x_center' in feature and 'y_center' in feature:
                mz_center = feature['x_center']
                rt_center = feature['y_center']
            else:
                mz_center = (mz_start + mz_end) / 2
                rt_center = (rt_start + rt_end) / 2
            
            # Get charge state
            charge = feature.get('charge', 0)
            charge_color = charge_colors.get(charge, 'rgb(128, 128, 128)')
            
            # Create rectangle shape for feature box
            shape = dict(
                type="rect",
                x0=mz_start,
                y0=rt_start,
                x1=mz_end,
                y1=rt_end,
                line=dict(color=charge_color, width=1),
                fillcolor=charge_color,
                opacity=0.15,
                layer="below"
            )
            shapes.append(shape)
            
            # Store hover point data (center of each feature)
            hover_x.append(mz_center)
            hover_y.append(rt_center)
            
            # Hover info
            pep_ident = feature.get('pep_ident', 'N/A')
            hover_text.append(f"<b>Feature {idx}</b><br>m/z: {mz_center:.4f}<br>RT: {rt_center:.2f}s<br>Charge: +{charge}<br>m/z range: {mz_start:.4f}-{mz_end:.4f}<br>RT range: {rt_start:.2f}-{rt_end:.2f}s<br>ID: {pep_ident}")
        
        # Prepare JSON for shapes
        shapes_json = json.dumps(shapes)
        
        # Prepare hover data (invisible markers for text on hover)
        if hover_x:
            hover_marker_data = {
                'x': hover_x,
                'y': hover_y,
                'text': hover_text
            }
            feature_centers_json = json.dumps(hover_marker_data)
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Heatmap Report - {basename}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
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
            max-width: 1400px;
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
        
        .heatmap-section {{
            margin-bottom: 40px;
            text-align: center;
        }}
        
        .heatmap-section h2 {{
            color: #333;
            margin-bottom: 20px;
            font-size: 1.8em;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
        }}
        
        .heatmap-container {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            margin-bottom: 15px;
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
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .stat-box {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 25px;
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
            font-size: 0.95em;
            color: #555;
            font-weight: 500;
        }}
        
        .stat-unit {{
            font-size: 0.8em;
            color: #999;
            margin-top: 5px;
        }}
        
        .parameters-section {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        
        .parameters-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
        }}
        
        .param-item {{
            background: rgba(255,255,255,0.1);
            padding: 15px;
            border-radius: 8px;
            backdrop-filter: blur(10px);
        }}
        
        .param-label {{
            font-size: 0.9em;
            opacity: 0.9;
            margin-bottom: 5px;
        }}
        
        .param-value {{
            font-size: 1.1em;
            font-weight: 600;
        }}
        
        .footer {{
            background: #f8f9fa;
            padding: 20px 40px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
            border-top: 1px solid #ddd;
        }}
        
        .feature-mode-badge {{
            display: inline-block;
            background: #ff6b6b;
            color: white;
            padding: 5px 15px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
            margin-top: 10px;
        }}
        
        .feature-mode-badge.normal {{
            background: #51cf66;
        }}
        
        .feature-legend {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 10px;
            margin-top: 20px;
            display: inline-block;
        }}
        
        .feature-legend h4 {{
            color: #333;
            margin-bottom: 15px;
            font-size: 1.1em;
        }}
        
        .legend-items {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
        }}
        
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        
        .legend-color {{
            width: 30px;
            height: 30px;
            border: 2px solid;
            border-radius: 4px;
        }}
        
        .legend-label {{
            font-size: 0.9em;
            color: #555;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Heatmap Visualization Report</h1>
            <p>Spectral Data Intensity Map</p>
            <p style="font-size: 0.9em; margin-top: 10px;">File: {basename}</p>
        </div>
        
        <div class="content">
            <!-- Heatmap Display Section -->
            <div class="heatmap-section">
                <h2>Spectral Heatmap{' - Feature Visualization' if features_only else ''}</h2>
                <div class="heatmap-container" id="heatmapPlotly" style="width:100%; height:600px;"></div>
                {'<div class="feature-mode-badge">Features Only (No Background)</div>' if features_only else '<div class="feature-mode-badge normal">Full Heatmap with Background</div>'}
                
                {'<div class="feature-legend"><h4>Feature Visualization Legend</h4><div class="legend-items"><div class="legend-item"><div class="legend-color" style="border-color: #ff0000; background: rgba(255,0,0,0.1)"></div><div class="legend-label">Charge +1</div></div><div class="legend-item"><div class="legend-color" style="border-color: #0000ff; background: rgba(0,0,255,0.1)"></div><div class="legend-label">Charge +2</div></div><div class="legend-item"><div class="legend-color" style="border-color: #008000; background: rgba(0,128,0,0.1)"></div><div class="legend-label">Charge +3</div></div><div class="legend-item"><div class="legend-color" style="border-color: #ffa500; background: rgba(255,165,0,0.1)"></div><div class="legend-label">Charge +4</div></div><div class="legend-item"><div class="legend-color" style="border-color: #808080; background: rgba(128,128,128,0.1)"></div><div class="legend-label">Unknown Charge</div></div><div class="legend-item"><div class="legend-color" style="border-color: #667eea; background: rgba(102,126,234,0.1)"></div><div class="legend-label">Feature Center</div></div></div><p style="margin-top: 10px; font-size: 0.85em; color: #666;"><strong>Boxes:</strong> Feature boundaries (m/z × RT) | <strong>Markers:</strong> Feature center points with charge-state coloring</p></div>' if features_only else ''}
            </div>
            
            {f'<!-- Feature Boxes Display Section --><div class="heatmap-section"><h2>Feature Visualization</h2><div class="heatmap-container" id="featureBoxesPlotly" style="width:100%; height:600px;"></div><div id="featureInfoTooltip" class="feature-info-tooltip" style="display:none; position: absolute; background: white; border: 2px solid #667eea; border-radius: 8px; padding: 15px; max-width: 300px; box-shadow: 0 4px 15px rgba(0,0,0,0.2); z-index: 1000; font-size: 0.9em;"></div><p style="font-size: 0.9em; color: #666; margin-top: 15px;"><strong>Click on features</strong> to display detailed information. Box color indicates charge state: <span style="color: #ff0000;">●</span> +1, <span style="color: #0000ff;">●</span> +2, <span style="color: #008000;">●</span> +3, <span style="color: #ffa500;">●</span> +4</p></div>' if (features_only and feature_centers_json) else ''}
            
            <!-- Feature Statistics -->
            <div class="stats-section">
                <h3>Feature Statistics</h3>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{num_features:,}</div>
                        <div class="stat-label">Total Features</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{rt_min:.2f}</div>
                        <div class="stat-label">Min RT</div>
                        <div class="stat-unit">seconds</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{rt_max:.2f}</div>
                        <div class="stat-label">Max RT</div>
                        <div class="stat-unit">seconds</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{rt_max - rt_min:.2f}</div>
                        <div class="stat-label">RT Range</div>
                        <div class="stat-unit">seconds</div>
                    </div>
                </div>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{mz_min:.4f}</div>
                        <div class="stat-label">Min m/z</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{mz_max:.4f}</div>
                        <div class="stat-label">Max m/z</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{mz_max - mz_min:.4f}</div>
                        <div class="stat-label">m/z Range</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{(mz_max - mz_min):.4f}</div>
                        <div class="stat-label">Coverage</div>
                    </div>
                </div>
            </div>
            
            <!-- Intensity Statistics -->
            <div class="stats-section">
                <h3>Intensity Statistics</h3>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{intensity_min:,.0f}</div>
                        <div class="stat-label">Min Intensity</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{intensity_max:,.0f}</div>
                        <div class="stat-label">Max Intensity</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{intensity_mean:,.0f}</div>
                        <div class="stat-label">Mean Intensity</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{intensity_max / intensity_min if intensity_min > 0 else 0:.2f}x</div>
                        <div class="stat-label">Dynamic Range</div>
                    </div>
                </div>
            </div>
            
            <!-- Processing Parameters -->
            <div class="parameters-section">
                <h3 style="margin-bottom: 20px;">Processing Parameters</h3>
                <div class="parameters-grid">
                    <div class="param-item">
                        <div class="param-label">Log Scale</div>
                        <div class="param-value">{'✓ Enabled' if log_scale else '✗ Disabled'}</div>
                    </div>
                    <div class="param-item">
                        <div class="param-label">Scaled Colors</div>
                        <div class="param-value">{'✓ Enabled' if scale_colors else '✗ Disabled'}</div>
                    </div>
                    <div class="param-item">
                        <div class="param-label">Inverted Colors</div>
                        <div class="param-value">{'✓ Enabled' if invert_colors else '✗ Disabled'}</div>
                    </div>
                    <div class="param-item">
                        <div class="param-label">Features Only</div>
                        <div class="param-value">{'✓ Enabled' if features_only else '✗ Disabled'}</div>
                    </div>
                    <div class="param-item">
                        <div class="param-label">PNG File Size</div>
                        <div class="param-value">{png_size:.2f} MB</div>
                    </div>
                    <div class="param-item">
                        <div class="param-label">HDF5 File Size</div>
                        <div class="param-value">{hdf5_size:.2f} MB</div>
                    </div>
                </div>
            </div>
        </div>
        
        <script>
            // Prepare heatmap data
            const heatmapData = {heatmap_json};
            const rtRange = {json.dumps(rt_range.tolist())};
            const mzRange = {json.dumps(mz_range.tolist())};
            
            // Create Plotly heatmap
            const trace = {{
                z: heatmapData,
                x: mzRange,
                y: rtRange,
                type: 'heatmap',
                colorscale: 'Viridis',
                colorbar: {{
                    title: 'Intensity{"(log10)" if log_scale else ""}',
                    thickness: 20,
                    len: 0.7
                }},
                hovertemplate: '<b>m/z:</b> %{{x:.4f}}<br><b>RT:</b> %{{y:.2f}}s<br><b>Intensity:</b> %{{z:.2e}}<extra></extra>'
            }};
            
            const layout = {{
                title: '{basename} - Spectral Heatmap',
                xaxis: {{
                    title: 'm/z (mass-to-charge ratio)',
                    showgrid: true,
                    zeroline: false
                }},
                yaxis: {{
                    title: 'RT (Retention Time, seconds)',
                    showgrid: true,
                    zeroline: false
                }},
                width: null,
                height: 600,
                margin: {{ l: 80, r: 100, t: 80, b: 80 }},
                hovermode: 'closest',
                autosize: true
            }};
            
            const config = {{
                responsive: true,
                displayModeBar: true,
                displaylogo: false,
                modeBarButtonsToRemove: ['lasso2d', 'select2d']
            }};
            
            Plotly.newPlot('heatmapPlotly', [trace], layout, config);
            
            // Handle window resize
            window.addEventListener('resize', function() {{
                Plotly.Plots.resize('heatmapPlotly');
            }});
            
            {'// Render feature boxes visualization with click interaction' if feature_centers_json else ''}
            {f'''
            const featureShapes = {shapes_json};
            const clickMarkerData = {feature_centers_json};
            
            // Create invisible marker layer for click information
            const featureClickTrace = {{
                x: clickMarkerData.x,
                y: clickMarkerData.y,
                mode: 'markers',
                type: 'scatter',
                marker: {{
                    size: 0,
                    opacity: 0
                }},
                customdata: clickMarkerData.text,
                hovertemplate: '',
                textposition: 'none',
                name: 'Features',
                showlegend: false
            }};
            
            const featureLayout = {{
                title: '{basename} - Feature Boxes (Click for Details)',
                xaxis: {{
                    title: 'm/z (mass-to-charge ratio)',
                    showgrid: true,
                    zeroline: false
                }},
                yaxis: {{
                    title: 'RT (Retention Time, seconds)',
                    showgrid: true,
                    zeroline: false,
                    autorange: 'reversed'
                }},
                width: null,
                height: 600,
                margin: {{ l: 80, r: 100, t: 80, b: 80 }},
                hovermode: false,
                autosize: true,
                shapes: featureShapes,
                plot_bgcolor: 'rgba(240, 240, 245, 0.5)',
                paper_bgcolor: 'white'
            }};
            
            Plotly.newPlot('featureBoxesPlotly', [featureClickTrace], featureLayout, config);
            
            // Handle click events on features
            document.getElementById('featureBoxesPlotly').on('plotly_click', function(eventData) {{
                const point = eventData.points[0];
                const infoDiv = document.getElementById('featureInfoTooltip');
                let btnHtml = '<button onclick="';
                btnHtml += 'document.getElementById(' + String.fromCharCode(39) + 'featureInfoTooltip' + String.fromCharCode(39) + ').style.display=' + String.fromCharCode(39) + 'none' + String.fromCharCode(39);
                btnHtml += ';" style="float: right; background: none; border: none; font-size: 1.2em; cursor: pointer; padding: 0; margin: 0;">×</button>';
                const contentDiv = '<div style="clear: both; margin-top: 8px;">' + point.customdata + '</div>';
                infoDiv.innerHTML = btnHtml + contentDiv;
                infoDiv.style.display = 'block';
            }});
            
            window.addEventListener('resize', function() {{
                Plotly.Plots.resize('featureBoxesPlotly');
            }});
            ''' if feature_centers_json else ''}
        </script>
        
        <div class="footer">
            <p>Generated by unbeQuant - Heatmap Image Generator v1.0</p>
            <p style="margin-top: 10px; font-size: 0.85em;">This report visualizes the spectral intensity distribution across retention time and mass-to-charge ratio dimensions.</p>
        </div>
    </div>
</body>
</html>"""
    
    # Write HTML report
    print(f"  Generating HTML report: {os.path.basename(report_path)}")
    with open(report_path, 'w') as f:
        f.write(html_content)
    
    return report_path


def create_heatmap_image(spec_total: np.ndarray, RTINSECONDS_arr: List[float], 
                         mz_total_arr: List[float], mz_dict: Dict, rt_dict: Dict,
                         output_path: str, log_scale: bool = True, scale_colors: bool = True,
                         invert_colors: bool = False, save_png: bool = True, batch_size: int = 100000,
                         row_batch_size: int = 10, compression_level: int = None, features_only: bool = False,
                         save_report: bool = True, feature_data_json_path: str = None):
    """
    Create and save heatmap image directly to HDF5 datasets for minimal memory usage.
    Writes data in batches to HDF5 for optimal performance.
    Args:
        batch_size: Deprecated (kept for backwards compatibility)
        row_batch_size: Number of rows to prepare and write together (default: 10)
                       Higher values = faster writes but more RAM (each row ~3.2MB)
                       Lower values = slower writes but less RAM
        compression_level: GZIP compression level (None=no compression, 1=light, 4=default, 9=max)
                         None=fastest (1min, 36MB), 1=balanced (7min, 0.36MB), 4=original (13min, 0.2MB)
        save_png: If True, also exports a PNG copy from the HDF5 file
        features_only: If True, saves only the mapped feature points without the heatmap background (creates much smaller file)
        save_report: If True, generates an HTML report with the heatmap and statistics
    Returns: HDF5 file path
    """
    print(f"  Creating heatmap image...")
    
    rgb_array = initialize_rgb_dict(scale_colors)
    
    img_height = int(len(RTINSECONDS_arr) + 10)
    mz_total_num = (max(mz_total_arr) - min(mz_total_arr)) * (10 ** 2)  # Assuming round_up_to=2
    img_width = int(mz_total_num) + 10
    
    # Prepare HDF5 output path
    hdf5_output = output_path if output_path.endswith('.h5') else output_path.replace('.png', '.h5')
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(hdf5_output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"  Created output directory: {output_dir}")
    
    # Remove existing file if it exists to avoid lock issues
    if os.path.exists(hdf5_output):
        try:
            os.remove(hdf5_output)
            print(f"  Removed existing file: {hdf5_output}")
        except Exception as e:
            print(f"  Warning: Could not remove existing file: {e}")
    
    # If features_only mode, save mapped feature points and create visualization
    if features_only:
        print(f"  Features-only mode: saving mapped feature points and creating visualization...")
        try:
            import json
            with h5py.File(hdf5_output, 'w', libver='latest') as hf:
                print(f"  Mapping {len(spec_total)} spectrum points...")
                
                # Vectorize dictionary lookups
                try:
                    index_mz_arr = np.array([mz_dict.get(round(mz, 2), None) for mz in spec_total[:, 1]])
                    index_rt_arr = np.array([rt_dict.get(round(rt, 8), None) for rt in spec_total[:, 0]])
                except:
                    index_mz_arr = np.full(len(spec_total), None, dtype=object)
                    index_rt_arr = np.full(len(spec_total), None, dtype=object)
                    for i in range(len(spec_total)):
                        index_mz_arr[i] = mz_dict.get(round(spec_total[i][1], 2), None)
                        index_rt_arr[i] = rt_dict.get(round(spec_total[i][0], 8), None)
                
                # Filter valid entries
                valid_mask = (index_mz_arr != None) & (index_rt_arr != None)
                valid_indices = np.where(valid_mask)[0]
                print(f"  Found {len(valid_indices)} valid feature points to save")
                
                if len(valid_indices) > 0:
                    # Extract feature data
                    features_data = np.column_stack([
                        index_rt_arr[valid_indices].astype(int),
                        index_mz_arr[valid_indices].astype(int),
                        spec_total[valid_indices, 0],  # RT in seconds
                        spec_total[valid_indices, 1],  # m/z
                        spec_total[valid_indices, 2]   # intensity
                    ])
                    
                    # Create dataset for features
                    feature_dset = hf.create_dataset(
                        'features',
                        data=features_data,
                        compression='gzip',
                        compression_opts=4,
                        dtype=np.float64
                    )
                    feature_dset.attrs['description'] = 'Mapped feature points: [rt_idx, mz_idx, rt_seconds, mz, intensity]'
                    
                    # Store metadata
                    hf.attrs['format'] = 'features_only_v1'
                    hf.attrs['num_features'] = len(valid_indices)
                    hf.attrs['img_height'] = img_height
                    hf.attrs['img_width'] = img_width
                    hf.attrs['log_scale'] = log_scale
                    hf.attrs['scale_colors'] = scale_colors
                    hf.attrs['invert_colors'] = invert_colors
                    
                    print(f"  Successfully saved {len(valid_indices)} feature points to HDF5")
                else:
                    print(f"  Warning: No valid feature points found")
            
            # Create visualization from feature data if provided
            if save_png:
                if feature_data_json_path and os.path.exists(feature_data_json_path):
                    print(f"  Loading feature data from JSON: {os.path.basename(feature_data_json_path)}")
                    with open(feature_data_json_path, 'r') as f:
                        feature_list = json.load(f)
                    
                    # Create feature visualization
                    png_output = output_path
                    create_feature_visualization(spec_total, RTINSECONDS_arr, mz_total_arr, 
                                                mz_dict, rt_dict, feature_list, png_output)
                    
                    # Generate HTML report
                    if save_report:
                        basename = os.path.splitext(os.path.basename(png_output))[0]
                        generate_heatmap_report(png_output, hdf5_output, basename, spec_total, 
                                               log_scale, scale_colors, invert_colors, features_only=True,
                                               feature_list=feature_list)
                else:
                    print(f"  Note: Feature data JSON not provided, skipping visualization")
            
            return hdf5_output
        except Exception as e:
            print(f"  Error creating features-only HDF5 file: {e}")
            if os.path.exists(hdf5_output):
                try:
                    os.remove(hdf5_output)
                except:
                    pass
            raise
    
    # Standard mode: create full heatmap image
    
    # Create HDF5 file and datasets directly (no large arrays in memory)
    print(f"  Initializing HDF5 file: {os.path.basename(hdf5_output)}")
    try:
        with h5py.File(hdf5_output, 'w', libver='latest', swmr=False) as hf:
            # Create datasets with optional compression
            print(f"  Creating HDF5 datasets with shape ({img_height}, {img_width})...")
            if compression_level is not None:
                print(f"  Using GZIP compression level {compression_level}")
                img_frame_dset = hf.create_dataset(
                    'img_frame', 
                    shape=(img_height, img_width, 3), 
                    dtype=np.uint8, 
                    compression='gzip',
                    compression_opts=compression_level,
                    chunks=True,
                    fillvalue=255
                )
            else:
                print(f"  No compression (fastest write speed)")
                img_frame_dset = hf.create_dataset(
                    'img_frame', 
                    shape=(img_height, img_width, 3), 
                    dtype=np.uint8, 
                    chunks=True,
                    fillvalue=255
                )
            print(f"  Created 'img_frame' dataset")
            # Store metadata
            hf.attrs['format'] = 'heatmap_image_v1'
            hf.attrs['img_height'] = img_height
            hf.attrs['img_width'] = img_width
            print(f"  HDF5 file initialized and ready for data writing")
            # Normalize intensity values
            min_spec = spec_total[:, 2].min()
            max_spec = spec_total[:, 2].max()
            
            if log_scale:
                log_min = np.log1p(min_spec)
                log_max = np.log1p(max_spec)
                log_diff = log_max - log_min
            else:
                log_min = min_spec
                log_max = max_spec
                log_diff = max_spec - min_spec
            print(f"  Normalized intensity range: min={log_min}, max={log_max}")
            # Vectorized normalization: map spectrum points
            print(f"  Mapping {len(spec_total)} spectrum points to HDF5 datasets...")
            
            # Vectorize dictionary lookups
            try:
                index_mz_arr = np.array([mz_dict.get(round(mz, 2), None) for mz in spec_total[:, 1]])
                index_rt_arr = np.array([rt_dict.get(round(rt, 8), None) for rt in spec_total[:, 0]])
            except:
                # Fallback if dictionaries have issues
                index_mz_arr = np.full(len(spec_total), None, dtype=object)
                index_rt_arr = np.full(len(spec_total), None, dtype=object)
                for i in range(len(spec_total)):
                    index_mz_arr[i] = mz_dict.get(round(spec_total[i][1], 2), None)
                    index_rt_arr[i] = rt_dict.get(round(spec_total[i][0], 8), None)
            print(f"  Completed mapping spectrum points to indices, filtering valid entries...")
            # Filter valid entries
            valid_mask = (index_mz_arr != None) & (index_rt_arr != None)
            valid_indices = np.where(valid_mask)[0]
            print(f"  Found {len(valid_indices)} valid spectrum points to process")
            if len(valid_indices) == 0:
                print(f"  Warning: No valid spectrum points found")
                return hdf5_output
            
            # Pre-extract indices for chunked processing
            batch_index_mz = index_mz_arr[valid_indices].astype(int)
            batch_index_rt = index_rt_arr[valid_indices].astype(int)
            
            print(f"  Ready for chunked processing with vectorized operations")
            
            # All-at-once processing with row-wise HDF5 writes
            # One row at a time: extract, fill, write, discard (minimal memory)
            print(f"  Processing {len(valid_indices)} points with row-by-row strategy...")
            
            # Extract all valid data at once for efficient filtering
            all_rt = index_rt_arr[valid_indices].astype(int)
            all_mz = index_mz_arr[valid_indices].astype(int)
            all_intensity = spec_total[valid_indices, 2]
            
            # Compute RGB for all points (vectorized - very fast)
            all_log_intensity = np.log1p(all_intensity) if log_scale else all_intensity
            
            if invert_colors:
                all_norm = 1 - ((all_log_intensity - log_min) / log_diff)
            else:
                all_norm = (all_log_intensity - log_min) / log_diff
            
            all_norm = np.clip(all_norm, 0, 1)
            all_rgb_indices = (all_norm * 16777215).astype(np.uint32)
            all_rgb = intensity_to_rgb_vectorized(all_rgb_indices, rgb_array, scale_colors)
            
            # Get unique RT indices
            unique_rt = np.unique(all_rt)
            print(f"  Found {len(unique_rt)} RT indices to write")
            print(f"  Processing {len(unique_rt)} rows in batches of {row_batch_size}...")
            
            # Write rows in batches: prepare multiple rows, write together, discard
            total_written = 0
            num_rows = len(unique_rt)
            
            for batch_start_idx in range(0, num_rows, row_batch_size):
                batch_end_idx = min(batch_start_idx + row_batch_size, num_rows)
                batch_rt_indices = unique_rt[batch_start_idx:batch_end_idx]
                
                # Prepare all rows in this batch
                batch_rows_rgb = {}
                batch_rows_raw = {}
                
                for rt_idx in batch_rt_indices:
                    # Find indices of points for this RT
                    rt_mask = (all_rt == rt_idx)
                    rt_point_indices = np.where(rt_mask)[0]
                    
                    # Create row for this RT
                    row_rgb = np.full((img_width, 3), 255, dtype=np.uint8)
                    
                    # Fill m/z positions for this RT
                    for point_idx in rt_point_indices:
                        mz_idx = all_mz[point_idx]
                        row_rgb[mz_idx] = all_rgb[point_idx]
                    
                    batch_rows_rgb[rt_idx] = row_rgb
                
                # Write all rows in batch at once
                for rt_idx in batch_rt_indices:
                    img_frame_dset[rt_idx, :] = batch_rows_rgb[rt_idx]
                    total_written += np.sum(all_rt == rt_idx)
                
                # Progress counter
                progress_pct = 100 * batch_end_idx / num_rows
                print(f"  Writing rows: {batch_end_idx}/{num_rows} ({progress_pct:.1f}%) - {total_written} points written", end='\r')
                
                # Explicitly delete batch to free memory immediately
                del batch_rows_rgb
            
            print(f"\n  Successfully written all {total_written} points to HDF5")
            
            print(f"\n  Successfully saved heatmap to HDF5: {os.path.basename(hdf5_output)}")
    except Exception as e:
        print(f"  Error creating HDF5 file: {e}")
        print(f"  Attempting to clean up...")
        if os.path.exists(hdf5_output):
            try:
                os.remove(hdf5_output)
            except:
                pass
        raise
    
    # Optionally export PNG from HDF5
    if save_png:
        png_output = output_path if output_path.endswith('.png') else output_path.replace('.h5', '.png')
        export_heatmap_to_png(hdf5_output, png_output)
        
        # Generate HTML report
        if save_report:
            basename = os.path.splitext(os.path.basename(png_output))[0]
            generate_heatmap_report(png_output, hdf5_output, basename, spec_total, 
                                   log_scale, scale_colors, invert_colors, features_only=features_only)
    
    gc.collect()
    
    return hdf5_output


def main():
    parser = argparse.ArgumentParser(
        description="Create heatmap image from processed mzML spectral data"
    )
    parser.add_argument("--spectrum_pkl", help="Path to spectrum pickle file")
    parser.add_argument("--spectrum_hdf5", help="Path to spectrum HDF5 file")
    parser.add_argument("--output_h5", help="Output HDF5 heatmap file path (recommended)")
    parser.add_argument("--output_png", help="Output PNG heatmap file path (optional)")
    parser.add_argument("--export_png", action="store_true", help="Also export PNG from existing HDF5")
    parser.add_argument("--hdf5_input", help="Input HDF5 to export as PNG")
    parser.add_argument("--log_scale", type=bool, default=True, help="Apply logarithmic scaling")
    parser.add_argument("--scale_colors", type=bool, default=True, help="Scale colors based on intensities")
    parser.add_argument("--invert_colors", type=bool, default=False, help="Invert color scale")
    parser.add_argument("--round_up_to", type=int, default=2, help="Decimal places for m/z rounding")
    parser.add_argument("--batch_size", type=int, default=100000, help="Deprecated parameter (ignored)")
    parser.add_argument("--row_batch_size", type=int, default=10, help="Number of rows to prepare and write together (higher=faster but more RAM, each row ~3.2MB)")
    parser.add_argument("--compression_level", type=int, default=None, help="GZIP compression level: None=no compression (fastest, ~1min), 1=light (7min, 0.36MB), 4=balanced (13min, 0.2MB), 9=max (29min, 0.2MB)")
    parser.add_argument("--disable_batching", action="store_true", help="Deprecated parameter (ignored)")
    parser.add_argument("--features_only", action="store_true", help="Save only mapped feature points without heatmap background (creates much smaller file)")
    parser.add_argument("--save_report", action="store_true", default=True, help="Generate HTML report with heatmap and statistics (default: enabled)")
    parser.add_argument("--feature_data_json", help="Path to feature data JSON file for features-only visualization")
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print(f"Creating Heatmap Image")
    print(f"{'='*70}")
    
    # Mode 1: Export existing HDF5 to PNG
    if args.export_png and args.hdf5_input and args.output_png:
        export_heatmap_to_png(args.hdf5_input, args.output_png)
        print(f"Heatmap successfully exported to {args.output_png}")
        return
    
    # Mode 2: Load spectrum data and create heatmap
    # Load spectrum data from pickle or HDF5
    if args.spectrum_pkl:
        print("Loading spectrum data from pickle file...")
        with open(args.spectrum_pkl, 'rb') as f:
            data = pickle.load(f)
    elif args.spectrum_hdf5:
        print("Loading spectrum data from HDF5 file...")
        with h5py.File(args.spectrum_hdf5, 'r') as hf:
            data = {
                'spec_total': hf['spec_total'][:],
                'RTINSECONDS_arr': hf['RTINSECONDS_arr'][:],
                'mz_total_arr': hf['mz_total_arr'][:],
                'mz_dict': dict(hf['mz_dict'].attrs),
                'rt_dict': dict(hf['rt_dict'].attrs)
            }
    else:
        print("Error: Please provide either --spectrum_pkl or --spectrum_hdf5")
        return
    
    spec_total = data['spec_total']
    RTINSECONDS_arr = data['RTINSECONDS_arr']
    mz_total_arr = data['mz_total_arr']
    mz_dict = data['mz_dict']
    rt_dict = data['rt_dict']
    
    # Determine output paths
    if not args.output_h5 and not args.output_png:
        print("Error: Please provide either --output_h5 or --output_png")
        return
    
    # Default to HDF5 if no format specified
    output_path = args.output_h5 if args.output_h5 else args.output_png
    save_png = bool(args.output_png)
    
    # Create and save heatmap
    output_h5 = create_heatmap_image(
        spec_total, RTINSECONDS_arr, mz_total_arr, mz_dict, rt_dict,
        output_path,
        log_scale=args.log_scale,
        scale_colors=args.scale_colors,
        invert_colors=args.invert_colors,
        save_png=save_png,
        row_batch_size=args.row_batch_size,
        compression_level=args.compression_level,
        features_only=args.features_only,
        save_report=args.save_report,
        feature_data_json_path=args.feature_data_json
    )
    
    print(f"Heatmap successfully created and saved to {output_h5}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
from pdb import set_trace as bp
import argparse
import os
import sys
import pickle
import pandas as pd
import numpy as np
from PIL import Image
import plotly.graph_objects as go
from ast import literal_eval
import io
import base64


def argparse_setup():
    parser = argparse.ArgumentParser(
        description="Create an interactive Plotly chart of the intensity map with feature annotations"
    )
    parser.add_argument(
        "-heatmap_image", 
        required=True,
        help="Path to the heatmap PNG file"
    )
    parser.add_argument(
        "-features_tsv", 
        required=True,
        help="Path to the features TSV file"
    )
    parser.add_argument(
        "-mz_dict", 
        required=True,
        help="Path to the m/z dictionary pickle file"
    )
    parser.add_argument(
        "-rt_dict", 
        required=True,
        help="Path to the RT dictionary pickle file"
    )
    parser.add_argument(
        "-output_dir", 
        default=".",
        help="Output directory for the interactive HTML file"
    )
    
    return parser.parse_args()


def load_heatmap_as_base64(image_path):
    """Load heatmap image and convert to base64 for embedding in HTML"""
    with Image.open(image_path) as img:
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        img_base64 = base64.b64encode(buffered.getvalue()).decode()
    return img_base64


def parse_list(s):
    """Parse list strings from TSV"""
    if isinstance(s, str):
        try:
            parsed = literal_eval(s)
            # Handle both string and numeric values in the list
            result = []
            for x in parsed:
                try:
                    result.append(float(x))
                except (ValueError, TypeError):
                    pass
            return result if result else []
        except:
            return []
    elif isinstance(s, (list, tuple)):
        try:
            return [float(x) for x in s]
        except (ValueError, TypeError):
            return []
    return []


def create_interactive_map(heatmap_image, features_tsv, mz_dict_file, rt_dict_file, output_dir):
    """Create interactive Plotly heatmap visualization"""
    
    # Load dictionaries
    print("Loading coordinate dictionaries...")
    with open(mz_dict_file, 'rb') as f:
        mz_dict = pickle.load(f)
    with open(rt_dict_file, 'rb') as f:
        rt_dict = pickle.load(f)
    
    print(f"Loaded mz_dict with {len(mz_dict)} entries")
    print(f"Loaded rt_dict with {len(rt_dict)} entries")
    
    # Round RT values to 8 decimal places to ensure consistency
    rt_dict = {round(rt, 8): y for rt, y in rt_dict.items()}
    
    # Create a helper function to find nearest RT key
    def find_nearest_rt(rt_value, rt_dict):
        """Find the nearest RT value in the dictionary"""
        rounded_rt = round(rt_value, 8)
        if rounded_rt in rt_dict:
            return rt_dict[rounded_rt]
        # If not found, find the closest key
        rt_keys = list(rt_dict.keys())
        closest_key = min(rt_keys, key=lambda k: abs(k - rounded_rt))
        return rt_dict[closest_key]
    
    # Get axis ranges
    min_mz = min(mz_dict.keys())
    max_mz = max(mz_dict.keys())
    min_rt = min(rt_dict.keys())
    max_rt = max(rt_dict.keys())
    min_x = min(mz_dict.values())
    max_x = max(mz_dict.values())
    min_y = min(rt_dict.values())
    max_y = max(rt_dict.values())
    
    print(f"m/z range: {min_mz:.2f} - {max_mz:.2f}")
    print(f"RT range: {min_rt:.2f} - {max_rt:.2f}")
    print(f"Image coordinates: x({min_x}-{max_x}), y({min_y}-{max_y})")
    
    # Load heatmap image
    print("Loading heatmap image...")
    heatmap = Image.open(heatmap_image)
    img_width, img_height = heatmap.size
    print(f"Heatmap size: {img_width} x {img_height}")
    bp()
    # Convert image to base64 for embedding
    img_base64 = load_heatmap_as_base64(heatmap_image)
    
    # Load features data
    print("Loading features data...")
    features_data = pd.read_csv(features_tsv, sep='\t')
    print(f"Loaded {len(features_data)} features")
    
    # Extract feature information
    feature_boxes = []
    skipped_count = 0
    
    print("Extracting feature coordinates...")
    for idx, row in features_data.iterrows():
        mz_starts = parse_list(row['l_mz_start'])
        mz_ends = parse_list(row['l_mz_end'])
        rt_starts = parse_list(row['l_rt_start'])
        rt_ends = parse_list(row['l_rt_end'])
        
        if mz_starts and mz_ends and rt_starts and rt_ends:
            mz_start = min(mz_starts + mz_ends)
            mz_end = max(mz_starts + mz_ends)
            rt_start = min(rt_starts + rt_ends)
            rt_end = max(rt_starts + rt_ends)
            
            # Round to matching precision
            mz_precision = len(str(list(mz_dict.keys())[0]).split('.')[-1]) if mz_dict else 2
            
            mz_start_rounded = round(mz_start, mz_precision)
            mz_end_rounded = round(mz_end, mz_precision)
            
            try:
                x_start = mz_dict[mz_start_rounded]
                x_end = mz_dict[mz_end_rounded]
                y_start = find_nearest_rt(rt_start, rt_dict)
                y_end = find_nearest_rt(rt_end, rt_dict)
                
                x_min = int(min(x_start, x_end))
                x_max = int(max(x_start, x_end))
                y_min = int(min(y_start, y_end))
                y_max = int(max(y_start, y_end))
                
                feature_boxes.append({
                    'idx': idx,
                    'x_min': x_min,
                    'x_max': x_max,
                    'y_min': y_min,
                    'y_max': y_max,
                    'mz_start': mz_start,
                    'mz_end': mz_end,
                    'rt_start': rt_start,
                    'rt_end': rt_end,
                    'center_x': (x_min + x_max) / 2,
                    'center_y': (y_min + y_max) / 2
                })
                
                if len(feature_boxes) % 1000 == 0:
                    print(f"  Extracted {len(feature_boxes)} features...", end='\r')
            except KeyError:
                skipped_count += 1
                continue
        else:
            skipped_count += 1
    
    print(f"\nExtracted {len(feature_boxes)} features (skipped {skipped_count})")
    
    # Create hover text for features
    hover_texts = []
    for box in feature_boxes:
        text = (
            f"<b>Feature {box['idx']}</b><br>"
            f"m/z: {box['mz_start']:.4f} - {box['mz_end']:.4f}<br>"
            f"RT: {box['rt_start']:.4f} - {box['rt_end']:.4f}<br>"
            f"<extra></extra>"
        )
        hover_texts.append(text)
    
    # Create Plotly figure
    print("Creating interactive Plotly chart...")
    fig = go.Figure()
    
    # Add heatmap image as background
    fig.add_layout_image(
        dict(
            source=f"data:image/png;base64,{img_base64}",
            xref="x",
            yref="y",
            x=min_x,
            y=max_y,
            sizex=max_x - min_x,
            sizey=max_y - min_y,
            sizing="stretch",
            opacity=1.0,
            layer="below"
        )
    )
    
    # Add scatter plot for feature markers
    if feature_boxes:
        center_x = [box['center_x'] for box in feature_boxes]
        center_y = [box['center_y'] for box in feature_boxes]
        
        fig.add_trace(
            go.Scatter(
                x=center_x,
                y=center_y,
                mode='markers',
                marker=dict(
                    size=3,
                    color='rgba(255, 0, 0, 0.3)',
                    line=dict(width=0)
                ),
                hovertext=hover_texts,
                hoverinfo='text',
                showlegend=False,
                name='Features'
            )
        )
    
    # Update layout for interactivity
    fig.update_layout(
        title={
            'text': f'Interactive Intensity Map - {os.path.basename(heatmap_image).replace("_with_features_3.png", "").replace(".png", "")}',
            'x': 0.5,
            'xanchor': 'center'
        },
        xaxis=dict(
            title='m/z',
            range=[min_x, max_x],
            showgrid=False,
            scaleanchor="y",
            scaleratio=1,
            tickfont=dict(size=8)
        ),
        yaxis=dict(
            title='RT (min)',
            range=[max_y, min_y],  # Invert y-axis to match image orientation
            showgrid=False,
            scaleanchor="x",
            scaleratio=1,
            tickfont=dict(size=8)
        ),
        width=1400,
        height=1000,
        hovermode='closest',
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='white',
        margin=dict(l=80, r=20, t=60, b=80)
    )
    
    # Save the figure
    output_filename = os.path.join(
        output_dir, 
        os.path.basename(heatmap_image).replace('_with_features_3.png', '_interactive_plotly.html')
    )
    
    print(f"Saving interactive chart to {output_filename}")
    fig.write_html(output_filename)
    print(f"Done! Interactive map saved to {output_filename}")


if __name__ == "__main__":
    args = argparse_setup()
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    try:
        create_interactive_map(
            args.heatmap_image,
            args.features_tsv,
            args.mz_dict,
            args.rt_dict,
            args.output_dir
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

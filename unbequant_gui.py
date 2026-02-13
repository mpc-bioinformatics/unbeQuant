#!/usr/bin/env python3
"""
unbeQuant GUI - Interactive visualization and analysis tool
Integrates scripts from bin/ folder for heatmap display, feature overlay, and network analysis
"""

import os
import sys
import pickle
import json
import h5py
import numpy as np
import pandas as pd
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import base64
from io import BytesIO
import tempfile
import shutil

import dash
from dash import dcc, html, Input, Output, State, ctx, ALL, MATCH
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from PIL import Image

# Add bin to path for importing modules
BIN_DIR = Path(__file__).parent / 'bin'
sys.path.insert(0, str(BIN_DIR))


class UnbeQuantGUI:
    """Main GUI application class for unbeQuant"""
    
    def __init__(self, data_dir: str = "./gui_data"):
        """Initialize the GUI application"""
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(exist_ok=True)
        
        # Storage for loaded data
        self.heatmap_files = {}  # {filename: hdf5_path}
        self.feature_data = {}   # {filename: feature_list}
        self.spectrum_data = {}  # {filename: spectrum_pkl_data}
        self.paired_edges = {}   # {filename: edges_dict}
        
        # Current state
        self.current_heatmap = None
        self.current_features = None
        self.overlay_features = []
        
        # Cutoff parameters
        self.cutoff_mode = "euclidean"  # or "coordinate"
        self.euclidean_cutoff = 10.0
        self.mz_cutoff = 0.01
        self.rt_cutoff = 30.0
        
        # Initialize Dash app
        self.app = dash.Dash(
            __name__, 
            external_stylesheets=[dbc.themes.BOOTSTRAP],
            suppress_callback_exceptions=True
        )
        self.setup_layout()
        self.setup_callbacks()
    
    def setup_layout(self):
        """Setup the GUI layout"""
        self.app.layout = dbc.Container([
            dbc.Row([
                dbc.Col([
                    html.H1("unbeQuant - Interactive Feature Analysis", className="text-center mb-4"),
                ], width=12)
            ]),
            
            # Control panel
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader("Data Loading & Processing"),
                        dbc.CardBody([
                            # File upload section
                            dbc.Label("1. Upload mzML file:"),
                            dcc.Upload(
                                id='upload-mzml',
                                children=dbc.Button('Select mzML File', color='primary', size='sm'),
                                multiple=False
                            ),
                            html.Div(id='mzml-status', className='mt-2 mb-3'),
                            
                            dbc.Label("2. Upload TSV features file:"),
                            dcc.Upload(
                                id='upload-tsv',
                                children=dbc.Button('Select TSV File', color='primary', size='sm'),
                                multiple=False
                            ),
                            html.Div(id='tsv-status', className='mt-2 mb-3'),
                            
                            dbc.Button('Process Files', id='btn-process', color='success', className='w-100 mb-2'),
                            html.Div(id='process-status'),
                            
                            html.Hr(),
                            
                            # Heatmap selection
                            dbc.Label("Select Heatmap:"),
                            dcc.Dropdown(
                                id='heatmap-selector',
                                options=[],
                                placeholder="No heatmaps loaded"
                            ),
                            
                            html.Hr(),
                            
                            # Feature overlay controls
                            dbc.Label("Feature Overlays:"),
                            dcc.Checklist(
                                id='feature-overlay-checklist',
                                options=[],
                                value=[],
                                labelStyle={'display': 'block'}
                            ),
                        ])
                    ], className='mb-3'),
                    
                    dbc.Card([
                        dbc.CardHeader("Pairing & Network Settings"),
                        dbc.CardBody([
                            dbc.Label("Cutoff Mode:"),
                            dbc.RadioItems(
                                id='cutoff-mode',
                                options=[
                                    {'label': 'Euclidean Distance', 'value': 'euclidean'},
                                    {'label': 'Coordinate-based (m/z + RT)', 'value': 'coordinate'}
                                ],
                                value='euclidean',
                                inline=False
                            ),
                            
                            # Euclidean cutoff
                            html.Div([
                                dbc.Label("Euclidean Distance Cutoff:"),
                                dbc.Input(
                                    id='euclidean-cutoff',
                                    type='number',
                                    value=10.0,
                                    step=0.1,
                                    min=0.1
                                ),
                            ], id='euclidean-controls'),
                            
                            # Coordinate cutoffs
                            html.Div([
                                dbc.Label("m/z Cutoff:"),
                                dbc.Input(
                                    id='mz-cutoff',
                                    type='number',
                                    value=0.01,
                                    step=0.001,
                                    min=0.001
                                ),
                                dbc.Label("RT Cutoff (seconds):"),
                                dbc.Input(
                                    id='rt-cutoff',
                                    type='number',
                                    value=30.0,
                                    step=1.0,
                                    min=1.0
                                ),
                            ], id='coordinate-controls', style={'display': 'none'}),
                            
                            dbc.Button('Apply Pairing Settings', id='btn-apply-pairing', color='info', className='w-100 mt-3'),
                        ])
                    ])
                ], width=3),
                
                # Main display area
                dbc.Col([
                    dbc.Tabs([
                        dbc.Tab(label="Heatmap Viewer", tab_id="tab-heatmap"),
                        dbc.Tab(label="Network Graph", tab_id="tab-network"),
                        dbc.Tab(label="Diagnostic Plot", tab_id="tab-diagnostic"),
                    ], id="tabs", active_tab="tab-heatmap"),
                    html.Div(id='tab-content', className='mt-3')
                ], width=9)
            ]),
            
            # Hidden data stores
            dcc.Store(id='uploaded-mzml-data'),
            dcc.Store(id='uploaded-tsv-data'),
            dcc.Store(id='processed-data-store'),
            dcc.Store(id='selected-feature-store'),
            dcc.Store(id='heatmap-data-store'),
        ], fluid=True)
    
    def setup_callbacks(self):
        """Setup all callbacks for interactivity"""
        
        # Callback for cutoff mode visibility
        @self.app.callback(
            [Output('euclidean-controls', 'style'),
             Output('coordinate-controls', 'style')],
            [Input('cutoff-mode', 'value')]
        )
        def toggle_cutoff_controls(mode):
            if mode == 'euclidean':
                return {'display': 'block'}, {'display': 'none'}
            else:
                return {'display': 'none'}, {'display': 'block'}
        
        # Callback for file uploads
        @self.app.callback(
            Output('mzml-status', 'children'),
            Input('upload-mzml', 'contents'),
            State('upload-mzml', 'filename')
        )
        def upload_mzml(contents, filename):
            if contents is not None:
                return dbc.Alert(f"✓ Loaded: {filename}", color="success", dismissable=True)
            return ""
        
        @self.app.callback(
            Output('tsv-status', 'children'),
            Input('upload-tsv', 'contents'),
            State('upload-tsv', 'filename')
        )
        def upload_tsv(contents, filename):
            if contents is not None:
                return dbc.Alert(f"✓ Loaded: {filename}", color="success", dismissable=True)
            return ""
        
        # Main processing callback
        @self.app.callback(
            [Output('process-status', 'children'),
             Output('heatmap-selector', 'options'),
             Output('feature-overlay-checklist', 'options')],
            Input('btn-process', 'n_clicks'),
            [State('upload-mzml', 'contents'),
             State('upload-mzml', 'filename'),
             State('upload-tsv', 'contents'),
             State('upload-tsv', 'filename')]
        )
        def process_files(n_clicks, mzml_contents, mzml_filename, tsv_contents, tsv_filename):
            if n_clicks is None or mzml_contents is None or tsv_contents is None:
                return "", [], []
            
            try:
                # Save uploaded files
                mzml_path = self.save_uploaded_file(mzml_contents, mzml_filename)
                tsv_path = self.save_uploaded_file(tsv_contents, tsv_filename)
                
                # Process mzML file using script
                base_name = Path(mzml_filename).stem
                spectrum_pkl = self.data_dir / f"{base_name}_spectrum_data.pkl"
                
                if not self.process_mzml_to_spectrum(mzml_path, spectrum_pkl):
                    raise Exception("Failed to process mzML file")
                
                # Load spectrum data for later use
                with open(spectrum_pkl, 'rb') as f:
                    spectrum_data = pickle.load(f)
                self.spectrum_data[base_name] = spectrum_data
                
                # Create heatmap using script
                heatmap_png = self.data_dir / f"{base_name}_heatmap.png"
                
                if not self.create_heatmap_from_spectrum(spectrum_pkl, heatmap_png):
                    raise Exception("Failed to create heatmap")
                
                self.heatmap_files[base_name] = str(heatmap_png)
                
                # Extract features using script
                feature_pkl = self.data_dir / f"{base_name}_feature_data.pkl"
                feature_json = self.data_dir / f"{base_name}_feature_data.json"
                
                if not self.extract_features_from_tsv(tsv_path, feature_pkl, feature_json):
                    raise Exception("Failed to extract features")
                
                # Load feature data
                with open(feature_pkl, 'rb') as f:
                    feature_data = pickle.load(f)
                self.feature_data[base_name] = feature_data
                
                # Update dropdowns
                heatmap_options = [{'label': name, 'value': name} for name in self.heatmap_files.keys()]
                feature_options = [{'label': name, 'value': name} for name in self.feature_data.keys()]
                
                status = dbc.Alert(f"✓ Successfully processed {base_name}", color="success")
                return status, heatmap_options, feature_options
                
            except Exception as e:
                import traceback
                error_details = traceback.format_exc()
                print(f"Error processing files: {error_details}")
                error_msg = dbc.Alert(f"✗ Error: {str(e)}", color="danger")
                return error_msg, [], []
        
        # Tab content callback
        @self.app.callback(
            Output('tab-content', 'children'),
            [Input('tabs', 'active_tab'),
             Input('heatmap-selector', 'value'),
             Input('feature-overlay-checklist', 'value'),
             Input('selected-feature-store', 'data')]
        )
        def render_tab_content(active_tab, selected_heatmap, overlay_features, selected_feature):
            if active_tab == "tab-heatmap":
                return self.render_heatmap_tab(selected_heatmap, overlay_features)
            elif active_tab == "tab-network":
                return self.render_network_tab(selected_feature)
            elif active_tab == "tab-diagnostic":
                return self.render_diagnostic_tab(selected_feature)
            return html.Div("Select a tab")
    
    def save_uploaded_file(self, contents, filename):
        """Save uploaded file to data directory"""
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        
        file_path = self.data_dir / filename
        with open(file_path, 'wb') as f:
            f.write(decoded)
        
        return file_path
    
    def run_script(self, script_name: str, args: List[str]) -> Tuple[bool, str]:
        """Run a bin script via subprocess"""
        script_path = BIN_DIR / script_name
        cmd = [sys.executable, str(script_path)] + args
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=600  # 10 minute timeout
            )
            return True, result.stdout
        except subprocess.CalledProcessError as e:
            error_msg = f"Script failed: {e.stderr}"
            return False, error_msg
        except subprocess.TimeoutExpired:
            return False, "Script timed out after 10 minutes"
    
    def process_mzml_to_spectrum(self, mzml_path: Path, output_pkl: Path) -> bool:
        """Process mzML file to spectrum pickle"""
        args = [
            '--mzml', str(mzml_path),
            '--output_pickle', str(output_pkl),
            '--round_up_to', '2'
        ]
        success, output = self.run_script('process_mzml_file.py', args)
        if success:
            print(f"✓ Processed mzML: {output_pkl.name}")
        return success
    
    def create_heatmap_from_spectrum(self, spectrum_pkl: Path, output_png: Path) -> bool:
        """Create heatmap image from spectrum data"""
        args = [
            '--spectrum_pkl', str(spectrum_pkl),
            '--output_png', str(output_png),
            '--log_scale', 'True',
            '--scale_colors', 'True',
            '--row_batch_size', '10'
        ]
        success, output = self.run_script('create_heatmap_image_hdf5.py', args)
        if success:
            print(f"✓ Created heatmap: {output_png.name}")
        return success
    
    def extract_features_from_tsv(self, tsv_path: Path, output_pkl: Path, output_json: Path) -> bool:
        """Extract feature data from TSV file"""
        args = [
            '--tsv', str(tsv_path),
            '--output_pkl', str(output_pkl),
            '--output_json', str(output_json),
            '--round_up_to', '2',
            '--feature_mode', 'CoM',
            '--generate_diagnostic', 'False'
        ]
        success, output = self.run_script('extract_feature_data.py', args)
        if success:
            print(f"✓ Extracted features: {output_pkl.name}")
        return success
    
    def pair_features_multi(self, feature_dir: Path, output_pkl: Path, output_json: Path, 
                           cutoff_mode: str = 'euclidean', euclidean_cutoff: float = None,
                           mz_cutoff: float = None, rt_cutoff: float = None) -> bool:
        """Pair features across multiple files"""
        args = [
            '--input_dir', str(feature_dir),
            '--output_pkl', str(output_pkl),
            '--output_json', str(output_json),
            '--skip-matchfinder'  # We want just edges for network graph
        ]
        
        # Add cutoff parameters
        if cutoff_mode == 'euclidean' and euclidean_cutoff is not None:
            args.extend(['--edges_cutoff', str(euclidean_cutoff)])
        elif cutoff_mode == 'coordinate' and mz_cutoff is not None and rt_cutoff is not None:
            args.extend(['--mz_cutoff', str(mz_cutoff)])
            args.extend(['--rt_cutoff', str(rt_cutoff)])
        
        success, output = self.run_script('pair_features.py', args)
        if success:
            print(f"✓ Paired features: {output_pkl.name}")
        return success
    
    def render_heatmap_tab(self, selected_heatmap, overlay_features):
        """Render the heatmap viewer tab"""
        if selected_heatmap is None or selected_heatmap not in self.heatmap_files:
            return html.Div("No heatmap selected", className="text-center text-muted mt-5")
        
        # Load heatmap PNG image
        heatmap_path = self.heatmap_files[selected_heatmap]
        img = Image.open(heatmap_path)
        img_array = np.array(img)
        
        # Create Plotly figure with heatmap as image
        fig = go.Figure()
        
        # Add image as background
        fig.add_layout_image(
            dict(
                source=img,
                xref="x",
                yref="y",
                x=0,
                y=img_array.shape[0],
                sizex=img_array.shape[1],
                sizey=img_array.shape[0],
                sizing="stretch",
                layer="below"
            )
        )
        
        # Add feature boxes
        if selected_heatmap in self.feature_data:
            features = self.feature_data[selected_heatmap]
            self.add_feature_boxes(fig, features, selected_heatmap, color='rgba(255, 0, 0, 0.3)', name=selected_heatmap)
        
        # Add overlay features from other files
        for overlay_name in (overlay_features or []):
            if overlay_name in self.feature_data and overlay_name != selected_heatmap:
                features = self.feature_data[overlay_name]
                self.add_feature_boxes(fig, features, overlay_name, color='rgba(0, 0, 255, 0.3)', name=overlay_name)
        
        # Configure layout for zoom and pan
        fig.update_xaxes(range=[0, img_array.shape[1]], showgrid=False)
        fig.update_yaxes(range=[0, img_array.shape[0]], showgrid=False, scaleanchor="x", scaleratio=1)
        
        fig.update_layout(
            height=800,
            hovermode='closest',
            dragmode='pan',
            xaxis_title="m/z (index)",
            yaxis_title="RT (index)",
            showlegend=True,
            template='plotly_white'
        )
        
        return dcc.Graph(
            id='heatmap-graph',
            figure=fig,
            config={
                'scrollZoom': True,
                'displayModeBar': True,
                'modeBarButtonsToAdd': ['drawopenpath', 'eraseshape']
            },
            style={'height': '800px'}
        )
    
    def add_feature_boxes(self, fig, features, spectrum_name, color='rgba(255, 0, 0, 0.3)', name='Features'):
        """Add feature boxes to the figure"""
        # Get spectrum data for coordinate mapping
        spectrum_data = self.spectrum_data.get(spectrum_name)
        if spectrum_data is None:
            print(f"Warning: No spectrum data found for {spectrum_name}")
            return
        
        mz_dict = spectrum_data['mz_dict']
        rt_dict = spectrum_data['rt_dict']
        
        for feature in features:
            # Map coordinates to pixel indices
            try:
                x_min = mz_dict.get(round(feature['mz_start'], 2), 0)
                x_max = mz_dict.get(round(feature['mz_end'], 2), 0)
                y_min = rt_dict.get(round(feature['rt_start'], 8), 0)
                y_max = rt_dict.get(round(feature['rt_end'], 8), 0)
                
                # Create hover text
                pep_ident_str = str(feature.get('pep_ident', 'N/A'))
                if len(pep_ident_str) > 50:
                    pep_ident_str = pep_ident_str[:47] + "..."
                
                hover_text = (
                    f"Feature {feature['idx']}<br>"
                    f"x_center: {feature['x_center']:.4f}<br>"
                    f"y_center: {feature['y_center']:.4f}<br>"
                    f"pep_ident: {pep_ident_str}<br>"
                    f"charge: {feature.get('charge', 'N/A')}"
                )
                
                # Add rectangle
                fig.add_shape(
                    type="rect",
                    x0=x_min, y0=y_min,
                    x1=x_max, y1=y_max,
                    line=dict(color=color, width=1),
                    fillcolor=color,
                    name=name,
                    layer='above'
                )
                
                # Add invisible scatter point for hover and click
                fig.add_trace(go.Scatter(
                    x=[(x_min + x_max) / 2],
                    y=[(y_min + y_max) / 2],
                    mode='markers',
                    marker=dict(size=10, opacity=0),
                    hovertext=hover_text,
                    hoverinfo='text',
                    showlegend=False,
                    name=name,
                    customdata=[feature['idx']]  # Store feature index for click events
                ))
            except (KeyError, TypeError) as e:
                print(f"Warning: Could not map feature {feature.get('idx', 'unknown')}: {e}")
                continue
    
    def render_network_tab(self, selected_feature):
        """Render the network graph tab"""
        if selected_feature is None:
            return html.Div(
                "Click on a feature in the heatmap to view its network",
                className="text-center text-muted mt-5"
            )
        
        # TODO: Generate network graph for selected feature
        return html.Div([
            html.H4("Network Graph"),
            html.P(f"Feature: {selected_feature}"),
            html.P("Network graph will be displayed here")
        ])
    
    def render_diagnostic_tab(self, selected_feature):
        """Render the diagnostic plot tab"""
        if selected_feature is None:
            return html.Div(
                "Click on a feature in the heatmap to view its diagnostic plot",
                className="text-center text-muted mt-5"
            )
        
        # TODO: Generate diagnostic plot for selected feature
        return html.Div([
            html.H4("Diagnostic Plot"),
            html.P(f"Feature: {selected_feature}"),
            html.P("Diagnostic plot will be displayed here")
        ])
    
    def run(self, debug=True, port=8050):
        """Run the Dash application"""
        print(f"\n{'='*70}")
        print("unbeQuant GUI - Starting Application")
        print(f"{'='*70}")
        print(f"Data directory: {self.data_dir.absolute()}")
        print(f"Access the GUI at: http://localhost:{port}")
        print(f"{'='*70}\n")
        
        self.app.run_server(debug=debug, port=port, host='0.0.0.0')


def main():
    """Main entry point for the GUI"""
    import argparse
    
    parser = argparse.ArgumentParser(description="unbeQuant Interactive GUI")
    parser.add_argument("--data-dir", default="./gui_data", help="Directory for storing processed data")
    parser.add_argument("--port", type=int, default=8050, help="Port to run the web server on")
    parser.add_argument("--debug", action="store_true", help="Run in debug mode")
    
    args = parser.parse_args()
    
    gui = UnbeQuantGUI(data_dir=args.data_dir)
    gui.run(debug=args.debug, port=args.port)


if __name__ == "__main__":
    main()

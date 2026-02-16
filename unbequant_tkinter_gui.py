#!/usr/bin/env python3
"""
unbeQuant Tkinter GUI - Interactive visualization tool for large heatmaps
Local executable using tkinter to avoid HTML limitations with large files (+150MB)

Layout:
- Left/Middle: Large interactive heatmap with zoom/pan
- Right Top: Feature network graph (using graphviz)
- Right Bottom: Diagnostic plot for selected feature
"""

import os
import sys
import pickle
import json
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import subprocess
import threading
import tempfile

import numpy as np
import pandas as pd
from PIL import Image, ImageTk
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

# Add bin to path for importing modules
BIN_DIR = Path(__file__).parent / 'bin'
sys.path.insert(0, str(BIN_DIR))


class HeatmapCanvas:
    """Interactive heatmap viewer with zoom and pan capabilities"""
    
    def __init__(self, parent, on_feature_click=None):
        """Initialize the heatmap canvas"""
        self.parent = parent
        self.on_feature_click = on_feature_click
        
        # Create matplotlib figure
        self.figure = Figure(figsize=(10, 8), dpi=100)
        self.ax = self.figure.add_subplot(111)
        
        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.figure, master=parent)
        self.canvas.draw()
        
        # Add toolbar
        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()
        
        # Pack canvas
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Data storage
        self.heatmap_image = None
        self.features = []
        self.overlay_features = []
        self.feature_patches = []
        self.selected_feature = None
        
        # Coordinate mappings
        self.mz_dict = {}
        self.rt_dict = {}
        
        # Connect click event
        self.canvas.mpl_connect('button_press_event', self._on_click)
    
    def load_heatmap(self, image_path: str, mz_dict: Dict, rt_dict: Dict):
        """Load and display a heatmap image"""
        try:
            # Load image
            img = Image.open(image_path)
            img_array = np.array(img)
            
            # Store coordinate mappings
            self.mz_dict = mz_dict
            self.rt_dict = rt_dict
            
            # Get extents
            min_x = min(mz_dict.values()) if mz_dict else 0
            max_x = max(mz_dict.values()) if mz_dict else img_array.shape[1]
            min_y = min(rt_dict.values()) if rt_dict else 0
            max_y = max(rt_dict.values()) if rt_dict else img_array.shape[0]
            
            # Clear and display
            self.ax.clear()
            self.ax.imshow(img_array, aspect='auto', 
                          extent=[min_x, max_x, max_y, min_y],
                          interpolation='nearest')
            
            self.ax.set_xlabel('m/z', fontsize=10)
            self.ax.set_ylabel('RT (min)', fontsize=10)
            self.ax.set_title('MS1 Intensity Heatmap', fontsize=12)
            
            self.heatmap_image = img_array
            self.canvas.draw()
            
            return True
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load heatmap: {str(e)}")
            return False
    
    def add_features(self, features: List[Dict], color='red', alpha=0.3, is_overlay=False):
        """Add feature boxes to the heatmap"""
        if is_overlay:
            self.overlay_features.extend(features)
        else:
            self.features = features
        
        for feature in features:
            if 'x_min' in feature and 'x_max' in feature:
                x_min = feature['x_min']
                x_max = feature['x_max']
                y_min = feature['y_min']
                y_max = feature['y_max']
                
                # Create rectangle
                rect = Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                               linewidth=1, edgecolor=color, facecolor='none',
                               alpha=alpha, picker=True)
                rect.feature_data = feature
                self.ax.add_patch(rect)
                self.feature_patches.append(rect)
        
        self.canvas.draw()
    
    def clear_features(self):
        """Clear all feature overlays"""
        for patch in self.feature_patches:
            patch.remove()
        self.feature_patches = []
        self.features = []
        self.overlay_features = []
        self.canvas.draw()
    
    def _on_click(self, event):
        """Handle click events on the canvas"""
        if event.inaxes != self.ax:
            return
        
        # Find clicked feature
        for patch in self.feature_patches:
            if hasattr(patch, 'feature_data'):
                feature = patch.feature_data
                if ('x_min' in feature and 'x_max' in feature and
                    feature['x_min'] <= event.xdata <= feature['x_max'] and
                    feature['y_min'] <= event.ydata <= feature['y_max']):
                    
                    # Highlight selected feature
                    if self.selected_feature:
                        self.selected_feature.set_edgecolor('red')
                        self.selected_feature.set_linewidth(1)
                    
                    patch.set_edgecolor('yellow')
                    patch.set_linewidth(3)
                    self.selected_feature = patch
                    self.canvas.draw()
                    
                    # Callback
                    if self.on_feature_click:
                        self.on_feature_click(feature)
                    break


class NetworkGraphPanel:
    """Network graph visualization panel using graphviz"""
    
    def __init__(self, parent):
        """Initialize the network graph panel"""
        self.parent = parent
        
        # Create frame with canvas for displaying graphviz output
        self.canvas_frame = ttk.Frame(parent)
        self.canvas_frame.pack(fill=tk.BOTH, expand=True)
        
        # Label for status
        self.status_label = ttk.Label(self.canvas_frame, text="No feature selected")
        self.status_label.pack(pady=10)
        
        # Canvas for image display
        self.canvas = tk.Canvas(self.canvas_frame, bg='white')
        self.canvas.pack(fill=tk.BOTH, expand=True)
        
        # Scrollbars
        self.h_scrollbar = ttk.Scrollbar(self.canvas_frame, orient=tk.HORIZONTAL, command=self.canvas.xview)
        self.v_scrollbar = ttk.Scrollbar(self.canvas_frame, orient=tk.VERTICAL, command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=self.h_scrollbar.set, yscrollcommand=self.v_scrollbar.set)
        
        self.selected_feature = None
        self.current_image = None
        
        # Check if graphviz is available
        try:
            import graphviz
            self.graphviz_available = True
        except ImportError:
            self.graphviz_available = False
            self.status_label.config(text="graphviz not installed (pip install graphviz)")
    
    def display_network(self, selected_feature: Dict, all_features: List[Dict], 
                       edges: List[Dict], temp_dir: Path):
        """Display network graph for a selected feature using graphviz"""
        try:
            if not self.graphviz_available:
                self.status_label.config(text="graphviz not available - cannot render network")
                return
            
            import graphviz
            
            self.selected_feature = selected_feature
            self.status_label.config(text=f"Building network for feature {selected_feature.get('idx', 'N/A')}...")
            
            # Find connected features (subgraph extraction)
            selected_key = (selected_feature.get('filename', 'current'), 
                          selected_feature.get('idx', selected_feature.get('feature_idx', 0)))
            
            connected_nodes = set([selected_key])
            connected_edges = []
            
            for edge in edges:
                # Handle different edge formats
                if 'file1' in edge:
                    node1 = (edge['file1']['filename'], edge['file1']['feature_idx'])
                    node2 = (edge['file2']['filename'], edge['file2']['feature_idx'])
                else:
                    node1 = (edge.get('current_file', {}).get('filename', ''), 
                            edge.get('current_file', {}).get('feature_idx', 0))
                    node2 = (edge.get('matched_file', {}).get('filename', ''), 
                            edge.get('matched_file', {}).get('feature_idx', 0))
                
                if node1 == selected_key or node2 == selected_key:
                    connected_nodes.add(node1)
                    connected_nodes.add(node2)
                    connected_edges.append((node1, node2, edge.get('distance', 1.0)))
            
            if len(connected_nodes) == 1:
                self.status_label.config(text='No connections found for selected feature')
                self.canvas.delete("all")
                return
            
            # Create graphviz graph
            dot = graphviz.Digraph(comment='Feature Network', format='png')
            dot.attr(rankdir='LR', splines='polyline', overlap='false', sep='+0.3')
            dot.attr('node', shape='circle', style='filled', fontsize='10')
            
            # Generate colors for different files
            unique_files = sorted(set(f for f, _ in connected_nodes))
            file_colors = self._get_file_colors(unique_files)
            
            # Add nodes
            for node in connected_nodes:
                filename, feature_idx = node
                label = f"{Path(filename).stem}\n#{feature_idx}"
                
                # Color: red for selected, file-specific color for others
                if node == selected_key:
                    color = '#FF0000'  # Red for selected
                else:
                    color = file_colors.get(filename, '#87CEEB')
                
                dot.node(str(node), label=label, color=color, fillcolor=color)
            
            # Add edges
            for node1, node2, distance in connected_edges:
                edge_label = f"{distance:.3f}"
                dot.edge(str(node1), str(node2), label=edge_label, fontsize='8')
            
            # Render to temporary file
            output_path = temp_dir / f"network_graph_{selected_key[1]}"
            dot.render(output_path, cleanup=True)
            
            # Load and display the image
            png_path = f"{output_path}.png"
            if Path(png_path).exists():
                img = Image.open(png_path)
                
                # Resize if too large
                max_size = 800
                if img.width > max_size or img.height > max_size:
                    img.thumbnail((max_size, max_size), Image.Resampling.LANCZOS)
                
                self.current_image = ImageTk.PhotoImage(img)
                
                # Display on canvas
                self.canvas.delete("all")
                self.canvas.create_image(0, 0, anchor=tk.NW, image=self.current_image)
                self.canvas.config(scrollregion=self.canvas.bbox("all"))
                
                self.status_label.config(
                    text=f"Feature Network: {len(connected_nodes)} nodes, {len(connected_edges)} edges"
                )
                
                # Clean up temp file
                try:
                    Path(png_path).unlink()
                except:
                    pass
            else:
                self.status_label.config(text="Failed to render network graph")
            
        except Exception as e:
            import traceback
            print(f"Error displaying network: {traceback.format_exc()}")
            self.status_label.config(text=f'Error: {str(e)}')
            self.canvas.delete("all")
    
    def _get_file_colors(self, filenames: List[str]) -> Dict[str, str]:
        """Generate distinct colors for each file (matching build_network_graph.py)"""
        palette = [
            '#FF6B6B',  # Red
            '#4ECDC4',  # Teal
            '#45B7D1',  # Blue
            '#FFA07A',  # Light Salmon
            '#98D8C8',  # Mint
            '#F7DC6F',  # Yellow
            '#BB8FCE',  # Purple
            '#85C1E2',  # Sky Blue
            '#F8B88B',  # Peach
            '#AED6F1',  # Light Blue
        ]
        
        colors = {}
        for idx, filename in enumerate(sorted(set(filenames))):
            color_idx = idx % len(palette)
            colors[filename] = palette[color_idx]
        
        return colors


class DiagnosticPanel:
    """Diagnostic plot panel for individual features"""
    
    def __init__(self, parent):
        """Initialize the diagnostic panel"""
        self.parent = parent
        
        # Create matplotlib figure
        self.figure = Figure(figsize=(6, 5), dpi=80)
        self.ax = self.figure.add_subplot(111)
        
        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.figure, master=parent)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    def display_feature(self, feature: Dict):
        """Display diagnostic plot for a feature"""
        try:
            self.ax.clear()
            
            # Extract feature bounds
            x_min = feature.get('x_min', 0)
            x_max = feature.get('x_max', 100)
            y_min = feature.get('y_min', 0)
            y_max = feature.get('y_max', 100)
            
            # Draw bounding box
            rect = Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                           linewidth=2, edgecolor='blue', facecolor='lightblue',
                           alpha=0.3)
            self.ax.add_patch(rect)
            
            # Plot centers if available
            if 'center_x' in feature and 'center_y' in feature:
                self.ax.plot(feature['center_x'], feature['center_y'], 
                           'go', markersize=10, label='Geometric Center')
            
            if 'weighted_center_x' in feature and 'weighted_center_y' in feature:
                self.ax.plot(feature['weighted_center_x'], feature['weighted_center_y'],
                           'rx', markersize=10, markeredgewidth=2, 
                           label='Intensity-Weighted Center')
            
            # Set limits with padding
            padding = max((x_max - x_min) * 0.2, (y_max - y_min) * 0.2)
            self.ax.set_xlim(x_min - padding, x_max + padding)
            self.ax.set_ylim(y_min - padding, y_max + padding)
            
            # Labels
            self.ax.set_xlabel('m/z (pixels)', fontsize=10)
            self.ax.set_ylabel('RT (pixels)', fontsize=10)
            
            # Title with metadata
            title = f"Feature #{feature.get('idx', feature.get('feature_idx', 'N/A'))}"
            if 'charge' in feature:
                title += f" | Charge: {feature['charge']}"
            if 'pep_ident' in feature and feature['pep_ident']:
                pep = str(feature['pep_ident'])[:20]
                title += f"\nPeptide: {pep}"
            
            self.ax.set_title(title, fontsize=10)
            
            if 'center_x' in feature or 'weighted_center_x' in feature:
                self.ax.legend(fontsize=8, loc='upper right')
            
            self.ax.grid(True, alpha=0.3)
            
            self.canvas.draw()
            
        except Exception as e:
            self.ax.clear()
            self.ax.text(0.5, 0.5, f'Error displaying feature:\n{str(e)}',
                       ha='center', va='center', fontsize=10)
            self.canvas.draw()


class UnbeQuantTkinterGUI:
    """Main Tkinter GUI application"""
    
    def __init__(self, root, data_dir="./gui_data"):
        """Initialize the GUI"""
        self.root = root
        self.root.title("unbeQuant - Interactive Feature Analysis (Tkinter)")
        self.root.geometry("1600x900")
        
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(exist_ok=True)
        
        # Create temp directory for graphviz outputs
        self.temp_dir = self.data_dir / "temp"
        self.temp_dir.mkdir(exist_ok=True)
        
        # Data storage
        self.heatmap_files = {}
        self.feature_data = {}
        self.spectrum_data = {}
        self.paired_edges = None
        
        # Current state
        self.current_file = None
        self.selected_feature = None
        
        # Create UI
        self._create_menu()
        self._create_layout()
        
        # Status
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(root, textvariable=self.status_var, 
                              relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
    
    def _create_menu(self):
        """Create menu bar"""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        
        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Load mzML + TSV Files...", command=self.load_files)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        # View menu
        view_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="View", menu=view_menu)
        view_menu.add_command(label="Clear Features", command=self.clear_features)
        view_menu.add_command(label="Refresh Display", command=self._display_current_heatmap)
        
        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="About", command=self.show_about)
    
    def _create_layout(self):
        """Create the main layout"""
        # Main container
        main_container = ttk.Frame(self.root)
        main_container.pack(fill=tk.BOTH, expand=True)
        
        # Top control panel for file selection and overlays
        control_frame = ttk.Frame(main_container)
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        
        # Heatmap selector
        ttk.Label(control_frame, text="Select Heatmap:").pack(side=tk.LEFT, padx=5)
        self.heatmap_selector = ttk.Combobox(control_frame, width=30, state='readonly')
        self.heatmap_selector.pack(side=tk.LEFT, padx=5)
        self.heatmap_selector.bind('<<ComboboxSelected>>', self._on_heatmap_selected)
        
        ttk.Label(control_frame, text="  |  Feature Overlays:").pack(side=tk.LEFT, padx=5)
        
        # Overlay checkboxes frame
        self.overlay_frame = ttk.Frame(control_frame)
        self.overlay_frame.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        self.overlay_checkboxes = {}
        
        # Main panel container with three panels
        main_paned = ttk.PanedWindow(main_container, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True)
        
        # Left/Middle: Heatmap (larger)
        heatmap_frame = ttk.Frame(main_paned, width=900)
        main_paned.add(heatmap_frame, weight=3)
        
        # Right: Network graph (top) + Diagnostic (bottom)
        right_paned = ttk.PanedWindow(main_paned, orient=tk.VERTICAL)
        main_paned.add(right_paned, weight=1)
        
        # Network graph panel (top right)
        network_frame = ttk.LabelFrame(right_paned, text="Feature Network Graph", 
                                       height=400)
        right_paned.add(network_frame, weight=1)
        
        # Diagnostic panel (bottom right)
        diagnostic_frame = ttk.LabelFrame(right_paned, text="Diagnostic Plot",
                                         height=400)
        right_paned.add(diagnostic_frame, weight=1)
        
        # Create components
        self.heatmap_canvas = HeatmapCanvas(heatmap_frame, 
                                           on_feature_click=self.on_feature_selected)
        self.network_panel = NetworkGraphPanel(network_frame)
        self.diagnostic_panel = DiagnosticPanel(diagnostic_frame)
    
    def load_files(self):
        """Load multiple mzML and TSV files"""
        # Show instructions for multi-file selection
        response = messagebox.showinfo(
            "Multi-File Selection",
            "You can load one or multiple file pairs.\n\n"
            "How to select multiple files:\n"
            "• Windows/Linux: Hold Ctrl and click files\n"
            "• Mac: Hold Cmd and click files\n"
            "• For range: Hold Shift and click first and last file\n\n"
            "You will select mzML files first, then TSV files in the same order."
        )
        
        # Ask for mzML files (multiple selection)
        mzml_files = filedialog.askopenfilenames(
            title="Select one or more mzML files (Ctrl+Click or Cmd+Click for multiple)",
            filetypes=[("mzML files", "*.mzML"), ("All files", "*.*")]
        )
        
        if not mzml_files:
            return
        
        # Ask for TSV files (multiple selection, same count as mzML)
        tsv_files = filedialog.askopenfilenames(
            title=f"Select {len(mzml_files)} corresponding TSV files in same order (Ctrl+Click or Cmd+Click)",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        
        if not tsv_files:
            return
        
        if len(mzml_files) != len(tsv_files):
            messagebox.showerror(
                "File Count Mismatch",
                f"Number of mzML files ({len(mzml_files)}) must match "
                f"number of TSV files ({len(tsv_files)})"
            )
            return
        
        # Process all files in background thread
        thread = threading.Thread(target=self._process_multiple_files,
                                 args=(mzml_files, tsv_files))
        thread.daemon = True
        thread.start()
    
    def _process_multiple_files(self, mzml_files: tuple, tsv_files: tuple):
        """Process multiple mzML and TSV file pairs (runs in background thread)"""
        total_files = len(mzml_files)
        
        for idx, (mzml_file, tsv_file) in enumerate(zip(mzml_files, tsv_files), 1):
            try:
                self.status_var.set(f"Processing file {idx}/{total_files}: {Path(mzml_file).name}...")
                self._process_single_file(mzml_file, tsv_file)
            except Exception as e:
                import traceback
                error_msg = f"Error processing {Path(mzml_file).name}: {str(e)}\n{traceback.format_exc()}"
                print(error_msg)
                self.root.after(0, messagebox.showerror, "Error", 
                              f"Failed to process {Path(mzml_file).name}:\n{str(e)}")
        
        # Update UI on main thread
        self.root.after(0, self._update_file_selectors)
        self.status_var.set(f"Loaded {total_files} file(s)")
    
    def _process_single_file(self, mzml_file: str, tsv_file: str):
        """Process a single mzML and TSV file pair"""
        base_name = Path(mzml_file).stem
        
        # 1. Process mzML file to spectrum pickle
        spectrum_pkl = self.data_dir / f"{base_name}_spectrum_data.pkl"
        
        if not spectrum_pkl.exists():
            self._run_script('process_mzml_file.py', 
                           ['--mzml', mzml_file,
                            '--output_pickle', str(spectrum_pkl),
                            '--round_up_to', '2'])
        
        # Load spectrum data to get mz_dict and rt_dict
        with open(spectrum_pkl, 'rb') as f:
            spectrum_data = pickle.load(f)
        
        self.spectrum_data[base_name] = spectrum_data
        mz_dict = spectrum_data['mz_dict']
        rt_dict = spectrum_data['rt_dict']
        
        # 2. Create heatmap image
        heatmap_png = self.data_dir / f"{base_name}_heatmap.png"
        if not heatmap_png.exists():
            self._run_script('create_heatmap_image_hdf5.py',
                           ['--spectrum_pkl', str(spectrum_pkl),
                            '--output_png', str(heatmap_png),
                            '--log_scale', 'True',
                            '--scale_colors', 'True',
                            '--row_batch_size', '10'])
        
        # 3. Extract features from TSV
        feature_pkl = self.data_dir / f"{base_name}_feature_data.pkl"
        feature_json = self.data_dir / f"{base_name}_feature_data.json"
        if not feature_pkl.exists():
            self._run_script('extract_feature_data.py',
                           ['--tsv', tsv_file,
                            '--output_pkl', str(feature_pkl),
                            '--output_json', str(feature_json),
                            '--round_up_to', '2',
                            '--feature_mode', 'CoM',
                            '--generate_diagnostic', 'False'])
        
        # Load feature data
        with open(feature_pkl, 'rb') as f:
            features = pickle.load(f)
        
        # Map features to pixel coordinates for display
        mapped_features = self._map_features_to_pixels(features, mz_dict, rt_dict)
        
        self.feature_data[base_name] = mapped_features
        self.heatmap_files[base_name] = str(heatmap_png)
    
    def _map_features_to_pixels(self, features: List[Dict], mz_dict: Dict, rt_dict: Dict) -> List[Dict]:
        """Map feature coordinates to pixel indices for display"""
        mapped_features = []
        
        for feature in features:
            try:
                # Map m/z and RT coordinates to pixel indices
                mz_start = round(feature['mz_start'], 2)
                mz_end = round(feature['mz_end'], 2)
                rt_start = round(feature['rt_start'], 8)
                rt_end = round(feature['rt_end'], 8)
                
                x_min = mz_dict.get(mz_start, None)
                x_max = mz_dict.get(mz_end, None)
                y_min = rt_dict.get(rt_start, None)
                y_max = rt_dict.get(rt_end, None)
                
                if x_min is not None and x_max is not None and y_min is not None and y_max is not None:
                    # Create mapped feature with both original and pixel coordinates
                    mapped_feature = feature.copy()
                    mapped_feature['x_min'] = int(x_min)
                    mapped_feature['x_max'] = int(x_max)
                    mapped_feature['y_min'] = int(y_min)
                    mapped_feature['y_max'] = int(y_max)
                    mapped_feature['center_x'] = (x_min + x_max) / 2
                    mapped_feature['center_y'] = (y_min + y_max) / 2
                    
                    # Map weighted centers if available
                    if 'x_center' in feature and 'y_center' in feature:
                        # These are the intensity-weighted centers in m/z, RT space
                        mapped_feature['weighted_center_x'] = mz_dict.get(round(feature['x_center'], 2), mapped_feature['center_x'])
                        mapped_feature['weighted_center_y'] = rt_dict.get(round(feature['y_center'], 8), mapped_feature['center_y'])
                    
                    mapped_features.append(mapped_feature)
            except (KeyError, TypeError) as e:
                print(f"Warning: Could not map feature {feature.get('idx', 'unknown')}: {e}")
                continue
        
        return mapped_features
    
    def _update_file_selectors(self):
        """Update heatmap selector and overlay checkboxes (runs on main thread)"""
        # Update heatmap selector dropdown
        file_names = list(self.heatmap_files.keys())
        self.heatmap_selector['values'] = file_names
        
        if file_names:
            # Set first file as selected if nothing is selected
            if not self.heatmap_selector.get() and file_names:
                self.heatmap_selector.current(0)
                self.current_file = file_names[0]
                self._display_current_heatmap()
        
        # Update overlay checkboxes
        # Clear existing checkboxes
        for widget in self.overlay_frame.winfo_children():
            widget.destroy()
        self.overlay_checkboxes.clear()
        
        # Create checkbox for each file except the currently displayed one
        for file_name in file_names:
            var = tk.BooleanVar(value=False)
            cb = ttk.Checkbutton(
                self.overlay_frame,
                text=file_name,
                variable=var,
                command=self._on_overlay_changed
            )
            cb.pack(side=tk.LEFT, padx=2)
            self.overlay_checkboxes[file_name] = var
    
    def _on_heatmap_selected(self, event=None):
        """Handle heatmap selection change"""
        selected = self.heatmap_selector.get()
        if selected and selected in self.heatmap_files:
            self.current_file = selected
            self._display_current_heatmap()
    
    def _display_current_heatmap(self):
        """Display the currently selected heatmap with overlays"""
        if not self.current_file or self.current_file not in self.heatmap_files:
            return
        
        heatmap_path = self.heatmap_files[self.current_file]
        spectrum_data = self.spectrum_data[self.current_file]
        features = self.feature_data[self.current_file]
        
        mz_dict = spectrum_data['mz_dict']
        rt_dict = spectrum_data['rt_dict']
        
        # Load and display heatmap
        if self.heatmap_canvas.load_heatmap(heatmap_path, mz_dict, rt_dict):
            # Clear previous features
            self.heatmap_canvas.clear_features()
            
            # Add features from current file (red)
            self.heatmap_canvas.add_features(features, color='red', alpha=0.5, is_overlay=False)
            
            # Add overlay features from other files
            self._update_overlays()
            
            self.status_var.set(f"Displaying {self.current_file}: {len(features)} features")
    
    def _on_overlay_changed(self):
        """Handle overlay checkbox changes"""
        self._update_overlays()
    
    def _update_overlays(self):
        """Update overlay features based on checkbox selections"""
        if not self.current_file:
            return
        
        # Clear existing overlays (keep main features)
        # Get current main features
        main_features = self.feature_data.get(self.current_file, [])
        
        # Clear all and re-add main features
        self.heatmap_canvas.clear_features()
        self.heatmap_canvas.add_features(main_features, color='red', alpha=0.5, is_overlay=False)
        
        # Get current heatmap's coordinate system
        spectrum_data = self.spectrum_data[self.current_file]
        current_mz_dict = spectrum_data['mz_dict']
        current_rt_dict = spectrum_data['rt_dict']
        
        # Add overlays from checked files
        overlay_colors = ['blue', 'green', 'orange', 'purple', 'cyan', 'magenta']
        color_idx = 0
        
        for file_name, var in self.overlay_checkboxes.items():
            if var.get() and file_name != self.current_file and file_name in self.feature_data:
                # Get features from overlay file
                overlay_features = self.feature_data[file_name]
                
                # Remap features to current heatmap's coordinate system if different
                if file_name in self.spectrum_data:
                    overlay_spectrum_data = self.spectrum_data[file_name]
                    overlay_mz_dict = overlay_spectrum_data['mz_dict']
                    overlay_rt_dict = overlay_spectrum_data['rt_dict']
                    
                    # If coordinate systems differ, remap
                    if overlay_mz_dict != current_mz_dict or overlay_rt_dict != current_rt_dict:
                        # Remap overlay features to current coordinate system
                        remapped_features = []
                        for feature in overlay_features:
                            try:
                                # Get original m/z and RT values
                                mz_start = feature.get('mz_start', None)
                                mz_end = feature.get('mz_end', None)
                                rt_start = feature.get('rt_start', None)
                                rt_end = feature.get('rt_end', None)
                                
                                if all(v is not None for v in [mz_start, mz_end, rt_start, rt_end]):
                                    # Map to current heatmap's pixel coordinates
                                    x_min = current_mz_dict.get(round(mz_start, 2), None)
                                    x_max = current_mz_dict.get(round(mz_end, 2), None)
                                    y_min = current_rt_dict.get(round(rt_start, 8), None)
                                    y_max = current_rt_dict.get(round(rt_end, 8), None)
                                    
                                    if all(v is not None for v in [x_min, x_max, y_min, y_max]):
                                        remapped_feature = feature.copy()
                                        remapped_feature['x_min'] = int(x_min)
                                        remapped_feature['x_max'] = int(x_max)
                                        remapped_feature['y_min'] = int(y_min)
                                        remapped_feature['y_max'] = int(y_max)
                                        remapped_features.append(remapped_feature)
                            except Exception as e:
                                print(f"Warning: Could not remap overlay feature: {e}")
                                continue
                        
                        overlay_features = remapped_features
                
                # Add overlay with distinct color
                color = overlay_colors[color_idx % len(overlay_colors)]
                self.heatmap_canvas.add_features(overlay_features, color=color, alpha=0.3, is_overlay=True)
                color_idx += 1
    
    def _run_script(self, script_name: str, args: List[str]):
        """Run a bin script"""
        script_path = BIN_DIR / script_name
        cmd = [sys.executable, str(script_path)] + args
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"Script {script_name} failed:\n{result.stderr}")
    
    def on_feature_selected(self, feature: Dict):
        """Handle feature selection"""
        self.selected_feature = feature
        
        # Update diagnostic panel
        self.diagnostic_panel.display_feature(feature)
        
        # Update network panel if we have edges
        if self.paired_edges:
            all_features = []
            for features_list in self.feature_data.values():
                all_features.extend(features_list)
            
            self.network_panel.display_network(feature, all_features, 
                                              self.paired_edges, self.temp_dir)
        else:
            # Try to load edges if multiple files exist
            if len(self.feature_data) > 1:
                self._load_or_create_edges()
    
    def _load_or_create_edges(self):
        """Load or create feature pairing edges"""
        edges_pkl = self.data_dir / "paired_features_edges.pkl"
        
        if edges_pkl.exists():
            with open(edges_pkl, 'rb') as f:
                self.paired_edges = pickle.load(f)
        else:
            # Need to create edges - show message
            response = messagebox.askyesno(
                "Create Pairings",
                "Feature pairings not found. Create them now?\n"
                "(This may take a few minutes)")
            
            if response:
                thread = threading.Thread(target=self._create_pairings)
                thread.daemon = True
                thread.start()
    
    def _create_pairings(self):
        """Create feature pairings in background"""
        try:
            self.status_var.set("Creating feature pairings...")
            
            # Get all feature pkl files
            if len(self.feature_data) < 2:
                self.root.after(0, messagebox.showinfo, "Info",
                              "Need at least 2 files for pairing")
                self.status_var.set("Need at least 2 files")
                return
            
            edges_pkl = self.data_dir / "paired_features_edges.pkl"
            edges_json = self.data_dir / "paired_features_edges.json"
            
            # Run pairing script with input_dir pointing to data_dir
            # The script will find all feature_data.pkl files there
            args = ['--input_dir', str(self.data_dir),
                    '--output_pkl', str(edges_pkl),
                    '--output_json', str(edges_json),
                    '--skip-matchfinder']  # We only need edges for network graph
            
            self._run_script('pair_features.py', args)
            
            # Load edges
            with open(edges_pkl, 'rb') as f:
                paired_data = pickle.load(f)
            
            # Extract edges from paired data
            # The format can be dict with file indices as keys
            if isinstance(paired_data, dict):
                self.paired_edges = []
                for file_idx, file_edges in paired_data.items():
                    if isinstance(file_edges, list):
                        self.paired_edges.extend(file_edges)
                    else:
                        self.paired_edges.append(file_edges)
            elif isinstance(paired_data, list):
                self.paired_edges = paired_data
            else:
                self.paired_edges = []
            
            self.status_var.set(f"Loaded {len(self.paired_edges)} feature pairings")
            
            # Refresh network if feature selected
            if self.selected_feature:
                self.root.after(0, self.on_feature_selected, self.selected_feature)
            
        except Exception as e:
            import traceback
            error_msg = f"Error creating pairings: {str(e)}\n{traceback.format_exc()}"
            print(error_msg)
            self.root.after(0, messagebox.showerror, "Error", str(e))
            self.status_var.set("Error creating pairings")
    
    def clear_features(self):
        """Clear all features from display"""
        self.heatmap_canvas.clear_features()
        self.selected_feature = None
        self.status_var.set("Features cleared")
    
    def show_about(self):
        """Show about dialog"""
        messagebox.showinfo(
            "About unbeQuant GUI",
            "unbeQuant Interactive Feature Analysis\n"
            "Tkinter-based local GUI for large heatmap visualization\n\n"
            "Features:\n"
            "- Load and compare multiple mzML/TSV file pairs\n"
            "- Interactive heatmap viewer with zoom/pan\n"
            "- Switch between different heatmaps\n"
            "- Overlay features from multiple files\n"
            "- Feature network graph visualization\n"
            "- Diagnostic plots for selected features\n\n"
            "Click on features in the heatmap to view details."
        )


def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="unbeQuant Tkinter GUI - Local executable for large heatmap visualization"
    )
    parser.add_argument('--data-dir', default='./gui_data',
                       help='Directory for storing processed data')
    
    args = parser.parse_args()
    
    # Create and run GUI
    root = tk.Tk()
    app = UnbeQuantTkinterGUI(root, data_dir=args.data_dir)
    root.mainloop()


if __name__ == "__main__":
    main()

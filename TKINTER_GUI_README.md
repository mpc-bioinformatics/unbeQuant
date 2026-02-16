# unbeQuant Tkinter GUI - Local Executable for Large Heatmaps

A desktop application using Tkinter for visualizing and analyzing mass spectrometry feature data from the unbeQuant workflow. This implementation is designed to handle very large heatmaps (+150MB) that are problematic for web-based solutions.

## Why Tkinter Instead of HTML/Dash?

The previous HTML-based GUI using Dash struggled with large heatmaps (+150MB each) because:
- Browsers have memory limitations when embedding large base64-encoded images
- HTML rendering performance degrades with very large images
- Network transfer overhead for web-based solutions

The Tkinter GUI solves these issues by:
- Running as a native local executable
- Using matplotlib's efficient image rendering
- Leveraging native OS graphics capabilities
- No browser memory constraints

## Features

### Layout
- **Left/Middle Panel**: Large interactive heatmap with zoom and pan
  - High-resolution MS1 spectra visualization
  - Feature overlay with bounding boxes
  - Click to select features
  - Matplotlib navigation toolbar for zoom/pan
  
- **Right Top Panel**: Feature network graph
  - Shows connections between features across files
  - Selected feature highlighted in red
  - Connected features in light blue
  - NetworkX-based graph layout
  
- **Right Bottom Panel**: Diagnostic plot
  - Feature bounding box visualization
  - Geometric center (green circle)
  - Intensity-weighted center (red X)
  - Feature metadata display

### Capabilities
- **Interactive Heatmap Viewing**: Zoom, pan, and explore large heatmaps efficiently
- **Feature Selection**: Click on features to view details
- **Network Analysis**: View feature connections across multiple files
- **Diagnostic Information**: Detailed view of individual feature characteristics
- **Multi-file Support**: Load and compare features from multiple runs

## Installation

### Prerequisites

- Python 3.9+
- All dependencies from the main unbeQuant workflow
- GUI-specific dependencies

### Install Dependencies

```bash
pip install -r gui_requirements.txt
```

Key dependencies for Tkinter GUI:
- **matplotlib** - For interactive heatmap rendering and plots
- **graphviz** - For network graph visualization (matching build_network_graph.py)
- **Pillow** - For image loading
- **numpy, pandas** - For data processing
- **tkinter** - Usually included with Python

**Note**: graphviz also requires the Graphviz system package:
- Ubuntu/Debian: `sudo apt-get install graphviz`
- macOS: `brew install graphviz`
- Windows: Download from https://graphviz.org/download/

## Usage

### Quick Start

1. **Launch the GUI**:
```bash
python3 unbequant_tkinter_gui.py
```

Or use the launcher script:
```bash
./launch_tkinter_gui.sh
```

2. **Load Data**:
   - Go to `File → Load mzML + TSV`
   - Select an mzML file from your unbeQuant workflow
   - Select the corresponding TSV features file
   - Click OK - processing will begin automatically

3. **Explore Features**:
   - Use the matplotlib toolbar to zoom and pan
   - Click on feature boxes to select them
   - View network connections in the top right panel
   - View diagnostic information in the bottom right panel

### Command-Line Options

```bash
python3 unbequant_tkinter_gui.py --help

Options:
  --data-dir DATA_DIR  Directory for storing processed data (default: ./gui_data)
```

### Example

```bash
# Run GUI with custom data directory
python3 unbequant_tkinter_gui.py --data-dir /path/to/my/gui_data
```

## Workflow

### 1. Data Processing Pipeline

When you load files, the GUI automatically:

1. **process_mzml_file.py**: Extracts spectral data
   - Reads MS1 spectra (retention time, m/z, intensity)
   - Creates indexed dictionaries for coordinate mapping
   - Outputs: `*_spectrum_data.pkl`, `*_mz_dict.pkl`, `*_rt_dict.pkl`

2. **create_heatmap_image_hdf5.py**: Generates heatmap
   - Creates high-resolution PNG from spectral data
   - Uses logarithmic intensity scaling
   - Outputs: `*_heatmap.png`

3. **extract_feature_data.py**: Extracts features
   - Parses feature boundaries and metadata
   - Maps features to pixel coordinates
   - Outputs: `*_feature_data.pkl`

4. **pair_features.py** (on-demand): Creates feature pairings
   - Only runs when viewing network graphs
   - Uses KD-tree spatial matching
   - Outputs: `paired_features_edges.pkl`

### 2. Viewing Heatmaps

- **Heatmap Display**: Main panel shows the MS1 intensity heatmap
- **Zoom/Pan**: Use matplotlib toolbar:
  - 🏠 Home - Reset view
  - ← → Pan/zoom controls
  - 🔍 Zoom to rectangle
  - 💾 Save figure
- **Feature Boxes**: Features overlaid as red rectangles
- **Click to Select**: Click on a feature box to select it (turns yellow)

### 3. Feature Selection

When you click a feature:
- **Heatmap**: Feature box highlights in yellow
- **Network Panel**: Shows all connected features across files
- **Diagnostic Panel**: Shows detailed feature information

### 4. Network Graph

- Displays connections between the selected feature and related features
- Selected feature: Red node
- Connected features: Light blue nodes
- Edges show pairing relationships
- Labels show file name and feature index

### 5. Diagnostic Plot

For the selected feature, shows:
- Bounding box (blue rectangle with light blue fill)
- Geometric center (green circle)
- Intensity-weighted center (red X)
- Feature metadata in title (charge, peptide identification)

## Data Directory Structure

Processed data is stored in the specified data directory (default: `./gui_data`):

```
gui_data/
├── <sample>_spectrum_data.pkl      # Extracted spectral data
├── <sample>_mz_dict.pkl            # m/z coordinate mapping
├── <sample>_rt_dict.pkl            # RT coordinate mapping
├── <sample>_heatmap.png            # Heatmap image
├── <sample>_feature_data.pkl       # Feature data
├── paired_features_edges.pkl       # Network edges (created on-demand)
└── paired_features_edges.json      # Network edges (human-readable)
```

## Performance Considerations

### Large Heatmaps
- PNG files can be very large (+150MB)
- Initial load may take a few seconds
- Once loaded, matplotlib handles rendering efficiently
- Zoom/pan operations are smooth
- Much better performance than browser-based rendering

### Memory Usage
- Base application: ~50-100MB
- Per loaded heatmap: ~200-300MB (in memory)
- Feature data: ~50-100MB per file
- Total: Plan for ~500MB-1GB per loaded dataset

### Processing Time
- mzML processing: 1-10 minutes depending on file size
- Heatmap generation: 2-5 minutes
- Feature extraction: <1 minute
- Processed files are cached - subsequent loads are instant

## Troubleshooting

### GUI won't start
- **Tkinter not found**: Tkinter is usually included with Python, but on some systems:
  - Ubuntu/Debian: `sudo apt-get install python3-tk`
  - macOS: Should be included with Python
  - Windows: Should be included with Python
- **Missing dependencies**: `pip install -r gui_requirements.txt`

### Matplotlib errors
- If you get display errors, try setting the backend:
  ```python
  export MPLBACKEND=TkAgg
  ```

### Processing fails
- Ensure mzML and TSV files are from unbeQuant workflow
- Check that bin scripts are executable
- Verify sufficient disk space in data directory

### Heatmap not displaying
- Check for errors in terminal output
- Verify PNG file was created in gui_data directory
- Try loading a smaller dataset first

### Network graph empty
- Load at least 2 files first
- Click a feature to trigger pairing
- Check pairing creation dialog

### Performance issues
- Close other applications to free memory
- Process one file at a time
- Clear old data from gui_data directory

## Differences from Dash GUI

| Feature | Dash GUI | Tkinter GUI |
|---------|----------|-------------|
| Platform | Web browser | Native desktop |
| Max heatmap size | ~50-100MB | 150MB+ |
| Performance | Limited by browser | Native rendering |
| Installation | pip install dash | pip install matplotlib |
| Network access | localhost:8050 | Local only |
| UI Framework | React/Plotly | Tkinter/Matplotlib |
| Deployment | Can be deployed | Local exe only |

## Integration with unbeQuant Workflow

The Tkinter GUI works with output from the main unbeQuant workflow:

1. **Run the main workflow** to generate:
   - mzML files: `results/mgfs_mzmls/*.mzML`
   - TSV features: `results/quantifications/features_with_annotated_identifications/*.tsv`

2. **Launch the Tkinter GUI**

3. **Load files** via File → Load mzML + TSV

4. **Explore** your results interactively

## Reference Data

Example data for testing the GUI:

- **mzML files**: https://1drv.ms/f/c/a41a7100db890e3c/IgBzMEItLCItRIWzNsDbiUEHATEaIxwe2ykr48puJZrPh4c?e=DfS5Vo
- **TSV feature files**: https://1drv.ms/f/c/a41a7100db890e3c/IgADiGBIEdXYTpGSyeVyzCRuAYB0sUKbbo9DRq6zqFdG9J4?e=6Fa1Xg

Download corresponding mzML and TSV files to test the GUI with real data.

## Future Enhancements

Potential improvements:
- [ ] Batch file loading
- [ ] Export network graphs as images
- [ ] 3D intensity visualization
- [ ] Feature annotation tools
- [ ] Session save/load
- [ ] Progress bars for processing
- [ ] Configurable pairing cutoffs in UI
- [ ] Feature comparison view
- [ ] Export selected features to CSV

## Technical Details

### Architecture
- **GUI Framework**: Tkinter (Python's standard GUI library)
- **Visualization**: Matplotlib with TkAgg backend
- **Network Graphs**: Graphviz (matching build_network_graph.py implementation)
- **Image Handling**: PIL/Pillow
- **Data Processing**: NumPy, Pandas

### Why This Stack?
- **Tkinter**: Built-in, no installation needed, cross-platform
- **Matplotlib**: Excellent for scientific visualization, handles large images
- **Graphviz**: Professional graph visualization library, same as bin/build_network_graph.py
- **Native Performance**: No browser overhead, direct OS rendering

### Thread Safety
- File processing runs in background threads
- UI updates are marshaled to main thread
- Prevents GUI freezing during long operations

## Contributing

This Tkinter GUI is part of the unbeQuant project. For issues, suggestions, or contributions, please refer to the main project repository.

## License

Same license as the main unbeQuant project.

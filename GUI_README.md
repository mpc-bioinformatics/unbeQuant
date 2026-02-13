# unbeQuant Interactive GUI

A web-based interactive GUI for visualizing and analyzing mass spectrometry feature data from the unbeQuant workflow.

## Features

- **Heatmap Visualization**: High-resolution, zoomable heatmap viewer for MS1 spectra
- **Feature Overlay**: Display and compare features from multiple files with color coding
- **Interactive Feature Boxes**: Hover over features to see metadata (m/z center, RT center, peptide identifications, charge state)
- **Network Graph**: Click on features to view connected features across files
- **Diagnostic Plots**: Visualize individual feature characteristics and centroid calculations
- **Flexible Pairing**: Configure feature pairing cutoffs using Euclidean distance or separate m/z and RT thresholds

## Installation

### Prerequisites

- Python 3.9+
- All dependencies from the main unbeQuant workflow
- GUI-specific dependencies (install with):

```bash
pip install -r gui_requirements.txt
```

### Dependencies

The GUI requires:
- Dash (web framework)
- Plotly (interactive visualizations)
- Dash Bootstrap Components (UI styling)
- NetworkX (network graph generation)
- Plus all standard unbeQuant dependencies (pyOpenMS, h5py, pandas, numpy, etc.)

## Usage

### Quick Start

1. **Start the GUI**:
```bash
python3 unbequant_gui.py
```

2. **Access the interface**: Open your web browser and navigate to:
```
http://localhost:8050
```

3. **Process your data**:
   - Upload an mzML file (from your unbeQuant workflow)
   - Upload the corresponding TSV features file (from `results/quantifications/features_with_annotated_identifications/`)
   - Click "Process Files" to generate heatmaps and extract features

### Command-Line Options

```bash
python3 unbequant_gui.py --help

Options:
  --data-dir DATA_DIR  Directory for storing processed data (default: ./gui_data)
  --port PORT          Port to run the web server on (default: 8050)
  --debug              Run in debug mode (enables auto-reload on code changes)
```

### Example

```bash
# Run GUI with custom data directory and port
python3 unbequant_gui.py --data-dir /path/to/my/gui_data --port 8080 --debug
```

## Workflow

### 1. Data Processing Pipeline

The GUI integrates the following scripts from the `bin/` folder:

1. **process_mzml_file.py**: Extracts spectral data from mzML files
   - Reads MS1 spectra (retention time, m/z, intensity)
   - Creates indexed dictionaries for fast coordinate lookup
   - Outputs: `*_spectrum_data.pkl`

2. **create_heatmap_image_hdf5.py**: Generates heatmap images
   - Creates high-resolution PNG heatmaps from spectral data
   - Uses logarithmic intensity scaling and color mapping
   - Outputs: `*_heatmap.png`

3. **extract_feature_data.py**: Extracts feature information from TSV files
   - Parses feature boundaries, centers, and metadata
   - Calculates both geometric and intensity-weighted centroids
   - Outputs: `*_feature_data.pkl`, `*_feature_data.json`

4. **pair_features.py**: Pairs features across multiple files
   - Uses KD-tree based spatial indexing for efficient matching
   - Creates network edges between connected features
   - Supports both Euclidean distance and coordinate-based filtering
   - Outputs: `paired_features_edges.pkl`

### 2. Viewing Heatmaps

- **Heatmap Selection**: Choose from loaded heatmaps in the dropdown
- **Zoom and Pan**: Use mouse to zoom (scroll wheel) and pan (drag)
- **Feature Boxes**: Features are overlaid as colored rectangles
  - Red boxes: Features from the current file
  - Blue boxes: Features from overlay files
- **Hover Information**: Hover over feature boxes to see:
  - Feature index
  - m/z center (x_center)
  - RT center (y_center)
  - Peptide identifications (pep_ident)
  - Charge state

### 3. Feature Overlay

- **Add Overlays**: Check boxes in "Feature Overlays" to overlay features from other files
- **Color Coding**: Features from different files use different colors
- **Comparison**: Visually compare feature locations across runs

### 4. Feature Pairing Settings

Configure how features are matched across files:

#### Euclidean Distance Mode
- Single cutoff parameter (default: 10.0)
- Features within this distance are considered matches
- Distance calculated in scaled m/z-RT space

#### Coordinate-based Mode
- Separate m/z cutoff (default: 0.01 Da)
- Separate RT cutoff (default: 30.0 seconds)
- Features must be within both thresholds to match

**Click "Apply Pairing Settings"** to regenerate feature pairings with new parameters.

### 5. Network Graph View

- **Click a feature** in the heatmap to select it
- Switch to the **"Network Graph" tab**
- View all features connected to the selected feature
- The selected feature is highlighted in red
- Connected features are shown in light blue
- Edges show pairing relationships

### 6. Diagnostic Plots

- **Click a feature** in the heatmap to select it
- Switch to the **"Diagnostic Plot" tab**
- View:
  - Feature bounding box
  - Geometric center (green circle)
  - Intensity-weighted center (red X)
  - Charge state and peptide identification

## Data Directory Structure

The GUI stores processed data in the specified data directory (default: `./gui_data`):

```
gui_data/
├── <sample>_spectrum_data.pkl      # Extracted spectral data
├── <sample>_heatmap.png            # Heatmap image
├── <sample>_feature_data.pkl       # Feature data (binary)
├── <sample>_feature_data.json      # Feature data (human-readable)
├── paired_features_edges.pkl       # Network edges
└── paired_features_edges.json      # Network edges (human-readable)
```

## Performance Considerations

### Large Heatmaps
- Heatmaps can be ~150MB in size
- PNG format is used for fast loading and rendering
- Browser may take a few seconds to render very large images
- Zooming and panning are smooth after initial load

### Feature Pairing
- Pairing is only computed when needed (on-demand)
- First network graph generation may take a few seconds
- Results are cached for subsequent views
- Changing cutoff parameters clears the cache

### Multiple Files
- Processing multiple files sequentially
- Each file's data is cached in memory
- Consider available RAM when loading many files

## Troubleshooting

### GUI won't start
- Check that all dependencies are installed: `pip install -r gui_requirements.txt`
- Verify port 8050 is not in use (or use `--port` to specify a different port)

### Processing fails
- Ensure mzML and TSV files are from the unbeQuant workflow
- Check that file names match between mzML and TSV
- Verify bin scripts are executable: `chmod +x bin/*.py`

### Heatmap not displaying
- Large heatmaps may take time to render (wait 5-10 seconds)
- Check browser console for errors
- Try a smaller dataset first

### Features not clickable
- Make sure you're clicking on the scatter points (feature centers)
- Features may be too small at low zoom levels
- Try zooming in to the feature area

### Network graph empty
- Ensure at least 2 files are loaded
- Check that pairing cutoffs are not too restrictive
- Verify features actually match across files (similar m/z and RT)

## Integration with unbeQuant Workflow

The GUI is designed to work with output from the main unbeQuant workflow:

1. **Run the main workflow** to generate:
   - mzML files: `results/mgfs_mzmls/*.mzML`
   - TSV features: `results/quantifications/features_with_annotated_identifications/*.tsv`

2. **Launch the GUI** and upload the files

3. **Explore and analyze** your results interactively

## Reference Data

Example output and reference data can be found at:
https://1drv.ms/f/c/a41a7100db890e3c/IgAoErUC6eMmRogj7_ESRU7SAZF0W9WWLgZO5CnxktwVnuE?e=KUi9PO

## Technical Details

### Architecture
- **Framework**: Dash (Plotly)
- **Backend**: Python with subprocess-based script integration
- **Frontend**: React-based Dash components
- **Visualization**: Plotly.js for interactive plots

### Script Integration
The GUI uses subprocess calls to run bin scripts, matching the nextflow workflow pattern:
- Ensures consistency with workflow processing
- Isolates script execution
- Captures stdout/stderr for debugging

### Data Flow
1. User uploads files → Saved to data directory
2. GUI calls `process_mzml_file.py` → Spectrum data pickle
3. GUI calls `create_heatmap_image_hdf5.py` → Heatmap PNG
4. GUI calls `extract_feature_data.py` → Feature data
5. (On-demand) GUI calls `pair_features.py` → Network edges
6. Visualization rendered in browser

## Future Enhancements

Potential improvements:
- [ ] Batch processing of multiple files
- [ ] Export network graphs as images
- [ ] More sophisticated diagnostic plots (3D intensity maps)
- [ ] Feature annotation and manual curation
- [ ] Integration with database search results
- [ ] Session save/load functionality
- [ ] Real-time processing progress indicators

## Contributing

This GUI is part of the unbeQuant project. For issues, suggestions, or contributions, please refer to the main project repository.

## License

Same license as the main unbeQuant project.

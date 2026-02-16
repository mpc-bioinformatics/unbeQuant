# Multi-File Support Implementation Summary

## Changes Made (Commit c0b762b)

### Overview
Enhanced the tkinter GUI to accept and process multiple mzML/TSV file pairs simultaneously, with the ability to switch between heatmaps and overlay features from different files for comparison.

### Key Features Implemented

#### 1. Multi-File Loading
**File**: `unbequant_tkinter_gui.py` - `load_files()` method

- Changed from single file selection to multiple file selection
- Uses `filedialog.askopenfilenames()` for multi-select
- Validates that mzML and TSV file counts match
- Processes all files sequentially in background thread

**User Experience**:
```
File → Load mzML + TSV Files...
→ Select multiple mzML files (Ctrl+Click)
→ Select corresponding TSV files (in same order)
→ Processing begins automatically
```

#### 2. Heatmap Selector
**File**: `unbequant_tkinter_gui.py` - `_create_layout()` method

Added control panel at top of window with:
- **Dropdown (Combobox)**: Lists all loaded heatmap files
- **Event handler**: `_on_heatmap_selected()` - switches display when selection changes
- **Auto-select**: First loaded file is automatically selected

**Layout**:
```
┌────────────────────────────────────────────────────────┐
│ Select Heatmap: [File1_sample ▼] | Feature Overlays:  │
│                                   ☐ File2  ☑ File3     │
└────────────────────────────────────────────────────────┘
```

#### 3. Feature Overlay System
**File**: `unbequant_tkinter_gui.py` - `_update_overlays()` method

- **Checkboxes**: Dynamic checkboxes for each loaded file (except current)
- **Color coding**: Each overlay gets a distinct color (blue, green, orange, purple, cyan, magenta)
- **Coordinate remapping**: Features automatically mapped from their native coordinate system to the current heatmap's system
- **Real-time updates**: Overlays update immediately when checkboxes are toggled

**Color Scheme**:
- Red: Current/selected heatmap's features
- Blue, Green, Orange, etc.: Overlay features from other files

#### 4. Processing Pipeline
**File**: `unbequant_tkinter_gui.py` - New methods

**New Methods**:
- `_process_multiple_files()`: Batch processor for multiple file pairs
- `_process_single_file()`: Processes one mzML+TSV pair
- `_update_file_selectors()`: Updates dropdown and checkboxes after loading
- `_display_current_heatmap()`: Displays selected heatmap with overlays
- `_on_overlay_changed()`: Handles checkbox state changes

**Processing Flow**:
```
User selects files
    ↓
_process_multiple_files() [background thread]
    ↓
For each file pair:
    _process_single_file()
        → process_mzml_file.py (extract spectra)
        → create_heatmap_image_hdf5.py (generate PNG)
        → extract_feature_data.py (parse features)
        → Map features to pixels
    ↓
_update_file_selectors() [main thread]
    → Update dropdown with file names
    → Create overlay checkboxes
    → Display first file automatically
```

#### 5. Coordinate System Remapping
**File**: `unbequant_tkinter_gui.py` - `_update_overlays()` method

When overlaying features from different files:
1. Get overlay file's features (with original m/z, RT values)
2. Check if overlay's coordinate system matches current heatmap's
3. If different, remap using current heatmap's `mz_dict` and `rt_dict`
4. Convert to pixel coordinates for display
5. Add to canvas with distinct color

**Example**:
```python
# Overlay feature has: mz_start=450.123, rt_start=12.34
# Current heatmap's mz_dict: {450.12: 1234, ...}
# Remap: x_min = current_mz_dict[round(450.123, 2)] = 1234
# Display overlay feature at pixel 1234
```

### Modified Files

1. **unbequant_tkinter_gui.py**
   - Added control panel with dropdown and checkboxes
   - Modified `load_files()` for multi-file selection
   - Added batch processing methods
   - Added heatmap switching logic
   - Added overlay remapping logic
   - Updated menu labels
   - Updated About dialog

2. **TKINTER_GUI_README.md**
   - Updated feature list
   - Added multi-file usage instructions
   - Added overlay instructions
   - Updated workflow section

### Technical Details

#### Data Storage
```python
self.heatmap_files = {
    'file1_sample': '/path/to/file1_heatmap.png',
    'file2_sample': '/path/to/file2_heatmap.png',
    ...
}

self.feature_data = {
    'file1_sample': [list of mapped features],
    'file2_sample': [list of mapped features],
    ...
}

self.spectrum_data = {
    'file1_sample': {
        'mz_dict': {...},
        'rt_dict': {...},
        ...
    },
    ...
}
```

#### UI State
```python
self.current_file = 'file1_sample'  # Currently displayed heatmap
self.heatmap_selector = Combobox(['file1_sample', 'file2_sample', ...])
self.overlay_checkboxes = {
    'file2_sample': BooleanVar(False),
    'file3_sample': BooleanVar(True),
    ...
}
```

### Integration with Existing Scripts

No changes required to bin scripts. The GUI uses existing scripts as-is:
- `process_mzml_file.py` - processes each mzML individually
- `create_heatmap_image_hdf5.py` - creates heatmap for each file
- `extract_feature_data.py` - extracts features for each file
- `pair_features.py` - pairs features across all loaded files (called when needed for network graph)

### User Workflow

**Before (Single File)**:
1. Load one mzML + TSV pair
2. View heatmap
3. Click features

**After (Multi File)**:
1. Load multiple mzML + TSV pairs (Ctrl+Click to select multiple)
2. Select which heatmap to display from dropdown
3. Check boxes to overlay features from other files
4. Compare features across files visually
5. Click features to see network connections across all files
6. Switch between different heatmaps as needed

### Benefits

1. **Efficiency**: Load all files once, switch between them instantly
2. **Comparison**: Visual overlay of features from multiple runs
3. **Flexibility**: Toggle overlays on/off to reduce visual clutter
4. **Accuracy**: Automatic coordinate remapping ensures proper alignment
5. **Usability**: Intuitive controls (dropdown + checkboxes)

### Testing

Syntax validation completed:
```bash
python3 -m py_compile unbequant_tkinter_gui.py
✓ Syntax check passed
```

### Future Enhancements

Potential improvements:
- Batch export of all heatmaps
- Save/load session with multiple files
- Filter overlays by feature properties
- Adjustable overlay opacity
- Overlay statistics comparison

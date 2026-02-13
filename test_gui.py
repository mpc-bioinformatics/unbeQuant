#!/usr/bin/env python3
"""
Test script for unbeQuant GUI
Validates basic functionality without requiring actual data files
"""

import sys
import tempfile
from pathlib import Path

# Test imports
print("Testing imports...")
try:
    from unbequant_gui import UnbeQuantGUI
    import dash
    import plotly
    import networkx
    print("✓ All imports successful")
except ImportError as e:
    print(f"✗ Import failed: {e}")
    sys.exit(1)

# Test GUI initialization
print("\nTesting GUI initialization...")
try:
    with tempfile.TemporaryDirectory() as tmpdir:
        gui = UnbeQuantGUI(data_dir=tmpdir)
        print("✓ GUI object created successfully")
        
        # Check that app was created
        assert gui.app is not None, "Dash app not created"
        print("✓ Dash app initialized")
        
        # Check data structures
        assert isinstance(gui.heatmap_files, dict), "heatmap_files should be dict"
        assert isinstance(gui.feature_data, dict), "feature_data should be dict"
        assert isinstance(gui.spectrum_data, dict), "spectrum_data should be dict"
        assert isinstance(gui.paired_edges, dict), "paired_edges should be dict"
        print("✓ Data structures initialized")
        
        # Check cutoff parameters
        assert gui.cutoff_mode in ['euclidean', 'coordinate'], "Invalid cutoff mode"
        assert gui.euclidean_cutoff > 0, "Invalid euclidean cutoff"
        assert gui.mz_cutoff > 0, "Invalid mz cutoff"
        assert gui.rt_cutoff > 0, "Invalid rt cutoff"
        print("✓ Cutoff parameters initialized")
        
except Exception as e:
    print(f"✗ GUI initialization failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test helper methods
print("\nTesting helper methods...")
try:
    # Test BIN_DIR exists
    from unbequant_gui import BIN_DIR
    assert BIN_DIR.exists(), f"BIN_DIR does not exist: {BIN_DIR}"
    print(f"✓ BIN_DIR found: {BIN_DIR}")
    
    # Check required scripts exist
    required_scripts = [
        'process_mzml_file.py',
        'create_heatmap_image_hdf5.py',
        'extract_feature_data.py',
        'pair_features.py'
    ]
    
    for script in required_scripts:
        script_path = BIN_DIR / script
        assert script_path.exists(), f"Required script not found: {script}"
    print(f"✓ All required scripts found ({len(required_scripts)} scripts)")
    
except Exception as e:
    print(f"✗ Helper method test failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test layout structure
print("\nTesting GUI layout...")
try:
    layout = gui.app.layout
    assert layout is not None, "Layout not created"
    print("✓ Layout created successfully")
    
except Exception as e:
    print(f"✗ Layout test failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "="*50)
print("✓ All tests passed!")
print("="*50)
print("\nGUI is ready to use. Start with:")
print("  python3 unbequant_gui.py")
print("  or")
print("  ./launch_gui.sh")

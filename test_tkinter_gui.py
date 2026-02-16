#!/usr/bin/env python3
"""
Test script for unbequant_tkinter_gui.py
Validates the GUI structure and creates mockup visualization
"""

import sys
from pathlib import Path

# Add bin to path
BIN_DIR = Path(__file__).parent / 'bin'
sys.path.insert(0, str(BIN_DIR))

def test_imports():
    """Test all required imports"""
    print("Testing imports...")
    
    try:
        import pickle
        import json
        print("  ✓ pickle, json")
    except ImportError as e:
        print(f"  ✗ Standard library imports failed: {e}")
        return False
    
    try:
        import numpy as np
        import pandas as pd
        from PIL import Image
        print("  ✓ numpy, pandas, PIL")
    except ImportError as e:
        print(f"  ✗ Data processing imports failed: {e}")
        print("    Install with: pip install numpy pandas Pillow")
        return False
    
    try:
        import networkx as nx
        print("  ✓ networkx")
    except ImportError as e:
        print(f"  ✗ networkx import failed: {e}")
        print("    Install with: pip install networkx")
        return False
    
    # tkinter and matplotlib are tested separately since they're GUI-specific
    print("\nNote: tkinter and matplotlib GUI imports will be tested when GUI runs")
    
    return True


def test_bin_scripts():
    """Test that bin scripts exist and are accessible"""
    print("\nTesting bin scripts...")
    
    scripts_to_check = [
        'process_mzml_file.py',
        'create_heatmap_image_hdf5.py',
        'extract_feature_data.py',
        'pair_features.py'
    ]
    
    all_found = True
    for script in scripts_to_check:
        script_path = BIN_DIR / script
        if script_path.exists():
            print(f"  ✓ {script}")
        else:
            print(f"  ✗ {script} not found")
            all_found = False
    
    return all_found


def test_syntax():
    """Test Python syntax of the GUI script"""
    print("\nTesting GUI script syntax...")
    
    import py_compile
    
    try:
        py_compile.compile('unbequant_tkinter_gui.py', doraise=True)
        print("  ✓ unbequant_tkinter_gui.py syntax is valid")
        return True
    except py_compile.PyCompileError as e:
        print(f"  ✗ Syntax error in unbequant_tkinter_gui.py:")
        print(f"    {e}")
        return False


def create_mockup_visualization():
    """Create a mockup showing the GUI layout"""
    print("\nCreating GUI layout mockup...")
    
    try:
        import numpy as np
        from PIL import Image, ImageDraw, ImageFont
        
        # Create a mockup image showing the 3-panel layout
        width = 1600
        height = 900
        
        img = Image.new('RGB', (width, height), color='white')
        draw = ImageDraw.Draw(img)
        
        # Define panel boundaries
        # Left/Middle panel (heatmap) - 60% width
        heatmap_width = int(width * 0.6)
        
        # Right panels - 40% width, split vertically
        right_x = heatmap_width
        right_width = width - heatmap_width
        network_height = int(height * 0.5)
        
        # Draw panel outlines
        # Heatmap panel
        draw.rectangle([0, 0, heatmap_width-2, height-2], outline='blue', width=3)
        draw.text((heatmap_width//2 - 100, 20), "HEATMAP PANEL", fill='blue')
        draw.text((heatmap_width//2 - 150, height//2), "Interactive MS1 Heatmap", fill='gray')
        draw.text((heatmap_width//2 - 150, height//2 + 30), "with Feature Overlay", fill='gray')
        draw.text((heatmap_width//2 - 100, height//2 + 60), "(Zoom & Pan)", fill='gray')
        
        # Network graph panel (top right)
        draw.rectangle([right_x, 0, width-2, network_height-2], outline='green', width=3)
        draw.text((right_x + 50, 20), "NETWORK GRAPH PANEL", fill='green')
        draw.text((right_x + right_width//2 - 80, network_height//2), "Feature Network", fill='gray')
        draw.text((right_x + right_width//2 - 100, network_height//2 + 30), "Connected Features", fill='gray')
        
        # Diagnostic panel (bottom right)
        draw.rectangle([right_x, network_height, width-2, height-2], outline='red', width=3)
        draw.text((right_x + 50, network_height + 20), "DIAGNOSTIC PANEL", fill='red')
        draw.text((right_x + right_width//2 - 80, network_height + height//4), "Feature Details", fill='gray')
        draw.text((right_x + right_width//2 - 90, network_height + height//4 + 30), "Centroid Calculation", fill='gray')
        
        # Add annotations
        draw.text((10, height - 30), "Layout: Left (60%) = Heatmap | Right (40%) = Network (top) + Diagnostic (bottom)", fill='black')
        
        # Save mockup
        output_path = Path('gui_layout_mockup.png')
        img.save(output_path)
        print(f"  ✓ Created mockup: {output_path}")
        print(f"    - Heatmap panel: {heatmap_width}x{height} px")
        print(f"    - Network panel: {right_width}x{network_height} px")
        print(f"    - Diagnostic panel: {right_width}x{network_height} px")
        
        return True
        
    except Exception as e:
        print(f"  ✗ Failed to create mockup: {e}")
        return False


def print_usage_instructions():
    """Print usage instructions"""
    print("\n" + "="*70)
    print("USAGE INSTRUCTIONS")
    print("="*70)
    print("\n1. Install dependencies:")
    print("   pip install -r gui_requirements.txt")
    print("\n2. Ensure tkinter is installed:")
    print("   - Ubuntu/Debian: sudo apt-get install python3-tk")
    print("   - macOS/Windows: Included with Python")
    print("\n3. Run the GUI:")
    print("   python3 unbequant_tkinter_gui.py")
    print("   or")
    print("   ./launch_tkinter_gui.sh")
    print("\n4. Load data:")
    print("   - File → Load mzML + TSV")
    print("   - Select mzML file from: https://1drv.ms/f/c/a41a7100db890e3c/IgBzMEItLCItRIWzNsDbiUEHATEaIxwe2ykr48puJZrPh4c?e=DfS5Vo")
    print("   - Select TSV file from: https://1drv.ms/f/c/a41a7100db890e3c/IgADiGBIEdXYTpGSyeVyzCRuAYB0sUKbbo9DRq6zqFdG9J4?e=6Fa1Xg")
    print("\n5. Interact:")
    print("   - Use matplotlib toolbar to zoom/pan heatmap")
    print("   - Click feature boxes to select")
    print("   - View network graph (top right)")
    print("   - View diagnostics (bottom right)")
    print("\n" + "="*70)


def main():
    """Main test function"""
    print("="*70)
    print("unbeQuant Tkinter GUI - Validation Test")
    print("="*70)
    
    all_tests_passed = True
    
    # Test imports
    if not test_imports():
        all_tests_passed = False
    
    # Test bin scripts
    if not test_bin_scripts():
        all_tests_passed = False
    
    # Test syntax
    if not test_syntax():
        all_tests_passed = False
    
    # Create mockup
    if not create_mockup_visualization():
        all_tests_passed = False
    
    # Print usage
    print_usage_instructions()
    
    # Summary
    print("\n" + "="*70)
    if all_tests_passed:
        print("✓ ALL TESTS PASSED")
        print("\nThe tkinter GUI is ready to use!")
        print("Note: Actual GUI testing requires a display environment with tkinter.")
    else:
        print("✗ SOME TESTS FAILED")
        print("\nPlease install missing dependencies before running the GUI.")
    print("="*70)
    
    return 0 if all_tests_passed else 1


if __name__ == "__main__":
    sys.exit(main())

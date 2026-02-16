#!/bin/bash
# Launcher script for unbeQuant Tkinter GUI

# Set up Python environment
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "Error: python3 not found. Please install Python 3.9 or higher."
    exit 1
fi

# Check for tkinter
python3 -c "import tkinter" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Error: tkinter not found."
    echo "On Ubuntu/Debian: sudo apt-get install python3-tk"
    echo "On macOS/Windows: tkinter should be included with Python"
    exit 1
fi

# Check for matplotlib
python3 -c "import matplotlib" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Error: matplotlib not found."
    echo "Install GUI requirements: pip install -r gui_requirements.txt"
    exit 1
fi

# Set matplotlib backend for tkinter
export MPLBACKEND=TkAgg

# Launch the GUI
echo "Starting unbeQuant Tkinter GUI..."
python3 unbequant_tkinter_gui.py "$@"

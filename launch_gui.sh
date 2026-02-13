#!/bin/bash
# Quick launcher for unbeQuant GUI

# Set default values
DATA_DIR="${1:-./gui_data}"
PORT="${2:-8050}"

# Create data directory if it doesn't exist
mkdir -p "$DATA_DIR"

# Check if dependencies are installed
if ! python3 -c "import dash" 2>/dev/null; then
    echo "Installing GUI dependencies..."
    pip install -q -r gui_requirements.txt
fi

# Launch GUI
echo "================================"
echo "unbeQuant Interactive GUI"
echo "================================"
echo "Data directory: $DATA_DIR"
echo "Port: $PORT"
echo "================================"
echo ""
echo "Opening GUI at http://localhost:$PORT"
echo "Press Ctrl+C to stop the server"
echo ""

# Run the GUI
python3 unbequant_gui.py --data-dir "$DATA_DIR" --port "$PORT" --debug

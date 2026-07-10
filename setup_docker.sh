
#!/bin/bash

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "Docker is not installed. Installing Docker..."
    sudo apt update
    sudo apt install -y docker.io
    sudo systemctl start docker
    sudo systemctl enable docker
    echo "Docker installed. You may need to add your user to the docker group: sudo usermod -aG docker $USER"
else
    echo "Docker is already installed."
fi

# Install Java if not already installed
if ! command -v java &> /dev/null; then
    echo "Installing Java..."
    sudo apt update
    sudo apt install -y openjdk-17-jdk
else
    echo "Java is already installed."
fi

# Install CMake if not already installed
if ! command -v cmake &> /dev/null; then
    echo "Installing CMake..."
    sudo apt install -y cmake
else
    echo "CMake is already installed."
fi

# Install pip if not already installed
if ! command -v pip &> /dev/null; then
    echo "Installing pip..."
    sudo apt install -y python3-pip
else
    echo "pip is already installed."
fi

# Install Nextflow if not already installed
if ! command -v nextflow &> /dev/null; then
    echo "Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    if [ -f nextflow ]; then
        sudo mv nextflow /usr/local/bin/
        echo "Nextflow moved to /usr/local/bin"
    fi
else
    echo "Nextflow is already installed."
fi
bash compile_and_setup_depencies.sh


git clone --recursive https://github.com/mpc-bioinformatics/unbeQuant

#!/usr/bin/env bash
set -euo pipefail

clone_if_missing() {
    local url="$1"
    local dir="$2"

    if [ -d "$dir" ] && [ -n "$(ls -A "$dir" 2>/dev/null)" ]; then
        echo "Skipping $dir (already exists and is not empty)"
    else
        rm -rf "$dir"   # remove if it exists but is empty
        echo "Cloning $url into $dir"
        git clone "$url" "$dir"
    fi
}

clone_if_missing "https://github.com/mpc-bioinformatics/xic-extractor" "xic-extractor"
clone_if_missing "https://github.com/mpc-bioinformatics/ProGFASTAGen" "ProGFASTAGen"

#!/bin/bash

set -e 

##### External Dependencies (downloadable)

### Comet v2023.01.2 Download
# Delete previously downloaded Comet
rm -rf bin/comet.linux.exe
# Download Comet
wget -O bin/comet.linux.exe https://github.com/UWPR/Comet/releases/download/v2023.01.2/comet.linux.exe 


### ThermoRawFileParser 1.4.2 Download and extraction
# Delete previously downloaded ThermoRawFileParser
rm -rf bin/ThermoRawFileParser
# Download ThermoRawFileParser
wget -O trfp.zip https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.2/ThermoRawFileParser1.4.2.zip
# Extract archive
unzip trfp.zip -d bin/ThermoRawFileParser
# Delete downloaded archive
rm trfp.zip


### Percolator v3.0.6 Download and extraction
# Delete previously downloaded OPercolatorpenMS
rm -rf bin/percolator/
# Download percolator
wget -O percolator.deb https://github.com/percolator/percolator/releases/download/rel-3-06/percolator-v3-06-linux-amd64.deb
# Extract only the data-part of the deb package
ar x percolator.deb data.tar.gz
# Extract the data into the percolator-Folder
mkdir -p ./bin/percolator
tar -xf data.tar.gz -C ./bin/percolator
# Delete downloaded and extracted deb package
rm -rf percolator.deb data.tar.gz


### OpenMS Dependencies
# Delete previously downloaded openms
rm -rf bin/openms/
# Download openms
wget -O openms.deb https://github.com/OpenMS/OpenMS/releases/download/Release2.8.0/OpenMS-2.8.0-Debian-Linux-x86_64.deb
# Extract only the data-part of the deb package
ar x openms.deb data.tar.gz
# Extract the data into the openms_folder
mkdir -p ./bin/openms
# Delete downloaded and extracted deb package
tar -xf data.tar.gz -C ./bin/openms
rm -rf openms.deb data.tar.gz


### Python Dependencies (Python 3.9)
# Protein-Graph-Generation and other exports
pip install protgraph==0.3.10
# Needed for exporting SQLite via ProtGraph
pip install apsw==3.42.0.0
# Needed to generate compact Excel-Files (HeatMaps)
pip install XlsxWriter==3.0.3
# Used to load in comet-txt-files
pip install pandas
# Needed for generating various plots
pip install plotly
# Engine, needed to run all workflows
pip install nextflow
# Needed for parsing xmls (e.g. Percolator output)
pip install lxml
# For various things in python scripts
pip install pyopenms==3.1.0
# Exporting figures
pip install kaleido


##### Make all files in bin executable (excluding sub-directories) to be visible by processes in Nextflow
chmod +x bin/*

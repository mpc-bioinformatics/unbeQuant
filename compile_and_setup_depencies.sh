#!/bin/bash

set -e 

##### External Dependencies (downloadable)

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
# Needed for exporting SQLite via ProtGraph
pip install apsw==3.42.0.0
# Needed to generate compact Excel-Files (HeatMaps)
pip install XlsxWriter==3.0.3
# Used to load in comet-txt-files
pip install pandas
# tqdm, need by some python scripts
pip install tqdm
# Needed for generating various plots
pip install plotly
# Engine, needed to run all workflows
pip install nextflow
# Needed for parsing xmls (e.g. Percolator output)
pip install lxml
# For various things in python scripts
pip install numpy==1.26.4
pip install pyopenms==3.1.0
# Extracting XICs from hdf5
pip install h5py
# Exporting figures
pip install kaleido==0.2.1
# Install biosaur2
pip install biosaur2
# XLSX (final) reports 
pip install xlsxwriter


##### Make all files in bin executable (excluding sub-directories) to be visible by processes in Nextflow
chmod +x bin/*

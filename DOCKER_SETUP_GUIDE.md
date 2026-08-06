# Custom Docker Container for map_mzml_features

This guide explains how to build and use the custom Docker container for the `map_mzml_features_v2.nf` subworkflow.

## Overview

A dedicated Docker container has been created to isolate the map_mzml_features workflow with all its Python dependencies (scipy, numpy, pandas, python-igraph, matplotlib, plotly, h5py, etc.). This allows:

- **Isolation**: map_mzml_features runs in its own container with all dependencies
- **Portability**: Works across different systems without manual dependency installation
- **Integration**: Works seamlessly with the main `main_unbequant.nf` workflow
- **Flexibility**: Can still run with Docker disabled if scripts are run directly on host

## Building the Docker Image

### Prerequisites
- Docker installed and running
- ~1.5 GB disk space for the image

### Build Command

Run the provided build script:

```bash
bash build_docker_image.sh latest
```

Or build manually:

```bash
docker build -f Dockerfile.map_mzml_features -t unbequant/map_mzml_features:latest .
```

To build with a specific version tag:

```bash
bash build_docker_image.sh v1.0
docker build -f Dockerfile.map_mzml_features -t unbequant/map_mzml_features:v1.0 .
```

### Verify Build

```bash
docker images | grep map_mzml_features
```

You should see:
```
unbequant/map_mzml_features   latest    <IMAGE_ID>   <created>   <SIZE>
```

## Using the Container

### Option 1: With main_unbequant.nf (Recommended)

The container is automatically used for map_mzml_features processes when Docker is enabled:

```bash
# Enable Docker in nextflow.config
docker.enabled = true

# Run main workflow
nextflow run main_unbequant.nf
```

The main workflow will:
- Use ProGFASTAGen containers for identification and conversion
- Use the custom map_mzml_features container for feature mapping processes

### Option 2: Standalone map_mzml_features_v2.nf with Container

To run map_mzml_features as a standalone workflow with Docker:

```bash
# Update nextflow_standalone.config to enable Docker
# Change: docker.enabled = false → docker.enabled = true

nextflow run map_mzml_features_v2.nf -c nextflow_standalone.config \
    --mzml_files 'results/mgfs_mzmls/*.mzML' \
    --tsv_files 'results/quantifications/features_with_annotated_identifications/*.tsv'
```

### Option 3: Standalone without Docker (Host Scripts)

To run without Docker, using Python scripts directly on host:

```bash
# Keep Docker disabled in nextflow_standalone.config
docker.enabled = false

# Ensure all dependencies are installed on host
pip install numpy scipy pandas matplotlib plotly python-igraph h5py psutil kaleido

# Run workflow
nextflow run map_mzml_features_v2.nf -c nextflow_standalone.config \
    --mzml_files 'results/mgfs_mzmls/*.mzML' \
    --tsv_files 'results/quantifications/features_with_annotated_identifications/*.tsv'
```

## Configuration Details

### Process Labels

The map_mzml_features processes use these labels (defined in nextflow.config):

- `map_mzml_features` - General processes (extract, report, consensus conversion)
  - Memory: 8G, CPUs: 2

- `map_mzml_features_heatmap` - Heatmap generation process
  - Memory: 16G, CPUs: 2

- `map_mzml_features_pairing` - Feature pairing process
  - Memory: 12G, CPUs: 4

- `map_mzml_features_graph` - Network graph visualization
  - Memory: 12G, CPUs: 15

### Container Image Reference

```
unbequant/map_mzml_features:latest
```

To use a different version:

```
unbequant/map_mzml_features:v1.0
```

Update in:
- `nextflow.config` (line ~47-58)
- `nextflow_standalone.config` (line ~35-51)

## Installed Dependencies

The Docker image includes:

### System Packages
- build-essential (C/C++ compilation)
- graphviz (graph visualization)
- libigraph-dev (igraph library)
- git

### Python Packages
- numpy (numerical computing)
- scipy (scientific computing, KD-tree, optimization)
- pandas (data analysis)
- matplotlib (visualization)
- plotly (interactive plots)
- python-igraph (network analysis)
- h5py (HDF5 file handling)
- psutil (system monitoring)
- kaleido (static image export for plotly)

## Troubleshooting

### Image Not Found

If you get "docker: image not found" error:

1. Check image exists:
   ```bash
   docker images | grep map_mzml_features
   ```

2. Build the image:
   ```bash
   bash build_docker_image.sh latest
   ```

3. Verify Nextflow can access Docker:
   ```bash
   docker run --rm unbequant/map_mzml_features:latest echo "Docker is working"
   ```

### Memory Issues

If processes run out of memory:

1. Check current limits in nextflow.config (lines ~47-58)
2. Increase memory for specific label:
   ```nextflow
   withLabel: map_mzml_features_pairing {
       memory = '24G'  // Increased from 12G
   }
   ```

3. Reduce queue size to limit parallel processes:
   ```nextflow
   executor {
       queueSize = 2  // Reduce to 2 or 3
   }
   ```

### Permission Issues

If Docker daemon access fails:

```bash
# Add user to docker group (requires sudo)
sudo usermod -aG docker $USER

# Apply new group membership
newgrp docker

# Or restart docker service
sudo systemctl restart docker
```

## Integration with Main Workflow

The custom container integrates with main_unbequant.nf as follows:

```
main_unbequant.nf
├── ProGFASTAGen workflows
│   ├── convert_to_mgf       (uses: luxii/progfastagen:1.0)
│   ├── convert_to_mzml      (uses: thermorawfileparser & bruker containers)
│   ├── identification_via_comet (uses: luxii/progfastagen:1.0)
│   └── summarize_ident_results  (uses: luxii/progfastagen:1.0)
├── quantify_and_align       (uses: local OpenMS)
└── map_mzml_features        (uses: unbequant/map_mzml_features:latest) ← NEW
    ├── process_mzml_file
    ├── create_heatmap_image
    ├── extract_feature_data
    ├── pair_features
    ├── feature_pairing_report
    ├── build_network_graph_visualization
    ├── compare_rt_alignment
    ├── convert_clusters_to_consensusxml
    └── process_clusters_consensusxml
```

## Advanced: Customizing the Image

To add additional dependencies, edit `Dockerfile.map_mzml_features`:

```dockerfile
# Add system package
RUN apt-get update && apt-get install -y \
    my-new-package

# Add Python package
RUN pip install --no-cache-dir \
    my-python-package
```

Then rebuild:

```bash
docker build -f Dockerfile.map_mzml_features -t unbequant/map_mzml_features:custom .
```

Update configuration to use the new image:

```nextflow
withLabel: map_mzml_features {
    container = "unbequant/map_mzml_features:custom"
}
```

## Resources

- [Docker Documentation](https://docs.docker.com)
- [Nextflow Docker Integration](https://www.nextflow.io/docs/latest/container.html)
- [Python Dependencies](https://pypi.org)

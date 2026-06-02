#!/bin/bash
# Build Docker image for map_mzml_features_v2.nf workflow

set -e

IMAGE_NAME="unbequant/map_mzml_features"
IMAGE_TAG="${1:-latest}"
FULL_IMAGE_NAME="${IMAGE_NAME}:${IMAGE_TAG}"

echo "Building Docker image: ${FULL_IMAGE_NAME}"
echo "This may take a few minutes..."

docker build \
    -f Dockerfile.map_mzml_features \
    -t "${FULL_IMAGE_NAME}" \
    .

echo ""
echo "✓ Docker image built successfully!"
echo ""
echo "Image name: ${FULL_IMAGE_NAME}"
echo ""
echo "To use this image in Nextflow, update nextflow.config with:"
echo "  withLabel: map_mzml_features { container = \"${FULL_IMAGE_NAME}\" }"
echo ""
echo "Or tag processes in map_mzml_features_v2.nf with:"
echo "  label \"map_mzml_features\""

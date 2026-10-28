#!/bin/bash
# Build and upload container with version from pixi.toml (run in WSL; images
# cannot be built on DAIC).
#
#   ./build_and_upload.sh               build, then upload
#   ./build_and_upload.sh --skip-build  upload an existing .sif

set -e

SKIP_BUILD=0
[ "${1:-}" = "--skip-build" ] && SKIP_BUILD=1

# Extract version from pixi.toml
VERSION=$(grep '^version = ' pixi.toml | cut -d'"' -f2 | tr '.' '-')

if [ -z "$VERSION" ]; then
    echo "Error: Could not extract version from pixi.toml"
    exit 1
fi

CONTAINER_NAME="container_pixi_${VERSION}.sif"

echo "Building container version: ${VERSION}"
echo "Container file: ${CONTAINER_NAME}"
echo ""

# Build the container
if [ "$SKIP_BUILD" = 1 ]; then
    [ -f "${CONTAINER_NAME}" ] || { echo "Error: ${CONTAINER_NAME} not found"; exit 1; }
    echo "Skipping build, using existing ${CONTAINER_NAME}"
else
    apptainer build "${CONTAINER_NAME}" container_pixi.def
fi

echo ""
echo "Uploading to DAIC..."
echo ""

# Upload to DAIC
rsync --progress -vhz "${CONTAINER_NAME}" "daic:/tudelft.net/staff-umbrella/scGraphNN/rosmap-processing/${CONTAINER_NAME}"

echo ""
echo "Done! Container ${CONTAINER_NAME} uploaded successfully."

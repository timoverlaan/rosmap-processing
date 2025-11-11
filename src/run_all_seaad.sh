#!/bin/bash
# This script downloads the SeaAD dataset from the Amazon S3 bucket and prepares it for analysis.

# Exit on any error
set -e

echo "========================================="
echo "SeaAD Processing Pipeline"
echo "========================================="
echo ""

# Set up run directories and save config/git info
echo "Setting up run directories..."
RUN_INFO=$(pixi run python src/setup_run.py --config config.yaml --run-name seaad_processing)
eval "$RUN_INFO"
echo "Output directory: ${OUTPUT_DIR}"
echo "Log directory: ${LOG_DIR}"
echo ""

# Data folder (if it doesn't exist yet)
mkdir -p data/seaAD/PFC/RNAseq

# Download SeaAD from the amazon s3 bucket
echo "[Step 1/3] Downloading SeaAD data from S3..."
pixi run aws s3 cp --no-sign-request \
    s3://sea-ad-single-cell-profiling/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad \
    data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad

# Download the metadata file
pixi run aws s3 cp --no-sign-request \
    s3://sea-ad-single-cell-profiling/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv \
    data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv
echo "✓ Download complete"
echo ""

# Check the file (print a bunch of stuff)
# pixi run python -m rosmap_processing utils validate \
#     data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad

echo "[Step 2/3] Converting column format..."
pixi run python -m rosmap_processing core column-mapping \
    data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad --inplace --data-type SeaAD
echo "✓ Column format conversion complete"
echo ""

echo "[Step 3/3] Running scanpy pipeline..."
pixi run python -m rosmap_processing pipeline scanpy \
    data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad \
    ${OUTPUT_DIR}/SeaAD_2k_k30.h5ad \
    --layer UMIs --n-genes 2000 --k-neighbors 30
echo "✓ Scanpy pipeline complete"
echo ""

echo "========================================="
echo "Pipeline complete!"
echo "========================================="
echo "Output saved to: ${OUTPUT_DIR}/SeaAD_2k_k30.h5ad"

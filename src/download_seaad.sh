#!/bin/bash
# This script downloads the SeaAD dataset from the Amazon S3 bucket and prepares it for analysis.

# Data folder (if it doesn't exist yet)
mkdir -p data/seaAD/PFC/RNAseq

# Download SeaAD from the amazon s3 bucket
aws s3 cp --no-sign-request \
    s3://sea-ad-single-cell-profiling/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad \
    data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad

# Download the metadata file
aws s3 cp --no-sign-request \
    s3://sea-ad-single-cell-profiling/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv \
    data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv 

# Check the file (print a bunch of stuff)
python -u src/check_h5ad.py \
    data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad
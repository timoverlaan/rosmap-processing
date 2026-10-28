#!/bin/sh
# SeaAD (PFC) setup: download from the public S3 bucket, then add
# ROSMAP-style columns (braaksc, ceradsc, thalsc, Wang labels) into a copy in
# data/v2/counts/SeaAD/. The download itself is never modified. The counts are
# in the 'UMIs' layer; the processing itself is done by the slurm/*hvg*.sh jobs.
# Runs inside the container (submit slurm/run_all_seaad.sh).

set -eu

V2=data/v2
SRC_DIR=data/seaAD/PFC/RNAseq
H5AD=SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad
META=SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv

echo "========================================="
echo "SeaAD Processing Pipeline (v2)"
echo "========================================="

RUN_INFO=$(pixi run python src/setup_run.py --config config.yaml --run-name seaad_processing)
eval "$RUN_INFO"
echo "Run directory: ${OUTPUT_DIR}"

mkdir -p "$SRC_DIR" "$V2/counts/SeaAD"

echo "[Step 1/2] Downloading SeaAD data from S3 (skips existing files)..."
for f in "$H5AD" "$META"; do
    if [ -e "$SRC_DIR/$f" ]; then
        echo "exists, skipping: $SRC_DIR/$f"
    else
        pixi run aws s3 cp --no-sign-request "s3://sea-ad-single-cell-profiling/PFC/RNAseq/$f" "$SRC_DIR/$f"
    fi
done

echo "[Step 2/2] Adding ROSMAP-style columns and Wang labels (copy in $V2/counts/SeaAD)..."
if [ -e "$V2/counts/SeaAD/$H5AD" ]; then
    echo "exists, skipping: $V2/counts/SeaAD/$H5AD"
else
    pixi run python -m rosmap_processing core column-mapping "$SRC_DIR/$H5AD" \
        --data-type SeaAD --output "$V2/counts/SeaAD/$H5AD"
fi

echo "========================================="
echo "Done. SeaAD counts (layer 'UMIs'): $V2/counts/SeaAD/$H5AD"
echo "========================================="

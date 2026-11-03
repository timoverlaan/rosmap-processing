#!/bin/sh
# PsychAD: pinned download -> data/v2/counts/PsychAD/combined.h5ad, the counts
# of merged_final_clean.h5ad (layers['counts']) as int32 in X, with donor
# metadata, diagnoses, hashing set and channel in obs. Runs inside the
# container (submit slurm/run_all_psychad.sh).
#
# Inputs: data/v2/download/PsychAD (from `data download --dataset psychad`)
# and the ROSMAP clinical CSV (from run_all.sh) for the Rush donors.
# The step skips an output that already exists.

set -eu

V2=data/v2
DL=$V2/download/PsychAD
OUT=$V2/counts/PsychAD
CLINICAL=$V2/download/ROSMAP/ROSMAP_clinical.csv

rp() { pixi run python -m rosmap_processing "$@"; }

echo "========================================="
echo "PsychAD Processing Pipeline (v2)"
echo "========================================="

echo "[Step 1/2] Verifying the pinned downloads..."
rp data verify --dataset psychad --dir "$DL"

echo "[Step 2/2] Counts + donor metadata, streamed into one file..."
mkdir -p "$OUT"
if [ -e "$OUT/combined.h5ad" ]; then
    echo "exists, skipping: $OUT/combined.h5ad"
else
    rp data psychad "$DL/merged_final_clean.h5ad" --output "$OUT/combined.h5ad" --download-dir "$DL" \
        --clinical "$CLINICAL"
fi

echo "========================================="
echo "Done. Counts in $OUT/combined.h5ad"
echo "========================================="

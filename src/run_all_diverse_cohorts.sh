#!/bin/sh
# DiverseCohorts (AMP-AD) 10x multiome, gene expression: pinned download ->
# DLPFC count files in data/v2/counts/DiverseCohorts/. Runs inside the
# container (submit slurm/run_all_diverse_cohorts.sh).
#
# Inputs (data/v2/download/DiverseCohorts, from `data download --dataset
# diverse-cohorts`): gex_mat_<class>.rda (2025 release), the cell table
# rosmap_diversity_snmultiome_clusterdata_20241107.csv (donor, region and
# cell class of every cell; the libraries pool three donors), and the donor
# metadata. The ROSMAP clinical CSV (from run_all.sh) adds the ROSMAP columns.
# Outputs:
#   data/v2/intermediate/DiverseCohorts/gex_mat_<class>.h5   the R matrices as HDF5
#   data/v2/counts/DiverseCohorts/<class>.h5ad              DLPFC cells, int32 counts in X
#   data/v2/counts/DiverseCohorts/combined.h5ad             all classes
# Steps skip outputs that already exist, so a failed run can be resubmitted.

set -eu

V2=data/v2
DL=$V2/download/DiverseCohorts
MID=$V2/intermediate/DiverseCohorts
OUT=$V2/counts/DiverseCohorts
CLINICAL=$V2/download/ROSMAP/ROSMAP_clinical.csv
CLASSES="Astro Endo Epend Gad1 Gad2 Glut1 Glut2 Glut3 Glut4 Micro Olig1 Olig2 Opc"

rp() { pixi run python -m rosmap_processing "$@"; }

echo "========================================="
echo "DiverseCohorts Processing Pipeline (v2)"
echo "========================================="

echo "[Step 1/4] Verifying the pinned downloads..."
rp data verify --dataset diverse-cohorts --dir "$DL"

echo "[Step 2/4] gex_mat_*.rda -> HDF5 (skips existing)..."
mkdir -p "$MID" "$OUT"
for c in $CLASSES; do
    pixi run Rscript src/convert_R/rdata_to_h5.R "$DL/gex_mat_$c.rda" "$MID/gex_mat_$c.h5"
done

echo "[Step 3/4] Counts + cell and donor metadata, DLPFC cells only..."
for c in $CLASSES; do
    out=$OUT/$c.h5ad
    if [ -e "$out" ]; then echo "exists, skipping: $out"; continue; fi
    rp data diverse-cohorts "$MID/gex_mat_$c.h5" --output "$out" --download-dir "$DL" --clinical "$CLINICAL"
done

echo "[Step 4/4] Combining the classes..."
out=$OUT/combined.h5ad
if [ -e "$out" ]; then
    echo "exists, skipping: $out"
else
    inputs=""
    for c in $CLASSES; do inputs="$inputs $OUT/$c.h5ad"; done
    # shellcheck disable=SC2086
    rp core combine $inputs --output "$out"
fi

echo "========================================="
echo "Done. Counts in $OUT/"
echo "========================================="

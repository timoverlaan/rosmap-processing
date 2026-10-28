#!/bin/sh
# ROSMAP + ROSMAP-MIT pipeline: Synapse download -> raw-count h5ad files.
# Runs inside the container (submit slurm/run_all.sh).
#
# Inputs (all pinned to Synapse versions and checked by MD5, see SYNAPSE_FILES):
#   data/v2/download/ROSMAP/   Green 2024 objects + ROSMAP metadata, downloaded fresh
#                              (the copies in data/raw/ROSMAP were rewritten in 2025)
#   data/raw/ROSMAP_MIT/       Mathys 2023 RDS + metadata, verified in place; the
#                              .h5Seurat files made from the RDS are written here too
#                              (set MIT_RAW to use another directory)
# Outputs in data/v2/:
#   data/v2/counts/<DATASET>/<class>.h5ad     raw counts in X (int32), clinical metadata in obs
#                                             (RNA counts read from the .h5Seurat in Python)
#   data/v2/counts/<DATASET>/combined.h5ad    all classes of a dataset
# Steps skip outputs that already exist, so a failed run can be resubmitted.

set -eu

V2=data/v2
GREEN=$V2/download/ROSMAP
MIT=${MIT_RAW:-data/raw/ROSMAP_MIT}
CLINICAL=$GREEN/ROSMAP_clinical.csv
MIT_CLASSES="Astrocytes Excitatory_neurons_set1 Excitatory_neurons_set2 Excitatory_neurons_set3 Immune_cells Inhibitory_neurons OPCs Oligodendrocytes"
ROSMAP_CLASSES="astrocytes cux2+ cux2- inhibitory microglia oligodendroglia vascular.niche"

rp() { pixi run python -m rosmap_processing "$@"; }

echo "========================================="
echo "ROSMAP Processing Pipeline (v2)"
echo "========================================="

# Record config and git info for this run
RUN_INFO=$(pixi run python src/setup_run.py --config config.yaml --run-name rosmap_processing)
eval "$RUN_INFO"
echo "Run directory: ${OUTPUT_DIR}"

echo "[Step 1/4] Downloading / verifying pinned Synapse files..."
# Green objects + ROSMAP metadata into $GREEN (present files are MD5-checked, not re-downloaded)
rp data download --dataset rosmap --output-dir $V2/download
# Mathys 2023 files: checked in place; fails with a list of mismatches if any differ
rp data verify --dataset rosmap-mit --dir "$MIT"

echo "[Step 2/4] Converting ROSMAP-MIT RDS files to h5Seurat (skips existing)..."
# The April 2026 .h5Seurat files exist, so this normally does nothing. SeuratDisk
# may fail with the container's SeuratObject (`GetAssayData(slot = )` is defunct).
RDS=""
for c in $MIT_CLASSES; do RDS="$RDS $MIT/$c.rds"; done
# shellcheck disable=SC2086
pixi run Rscript src/convert_R/rds_to_h5Seurat.R $RDS

echo "[Step 3/4] RNA counts from the .h5Seurat into X, checked over the whole matrix, plus clinical metadata..."
for c in $MIT_CLASSES; do
    out=$V2/counts/ROSMAP_MIT/$c.h5ad
    if [ -e "$out" ]; then echo "exists, skipping: $out"; continue; fi
    rp data prepare $MIT/$c.h5Seurat --output "$out" \
        --columns cell_type_high_resolution --metadata $CLINICAL
done
for c in $ROSMAP_CLASSES; do
    out=$V2/counts/ROSMAP/$c.h5ad
    if [ -e "$out" ]; then echo "exists, skipping: $out"; continue; fi
    rp data prepare "$GREEN/$c.h5Seurat" --output "$out" --columns --metadata $CLINICAL
done

echo "[Step 4/4] Combining per-class files..."
for ds in ROSMAP_MIT ROSMAP; do
    if [ "$ds" = ROSMAP_MIT ]; then classes=$MIT_CLASSES; else classes=$ROSMAP_CLASSES; fi
    out=$V2/counts/$ds/combined.h5ad
    if [ -e "$out" ]; then echo "exists, skipping: $out"; continue; fi
    inputs=""
    for c in $classes; do inputs="$inputs $V2/counts/$ds/$c.h5ad"; done
    # shellcheck disable=SC2086
    rp core combine $inputs --output "$out"
done

echo "========================================="
echo "Pipeline complete. Counts in $V2/counts/"
echo "========================================="

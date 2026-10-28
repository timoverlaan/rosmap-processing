#!/bin/bash
# Green 2024 ROSMAP atlas on its own: per-donor HVGs (seurat_v3 on counts, and
# seurat after normalization) for 1k/2k/3k genes -> rosmap_green_hvg*_k30.h5ad.
# These have their own genes; for outputs on the same genes as ROSMAP-MIT and
# SeaAD, see aligned_hvg.sh. Donors are identified by individualID.
# Note: Green 2024 and Mathys 2023 (ROSMAP-MIT) largely share participants.
# Needs data/v2/counts/ROSMAP (slurm/run_all.sh) and data/v2/reference (slurm/gene_maps.sh).
# Submit from the repo root: sbatch slurm/green_hvg.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=500GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_green_hvg.out
#SBATCH --error=slurm/out/%j_green_hvg.out

source slurm/common.sh

for N in 1000 2000 3000; do
    process_green "hvg${N}_k30" --n-genes "$N"
    process_green "hvg${N}_postnorm_k30" --hvg-after-normalize "$N"
done

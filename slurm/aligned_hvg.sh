#!/bin/bash
# ROSMAP-MIT, SeaAD and Green 2024 on exactly the same 2000 genes: per-donor
# HVGs (seurat_v3 on counts) chosen on ROSMAP-MIT among the genes all three
# datasets have, then SeaAD and Green aligned to them ->
#   rosmap_mit_hvg2000_3way_k30.h5ad, seaad_..., rosmap_green_...
# Needs data/v2/counts (slurm/run_all.sh) and data/v2/reference (slurm/gene_maps.sh).
# Submit from the repo root: sbatch slurm/aligned_hvg.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_aligned_hvg.out
#SBATCH --error=slurm/out/%j_aligned_hvg.out

source slurm/common.sh

process_triplet "hvg2000_3way_k30" --n-genes 2000

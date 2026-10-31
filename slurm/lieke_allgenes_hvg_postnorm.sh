#!/bin/bash
# 2000 HVGs after normalization (seurat flavor, per donor) within Lieke's
# all-genes list, then SeaAD with the same genes. Compare with
# lieke_allgenes_hvg.sh.
# Needs data/v2/counts (slurm/run_all.sh) and data/v2/reference (slurm/gene_maps.sh).
# Submit from the repo root: sbatch slurm/lieke_allgenes_hvg_postnorm.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_lieke_allgenes_hvg_postnorm.out
#SBATCH --error=slurm/out/%j_lieke_allgenes_hvg_postnorm.out

source slurm/common.sh

process_pair "lieke_allgenes_hvg2k_postnorm_k30" \
    --import-genes data/lieke_allgenes_intersection.txt --hvg-after-normalize 2000

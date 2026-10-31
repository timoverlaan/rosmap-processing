#!/bin/bash
# Lieke's HVG list (intersected with ROSMAP-MIT x SeaAD), then SeaAD with the
# same genes.
# Needs data/v2/counts (slurm/run_all.sh) and data/v2/reference (slurm/gene_maps.sh).
# Submit from the repo root: sbatch slurm/lieke_hvg.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_lieke_hvg.out
#SBATCH --error=slurm/out/%j_lieke_hvg.out

source slurm/common.sh

process_pair "lieke_hvg_k30" --import-genes data/lieke_hvg_intersection.txt

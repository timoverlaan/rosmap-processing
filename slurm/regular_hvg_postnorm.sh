#!/bin/bash
# HVG selection after normalization (seurat flavor, per donor) for 1k/2k/3k
# genes, then SeaAD with the same genes. Compare with regular_hvg.sh.
# Needs data/v2/counts (slurm/run_all.sh) and data/v2/reference (slurm/gene_maps.sh).
# Submit from the repo root: sbatch slurm/regular_hvg_postnorm.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_regular_hvg_postnorm.out
#SBATCH --error=slurm/out/%j_regular_hvg_postnorm.out

source slurm/common.sh

for N in 1000 2000 3000; do
    process_pair "hvg${N}_postnorm_k30" --hvg-after-normalize "$N"
done

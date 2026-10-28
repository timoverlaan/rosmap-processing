#!/bin/bash
# Mathys 2019 DEG list, then SeaAD with the same genes.
# Needs data/v2/counts (slurm/run_all.sh) and data/v2/reference (slurm/gene_maps.sh).
# Submit from the repo root: sbatch slurm/rosmap_mathys_genes.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=500GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_mathys_genes.out
#SBATCH --error=slurm/out/%j_mathys_genes.out

source slurm/common.sh

process_pair "mathys2019degs_k30" --import-genes data/mathys2019_DEGs_genes.txt

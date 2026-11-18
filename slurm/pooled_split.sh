#!/bin/bash
# What splits one dataset's pooled pseudobulk samples into two groups
# (scripts/pooled_split.py; needs slurm/pooled_pseudobulk.sh). Aggregate
# output only: data/v2/pooled_pseudobulk/split_<dataset>/report.txt, umap.png.
#   sbatch slurm/pooled_split.sh                       # DiverseCohorts
#   sbatch slurm/pooled_split.sh --dataset rosmap_green --pool-column ... --class-column ...
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48GB
#SBATCH --output=slurm/out/%j_pooled_split.out
#SBATCH --error=slurm/out/%j_pooled_split.out

THREADS="${SLURM_CPUS_PER_TASK:-4}"
EXTRA_BINDS="--bind config/:/opt/app/config/ --bind scripts/:/opt/app/scripts/ --env OMP_NUM_THREADS=$THREADS --env OPENBLAS_NUM_THREADS=$THREADS --env MKL_NUM_THREADS=$THREADS"
source slurm/common.sh

in_container pixi run python scripts/pooled_split.py config/pooled_pseudobulk.yaml "$@"

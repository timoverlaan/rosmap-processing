#!/bin/bash
# Check that the v2 count files hold raw counts over the whole matrix. Prints
# aggregate statistics only (shape, dtype, number of non-integer values,
# min/max), never values or IDs.
# Submit from the repo root: sbatch slurm/inspect_counts.sh [files...]
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --output=slurm/out/%j_inspect_counts.out
#SBATCH --error=slurm/out/%j_inspect_counts.out

source slurm/common.sh

if [ "$#" -gt 0 ]; then
    FILES=("$@")
else
    FILES=("$V2"/counts/ROSMAP_MIT/*.h5ad "$V2"/counts/ROSMAP/*.h5ad)
fi
rp utils check-counts "${FILES[@]}"

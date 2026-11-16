#!/bin/bash
# Data QC, step 1: one pass over a count file (per-cell values and pseudobulk
# sums -> the spec's output_dir/restricted/, kept on the cluster).
# Submit from the repo root (slurm/submit_qc.sh submits all five):
#   sbatch slurm/qc_scan.sh config/data_qc/rosmap_mit.yaml [--overwrite]
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_qc_scan.out
#SBATCH --error=slurm/out/%j_qc_scan.out

SPEC="${1:?usage: sbatch slurm/qc_scan.sh config/data_qc/<name>.yaml [--overwrite]}"
shift
THREADS="${SLURM_CPUS_PER_TASK:-4}"
EXTRA_BINDS="--bind config/:/opt/app/config/ --env OMP_NUM_THREADS=$THREADS --env OPENBLAS_NUM_THREADS=$THREADS --env MKL_NUM_THREADS=$THREADS"
source slurm/common.sh

echo "spec: $SPEC $*"
rp pipeline qc-scan "$SPEC" "$@"

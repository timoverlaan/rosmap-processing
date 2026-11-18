#!/bin/bash
# Data QC, step 2: one check of `pipeline qc-check` (aggregate tables) for
# each given spec, or once over all of them for the checks that compare
# datasets (sex-de-compare, identity -> data/v2/qc/<check>/).
# A failure for one dataset does not stop the others (the job then fails).
#   sbatch slurm/qc_check.sh sex config/data_qc/*.yaml
#   sbatch slurm/qc_check.sh identity config/data_qc/{rosmap_mit,rosmap_green,diverse_cohorts,psychad}.yaml
#   sbatch slurm/qc_check.sh genotype config/data_qc/{rosmap_mit,rosmap_green,diverse_cohorts,psychad}.yaml
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_qc_check.out
#SBATCH --error=slurm/out/%j_qc_check.out

CHECK="${1:?usage: sbatch slurm/qc_check.sh CHECK SPEC [SPEC ...]}"
shift
THREADS="${SLURM_CPUS_PER_TASK:-4}"
EXTRA_BINDS="--bind config/:/opt/app/config/ --env OMP_NUM_THREADS=$THREADS --env OPENBLAS_NUM_THREADS=$THREADS --env MKL_NUM_THREADS=$THREADS"
source slurm/common.sh

echo "check: $CHECK on $*"
case "$CHECK" in
    sex-de-compare|identity)
        rp pipeline qc-check "$CHECK" "$@"
        ;;
    genotype)
        # Dosages at the eQTL variants (scripts/extract_dosages.py, restricted)
        # and the GTEx cortex eQTLs (scripts/gtex_eqtls.py)
        rp pipeline qc-check genotype "$@" \
            --genotypes data/v2/genotypes/restricted/rosmap_imputed_eqtl_dosages.tsv.gz \
            --eqtls data/v2/reference/gtex_cortex_eqtls.tsv
        ;;
    *)
        status=0
        for spec in "$@"; do
            echo "======================================== $CHECK: $spec"
            if [ "$CHECK" = ad-de ]; then
                rp pipeline qc-check "$CHECK" "$spec" --mathys2019 data/mathys2019_DEGs.csv || status=1
            else
                rp pipeline qc-check "$CHECK" "$spec" || status=1
            fi
        done
        exit "$status"
        ;;
esac

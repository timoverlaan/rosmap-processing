#!/bin/bash
# The genotype check against WGS, after slurm/fetch_wgs.sh amp1 and divco:
#   1. scripts/genotype_wgs_check.py per WGS set: array vs WGS per person, DNA
#      sex (chrX) vs the record, the flagged donors' genotype entries, and the
#      WGS dosages renamed by individualID;
#   2. qc-check genotype with the WGS instead of the array: MIT, Green and
#      DiverseCohorts against AMP-AD 1.0 WGS (GRCh37), DiverseCohorts against
#      its own WGS (GRCh38).
# Everything goes to data/v2/qc/genotype_wgs/. Submit from the repo root:
#   sbatch slurm/genotype_wgs.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=slurm/out/%j_genotype_wgs.out
#SBATCH --error=slurm/out/%j_genotype_wgs.out

EXTRA_BINDS="--bind scripts/:/opt/app/scripts/"
source slurm/common.sh

G="$V2/genotypes/restricted"
OUTDIR="$V2/qc/genotype_wgs"
D="$V2/download"
for name in amp1 divco; do
    build=37
    [ "$name" = divco ] && build=38
    echo "======================================== array vs WGS: $name"
    in_container pixi run python -u scripts/genotype_wgs_check.py \
        --array "$G/rosmap_imputed_eqtl_dosages.tsv.gz" --array-chrx "$G/rosmap_array_chrx_het.tsv" \
        --wgs "$G/wgs_${name}_eqtl_dosages.tsv.gz" --wgs-chrx "$G/wgs_${name}_eqtl_dosages.chrx.tsv" \
        --wgs-build "$build" --eqtls "$REF/gtex_cortex_eqtls.tsv" --clinical "$D/ROSMAP/ROSMAP_clinical.csv" \
        --biospecimen "$D/ROSMAP/ROSMAP_biospecimen_metadata.csv" \
        "$D/DiverseCohorts/AMP-AD_DiverseCohorts_biospecimen_metadata.csv" \
        --actions "$V2/qc/genotype_actions/restricted/genotype_actions.tsv" --output "$OUTDIR/$name"
done
echo "======================================== genotype check against AMP-AD 1.0 WGS"
rp pipeline qc-check genotype config/data_qc/rosmap_mit.yaml config/data_qc/rosmap_green.yaml \
    config/data_qc/diverse_cohorts.yaml \
    --genotypes "$OUTDIR/amp1/restricted/wgs_amp1_eqtl_dosages.by_individual.tsv.gz" \
    --eqtls "$REF/gtex_cortex_eqtls.tsv" --genotype-build 37 --sample-ids match --output "$OUTDIR/check_amp1"
echo "======================================== genotype check against DiverseCohorts WGS"
rp pipeline qc-check genotype config/data_qc/diverse_cohorts.yaml \
    --genotypes "$OUTDIR/divco/restricted/wgs_divco_eqtl_dosages.by_individual.tsv.gz" \
    --eqtls "$REF/gtex_cortex_eqtls.tsv" --genotype-build 38 --sample-ids match --output "$OUTDIR/check_divco"

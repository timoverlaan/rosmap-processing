#!/bin/bash
# Submit the whole data QC from the repo root (login node; only calls sbatch):
#   bash slurm/submit_qc.sh
# Gene lists -> one scan per dataset -> one job per check (after the scans),
# the shared-donor and sex-effect comparisons, the metadata check (obs only)
# and the SEA-AD hash, each its own job. Outputs: data/v2/qc/<dataset>/ and
# data/v2/qc/{identity,sex-de-compare}/.
set -euo pipefail
klist -s || { echo "no valid Kerberos ticket: kinit first"; exit 1; }

ALL="config/data_qc/rosmap_mit.yaml config/data_qc/rosmap_green.yaml config/data_qc/seaad.yaml config/data_qc/diverse_cohorts.yaml config/data_qc/psychad.yaml"
SHARED="config/data_qc/rosmap_mit.yaml config/data_qc/rosmap_green.yaml config/data_qc/diverse_cohorts.yaml config/data_qc/psychad.yaml"
submit() { sbatch --parsable --exclude=gpu23,gpu38 "$@"; }

genes=$(submit slurm/qc_genes.sh)
echo "qc_genes: $genes"
scans=()
for spec in $ALL; do
    mem=64GB
    case "$spec" in *psychad*) mem=128GB ;; esac
    id=$(submit --dependency=afterok:"$genes" --mem="$mem" --job-name="qc_scan_$(basename "$spec" .yaml)" \
        slurm/qc_scan.sh "$spec")
    echo "qc_scan $(basename "$spec" .yaml): $id"
    scans+=("$id")
done
after="afterok:$(IFS=:; echo "${scans[*]}")"
for check in totals cells sex markers ad-de; do
    # shellcheck disable=SC2086
    id=$(submit --dependency="$after" --job-name="qc_$check" slurm/qc_check.sh "$check" $ALL)
    echo "qc_check $check: $id"
done
# shellcheck disable=SC2086
sexde=$(submit --dependency="$after" --job-name=qc_sex-de slurm/qc_check.sh sex-de $ALL)
echo "qc_check sex-de: $sexde"
# shellcheck disable=SC2086
id=$(submit --dependency=afterok:"$sexde" --job-name=qc_sex-de-compare slurm/qc_check.sh sex-de-compare $ALL)
echo "qc_check sex-de-compare: $id"
# shellcheck disable=SC2086
id=$(submit --dependency="$after" --mem=96GB --job-name=qc_identity slurm/qc_check.sh identity $SHARED)
echo "qc_check identity: $id"
id=$(submit --job-name=qc_metadata slurm/qc_check.sh metadata \
    config/data_qc/diverse_cohorts.yaml config/data_qc/psychad.yaml)
echo "qc_check metadata: $id"
id=$(submit slurm/seaad_hash.sh)
echo "seaad_hash: $id"

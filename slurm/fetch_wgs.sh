#!/bin/bash
# Genotypes at the GTEx cortex eQTL sites (and chrX heterozygosity) from WGS
# joint calls on Synapse, through byte ranges of the indexed VCFs
# (scripts/fetch_wgs_sites.py): a few GB transferred, not the files. Needs
# token.txt in the repo root. Submit from the repo root:
#   sbatch slurm/fetch_wgs.sh amp1    # AMP-AD 1.0 WGS (syn11707420; ROSMAP, MSBB, Mayo; GRCh37)
#   sbatch slurm/fetch_wgs.sh divco   # DiverseCohorts WGS + AMP-AD 1.0 overlap (syn66709613; GRCh38)
# Extra arguments go to the script (e.g. --chromosomes 22 --max-gb 0.5).
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --output=slurm/out/%j_fetch_wgs.out
#SBATCH --error=slurm/out/%j_fetch_wgs.out

EXTRA_BINDS="--bind token.txt:/opt/app/token.txt --bind scripts/:/opt/app/scripts/"
source slurm/common.sh

name="${1:-}"
case "$name" in
    amp1) folder=syn11707420; build=37 ;;
    divco) folder=syn66709613; build=38 ;;
    *) echo "usage: sbatch slurm/fetch_wgs.sh amp1|divco [script arguments]"; exit 1 ;;
esac
shift
in_container pixi run python scripts/fetch_wgs_sites.py "$folder" "$REF/gtex_cortex_eqtls.tsv" \
    "$V2/genotypes/restricted/wgs_${name}_eqtl_dosages.tsv.gz" --build "$build" "$@"

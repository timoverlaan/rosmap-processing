#!/bin/bash
# Size, MD5 and S3 multipart ETag of the SEA-AD download, against the S3
# object version it came from. The 2024-02-13 file has had a delete marker in
# the public bucket since 2026-06-23 (replaced by a 2026-06-22 release); its
# old version is listed with ETag aee9da58948f244edb388e0abe1e96f5-7237 and
# 37938782321 bytes (ListObjectVersions, checked 2026-11-16).
#   sbatch slurm/seaad_hash.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --output=slurm/out/%j_seaad_hash.out
#SBATCH --error=slurm/out/%j_seaad_hash.out

EXTRA_BINDS="--bind scripts/:/opt/app/scripts/"
source slurm/common.sh

in_container pixi run python scripts/s3_etag.py "$SEAAD_ORIGINAL" \
    --expected-etag aee9da58948f244edb388e0abe1e96f5-7237 --expected-size 37938782321

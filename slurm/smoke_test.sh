#!/bin/bash
# Smoke test in the cluster container, to run before the real jobs after a code
# change: the unit tests that use synthetic data only (never the fixtures that
# load real files), including the seurat_v3 HVG path that needs the Linux
# scikit-misc build, and a check that R can load SeuratDisk.
# Submit from the repo root: sbatch slurm/smoke_test.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --output=slurm/out/%j_smoke_test.out
#SBATCH --error=slurm/out/%j_smoke_test.out

EXTRA_BINDS="--bind tests/:/opt/app/tests/ --bind pyproject.toml:/opt/app/pyproject.toml"
source slurm/common.sh

in_container pixi run Rscript -e 'suppressMessages(library(SeuratDisk)); cat("SeuratDisk", as.character(packageVersion("SeuratDisk")), "OK\n")'

in_container pixi run python -m pytest -q -p no:cacheprovider -rs \
    tests/test_utils_counts.py \
    tests/test_data_prepare.py \
    tests/test_data_h5seurat.py \
    tests/test_data_genes.py \
    tests/test_data_download.py \
    tests/test_core_scanpy_pipeline.py \
    tests/test_pipelines_hvg_compare.py \
    tests/test_pipelines_mathys2023_for_adprs.py \
    tests/test_utils_constants.py \
    tests/test_utils_logging.py

#!/bin/bash
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=insy,general # Request partition. Default is 'general'
#SBATCH --qos=long         # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=12:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1          # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=4   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem=500GB          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --mail-type=END     # Set mail type to 'END' to receive a mail when the job finishes.
#SBATCH --output=slurm/out/%j_lieke_allgenes_hvg_postnorm.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm/out/%j_lieke_allgenes_hvg_postnorm.out # Set name of error log. %j is the Slurm jobId
/usr/bin/scontrol show job -d "$SLURM_JOB_ID"  # check sbatch directives are working

# Scenario: intersection of Lieke's allgenes x ROSMAP_MIT x SeaAD (~27936 genes),
# then select 2000 HVGs AFTER normalization using the seurat (v1) flavor.
# Compare with lieke_allgenes_hvg.sh which does the same but pre-normalization (seurat_v3).
# Step 1: ROSMAP_MIT — intersect with allgenes.txt, normalize, then select HVGs within that subset
# Step 2: SeaAD    — uses the ROSMAP_MIT output as gene reference (exact same HVGs)

apptainer exec --writable-tmpfs --pwd /opt/app --containall \
	--bind src/:/opt/app/src/ \
	--bind data/:/opt/app/data/ \
	--env PYTHONPATH=/opt/app/src \
	./container_pixi_0-1-3.sif pixi run python -m rosmap_processing pipeline scanpy \
		data/raw/ROSMAP_MIT/combined.h5ad \
		data/processed/rosmap_mit_lieke_allgenes_hvg2k_postnorm_k30.h5ad \
		--import-genes data/lieke_allgenes_intersection.txt \
		--hvg-after-normalize 2000 \
		--k-neighbors 30

apptainer exec --writable-tmpfs --pwd /opt/app --containall \
	--bind src/:/opt/app/src/ \
	--bind data/:/opt/app/data/ \
	--env PYTHONPATH=/opt/app/src \
	./container_pixi_0-1-3.sif pixi run python -m rosmap_processing pipeline scanpy \
		data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad \
		data/seaAD/PFC/RNAseq/seaad_lieke_allgenes_hvg2k_postnorm_k30.h5ad \
		--layer UMIs \
		--import-genes data/processed/rosmap_mit_lieke_allgenes_hvg2k_postnorm_k30.h5ad \
		--k-neighbors 30

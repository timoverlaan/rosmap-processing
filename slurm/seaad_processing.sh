#!/bin/sh
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=insy,general # Request partition. Default is 'general' 
#SBATCH --qos=short         # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=1:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1          # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=16   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem=500GB          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --mail-type=END     # Set mail type to 'END' to receive a mail when the job finishes. 
#SBATCH --output=slurm/out/%j_all.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm/out/%j_all.out # Set name of error log. %j is the Slurm jobId
/usr/bin/scontrol show job -d "$SLURM_JOB_ID"  # check sbatch directives are working

apptainer exec --writable-tmpfs --pwd /opt/app --containall \
	--bind src/:/opt/app/src/ \
	--bind data/:/opt/app/data/ \
	--bind out/:/opt/app/out/ \
	./container_pixi_0-1-2.sif pixi run python -u src/scanpy_pipeline.py data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad data/seaAD/PFC/RNAseq/seaad1000_k30.h5ad --layer UMIs --n_genes 1000 k-neighbors 30

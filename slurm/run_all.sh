#!/bin/sh
#SBATCH --account=ewi-insy-prb-dbl
#SBATCH --partition=insy,general # Request partition. Default is 'general' 
#SBATCH --qos=medium         # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=36:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1          # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=8   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem=300GB          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --mail-type=END     # Set mail type to 'END' to receive a mail when the job finishes. 
#SBATCH --output=slurm/out/%j_download.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm/out/%j_download.out # Set name of error log. %j is the Slurm jobId
#SBATCH --gres=gpu
/usr/bin/scontrol show job -d "$SLURM_JOB_ID"  # check sbatch directives are working

module use /opt/insy/modulefiles
module load cuda/12.4

cd /tudelft.net/staff-umbrella/scGraphNN/rosmap-processing/

apptainer exec --nv --writable-tmpfs --pwd /opt/app --bind /tudelft.net/:/tudelft.net/ ./container_pixi.sif \
	sh ./src/run_all.sh
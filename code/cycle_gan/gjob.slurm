#!/bin/bash
#SBATCH --mem=32g
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1    # <- match to OMP_NUM_THREADS
#SBATCH --partition=gpuA100x4      # <- or one of: gpuA100x4 gpuA40x4 gpuA100x8 gpuMI100x8
#SBATCH --account=bbhg-delta-gpu
#SBATCH --time=06:00:00      # hh:mm:ss for the job

 
module purge
module load gcc anaconda3 cuda cudnn
module list
echo "job is starting on `hostname`"
srun python3 generator.py $idfn $model_type $exp_dir_fn $gaussian


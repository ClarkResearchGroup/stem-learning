#!/bin/bash
#SBATCH --mem=64g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1    # <- match to OMP_NUM_THREADS
#SBATCH --partition=cpu      # <- or one of: gpuA100x4 gpuA40x4 gpuA100x8 gpuMI100x8
#SBATCH --account=bbhg-delta-cpu
#SBATCH --time=02:00:00      # hh:mm:ss for the job
#SBATCH --output=slurm-%x.%j.out

 
module purge
module load gcc anaconda3
module list  
echo "job is starting on `hostname`"
srun python3 1_preprocess_express.py $idfn $model_type


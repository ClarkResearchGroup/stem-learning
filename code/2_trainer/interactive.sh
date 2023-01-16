srun \
 --mem=32g \
 --nodes=1 \
 --ntasks-per-node=1 \
 --cpus-per-task=1 \
 --partition=gpuA100x4 \
 --account=bbhg-delta-gpu \
 --gpus-per-node=1 \
 --gpus-per-task=1 \
 --gpu-bind=verbose,per_task:1 \
 --pty \
 singularity run --nv \
/sw/external/NGC/tensorflow:22.02-tf2-py3 /bin/bash

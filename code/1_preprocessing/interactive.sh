srun --account=bbhg-delta-cpu --partition=cpu-interactive \
  --nodes=1 --tasks=1 --tasks-per-node=1 \
  --cpus-per-task=1 --mem=16g \
  --pty bash

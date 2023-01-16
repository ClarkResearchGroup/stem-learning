#!/bin/bash
for mat in Gr STO
do
    for num in 1 2 3
    do
        echo $mat $num
        sbatch --job-name=${mat}_${num} --export=ALL,mat=${mat},num=${num} stogr_job.slurm
    done
done

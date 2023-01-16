#!/bin/bash

for idfn in 107_1840 107_1847 211_1528 211_1922 211_2124 107_1750 107 212 all
do
    for model_type in exp sim
    do
        echo $idfn $model_type
        sbatch --job-name=${idfn}_${model_type} --export=ALL,idfn=$idfn,model_type=$model_type job.slurm
    done
done

#idfn=107_1840
#model_type=sim
#sbatch --job-name=${idfn}_${model_type} --export=ALL,idfn=$idfn,model_type=$model_type job.slurm

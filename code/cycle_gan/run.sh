#!/bin/bash

#for idfn in 107  212 all
#do
#
#    exp_fn="${idfn}_normalized"
#    echo $exp_fn
#    sbatch --job-name=$idfn --export=ALL,exp_fn=$exp_fn,idfn=$idfn job.slurm
#done
#
#for idfn in 107_1840 107_1847 211_1528 211_1922 211_2124 107_1750 #
#do
#
#    exp_fn="mini_experiment_batch/RR_$idfn"
#    echo $exp_fn
#    sbatch --job-name=$idfn --export=ALL,exp_fn=$exp_fn,idfn=$idfn job.slurm
#done

for exp_fn in STO Gr
do
    for idfn in 0 1 2 3
    do
        for gaussian in 0.1 0.2 0.4 0.8
        do
            echo $exp_fn $idfn
            sbatch --job-name=${exp_fn}_${idfn}_${gaussian} --export=ALL,exp_fn=${exp_fn},idfn=${idfn},gaussian=${gaussian} job.slurm
        done
    done
done

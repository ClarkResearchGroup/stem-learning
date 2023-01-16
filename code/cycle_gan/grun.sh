#!/bin/bash

#for idfn in 107_1840 107_1847 211_1528 211_1922 211_2124 107_1750 107 212 all
#do
#    for model_type in exp sim
#    do
#        echo $idfn $model_type
#        sbatch --job-name=$exp_fn --export=ALL,idfn=$idfn,model_type=$model_type gjob.slurm
#    done
#done

for exp_dir_fn in STO Gr
do
    for idfn in 0 1 2 3
    do
        for gaussian in 0.1 0.2 0.4 0.8
        do
            for model_type in exp sim
            do
                echo ${exp_dir_fn}_${idfn}_${gaussian}_${model_type}
                sbatch --job-name=${exp_dir_fn}_${idfn}_${gaussian}_${model_type} --export=ALL,idfn=$idfn,model_type=$model_type,exp_dir_fn=$exp_dir_fn,gaussian=$gaussian gjob.slurm
            done
        done
    done
done

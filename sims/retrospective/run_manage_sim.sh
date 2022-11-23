#!/bin/bash

# num_combos is number of parameter combinations (see manage_sim.R)
# 2 is the number of reps per combos
# 3 is number of reps performed per job
# 1 is sim name
num_combos=48
njobs=`expr $2 / $3 \* $num_combos`

qsub -cwd -l h="biostat-b34|biostat-b35|biostat-b36|biostat-b37" -l h_vmem=16G -e iotrash/ -o iotrash/ -t 1-$njobs /home/users/cwolock/stack_supplementary/sims/retrospective/call_manage_sim.sh $1 $2 $3

#!/bin/bash

# num_combos is number of parameter combinations (see manage_sim.R)
#    this does not need to be changed
num_combos=16
njobs=`expr 100 \* $num_combos`

sbatch --array=1-$njobs -p short -t 6:00:00 -e ./iotrash/s-%A_%a.out -o ./iotrash/s-%A_%a.out /home/cwolock/stack_supplementary/sims/form_comparison/call_manage_sim.sh

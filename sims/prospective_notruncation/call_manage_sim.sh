#!/bin/bash

Rscript /home/users/cwolock/stack_supplementary/sims/prospective_notruncation/manage_sim.R --sim-name $1 --nreps-total $2 --nreps-per-job $3 --scheduler $4

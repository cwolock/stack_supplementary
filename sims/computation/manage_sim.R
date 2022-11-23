#!/usr/local/bin/Rscript
suppressMessages(library(survML))
suppressMessages(library(survSuperLearner))
suppressMessages(library(SuperLearner))
suppressMessages(library(survival))
suppressMessages(library(argparse))
suppressMessages(library(dplyr))

source("/home/ec2-user/stack_supplementary/sims/computation/do_one_compare.R")
source("/home/ec2-user/stack_supplementary/sims/generate_data.R")

sim_name <- "timing_sim"

nreps_total <- 100

n_trains <- c(500)
dgps <- c("leftskew")
estimators <- c("stackG_fine", "stackG_medium", "stackG_coarse",
                "stackL_fine", "stackL_medium", "stackL_coarse",
                "coxph", "survSL")

param_grid <- expand.grid(rep_id = 1:nreps_total,
			  estimator = estimators,
                          dgp = dgps,
                          n_train = n_trains)
param_grid$mc_id <- 1:nrow(param_grid)

f <- function(i){
  return(do_one(sim_name = sim_name,
		mc_id = param_grid$mc_id[i],
		estimator = param_grid$estimator[i],
                dgp = param_grid$dgp[i],
                n_train = param_grid$n_train[i]))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

res <- lapply(1:nrow(param_grid), f)

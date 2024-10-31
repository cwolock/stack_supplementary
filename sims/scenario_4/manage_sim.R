#!/usr/local/bin/Rscript
.libPaths(c(
  "/home/cwolock/R_lib",
  .libPaths()
))
suppressMessages(library(survML))
suppressMessages(library(SuperLearner))
suppressMessages(library(survival))
suppressMessages(library(dplyr))
suppressMessages(library(LTRCforests))

sim_name <- "scenario_4"
nreps_total <- 100
nreps_per_job <- 1

source("/home/cwolock/stack_supplementary/sims/scenario_4/do_one.R")
source("/home/cwolock/stack_supplementary/sims/generate_data.R")

n_trains <- c(250, 500, 750, 1000)
dgps <- c("leftskew", "rightskew")
estimators <- c("stackG_fine", "stackG_medium", "stackG_coarse",
                "stackL_fine", "stackL_medium", "stackL_coarse",
                "coxph", "LTRCforests")

njobs_per_combo <- nreps_total/nreps_per_job

param_grid <- expand.grid(mc_id = 1:njobs_per_combo,
                          dgp = dgps,
                          n_train = n_trains,
                          estimator = estimators)

job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

current_dynamic_args <- param_grid[job_id, ]

current_seed <- job_id
set.seed(current_seed)
output <- replicate(nreps_per_job,
                    do_one(n_train = current_dynamic_args$n_train,
                           estimator = current_dynamic_args$estimator,
                           dgp = current_dynamic_args$dgp),
                    simplify = FALSE)
sim_output <- lapply(as.list(1:length(output)),
                     function(x) tibble::add_column(output[[x]]))
sim_output_tib <- do.call(rbind.data.frame, sim_output)
file_name <- paste0("output/", sim_name, "_", job_id, ".rds")
saveRDS(sim_output_tib, file = file_name)

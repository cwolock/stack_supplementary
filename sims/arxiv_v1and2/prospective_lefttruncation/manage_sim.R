#!/usr/local/bin/Rscript
.libPaths(c(
  "/home/cwolock/R_lib",
  .libPaths()
))
suppressMessages(library(survML))
suppressMessages(library(survSuperLearner))
suppressMessages(library(SuperLearner))
suppressMessages(library(survival))
suppressMessages(library(argparse))
suppressMessages(library(dplyr))

parser <- ArgumentParser()
parser$add_argument("--sim-name",
                    default = "sim",
                    help = "Name of simulation")
parser$add_argument("--nreps-total",
                    type = "double",
                    default = 100,
                    help = "Number of replicates for each set of params")
parser$add_argument("--nreps-per-job",
                    type = "double",
                    default = 10,
                    help = "number of replicates per job")
parser$add_argument("--scheduler",
                    default = "sge",
                    help = "Job scheduler")
args <- parser$parse_args()

if (args$scheduler == "slurm"){
  source("/home/cwolock/stack_supplementary/sims/prospective_lefttruncation/do_one.R")
  source("/home/cwolock/stack_supplementary/sims/generate_data.R")
} else if (args$scheduler == "sge"){
  source("/home/users/cwolock/stack_supplementary/sims/prospective_lefttruncation/do_one.R")
  source("/home/users/cwolock/stack_supplementary/sims/generate_data.R")
}

n_trains <- c(250, 500, 750, 1000)
dgps <- c("leftskew", "rightskew")
estimators <- c("stackG_fine", "stackG_medium", "stackG_coarse",
                "stackL_fine", "stackL_medium", "stackL_coarse",
                "coxph")

njobs_per_combo <- args$nreps_total/args$nreps_per_job

param_grid <- expand.grid(mc_id = 1:njobs_per_combo,
                          dgp = dgps,
                          n_train = n_trains,
                          estimator = estimators)

if (args$scheduler == "slurm"){
  job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
} else if (args$scheduler == "sge"){
  job_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
}

current_dynamic_args <- param_grid[job_id, ]

current_seed <- job_id
set.seed(current_seed)
output <- replicate(args$nreps_per_job,
                    do_one(n_train = current_dynamic_args$n_train,
                           estimator = current_dynamic_args$estimator,
                           dgp = current_dynamic_args$dgp),
                    simplify = FALSE)
sim_output <- lapply(as.list(1:length(output)),
                     function(x) tibble::add_column(output[[x]]))
sim_output_tib <- do.call(rbind.data.frame, sim_output)
file_name <- paste0("output/", args$sim_name, "_", job_id, ".rds")
saveRDS(sim_output_tib, file = file_name)

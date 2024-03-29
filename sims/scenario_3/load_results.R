#!/usr/local/bin/Rscript

sim_name <- "scenario_3"
nreps_total <- 100
nreps_per_job <- 1

## set up directories for output, plots
output_dir <- "output/"

## set up parameter grid
n_trains <- c(250, 500, 750, 1000)
dgps <- c("leftskew", "rightskew")
estimators <- c("stackG_fine", "stackG_medium", "stackG_coarse",
                "stackL_fine", "stackL_medium", "stackL_coarse")

## number of monte-carlo iterations per job
nreps_per_combo <- nreps_total/nreps_per_job
## set up grid of parameters
param_grid <- expand.grid(mc_id = 1:nreps_per_combo,
                          n_train = n_trains,
                          estimator = estimators,
                          dgp = dgps)

## names of files to read in
output_nms <- paste0(sim_name, "_", 1:dim(param_grid)[1], ".rds")
avail_nms <- list.files(output_dir, pattern = paste0(sim_name, "_*"))
names_to_try <- output_nms[which(output_nms %in% avail_nms)]
print(output_nms[which(!(output_nms %in% avail_nms))])
## list of output
output_lst <- lapply(paste0(output_dir, names_to_try), readRDS)
## make it a matrix
output_df <- do.call(rbind.data.frame, output_lst)
saveRDS(output_df, paste0(sim_name, ".rds"))

#!/usr/local/bin/Rscript

suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "prospective_truncation",
                    help = "name of simulation")
parser$add_argument("--nreps-total", type = "double", default = 100,
                    help = "number of replicates for each set of params")
parser$add_argument("--nreps-per-job", type = "double", default = 1,
                    help = "number of replicates per job")
args <- parser$parse_args()

## set up directories for output, plots
output_dir <- "output/"

## set up parameter grid
n_trains <- c(250, 500, 750, 1000)
dgps <- c("leftskew", "rightskew")
estimators <- c(#"stackG_fine", "stackG_medium", "stackG_coarse",
                "stackL_fine", "stackL_medium", "stackL_coarse")#,
                #"coxph")
cens <- c(0.25)
## number of monte-carlo iterations per job
nreps_per_combo <- args$nreps_total/args$nreps_per_job
## set up grid of parameters
param_grid <- expand.grid(mc_id = 1:nreps_per_combo,
                          n_train = n_trains,
                          estimator = estimators,
                          dgp = dgps,
			  cens = cens)

## names of files to read in
output_nms <- paste0(args$sim_name, "_", 1:dim(param_grid)[1], ".rds")
avail_nms <- list.files(output_dir, pattern = paste0(args$sim_name, "_*"))
names_to_try <- output_nms[which(output_nms %in% avail_nms)]
print(length(output_nms) - length(avail_nms))
print(output_nms[which(!(output_nms %in% avail_nms))])
print(avail_nms[which(!(avail_nms %in% output_nms))])
## list of output
output_lst <- lapply(paste0(output_dir, names_to_try), readRDS)
## make it a matrix
output_df <- do.call(rbind.data.frame, output_lst)

saveRDS(output_df, paste0(args$sim_name, ".rds"))

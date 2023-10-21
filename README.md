# Supplementary materials for global survival stacking paper

This repository contains code to reproduce the analyses in ["A framework for leveraging machine learning tools to estimate personalized survival curves"](https://arxiv.org/abs/2211.03031) by Wolock, Gilbert, Simon, and Carone (2022+). All analyses were implemented using the `R` package `survML`, which is available on CRAN (stable release) and [here](https://github.com/cwolock/survML) (development version).

## Appendix

The appendix contains all technical details, as well as simulation results not included in the main text and details on the publicly available datasets. 

## Code

The code directory contains all code needed to replicate the simulations, analysis of publicly available data, and analysis of STEP data. 

The code depends on the following `R` packages: 

* `cowplot`: Available on CRAN.
* `ggpubr`: Available on CRAN.
* `LTRCforests`: Removed from CRAN as of March 2023. Available on Github at https://github.com/weichiyao/TimeVaryingData_LTRCforests/pkg/LTRCforests.
* `squash`: Available on CRAN.
* `SuperLearner`: Available on CRAN.
* `survival`: Available on CRAN.
* `survML`: Available on CRAN.
* `survSuperLearner`: Available on Github at https://github.com/tedwestling/survSuperLearner.
* `tidyverse`: Available on CRAN.

Note that `survSuperLearner` and `LTRCforests` are only used a comparator methods in the simulation studies. If you have difficulty installing them, you can simply leave those methods out. 

### Simulations

These simulations were performed on a Linux cluster using the Slurm job scheduler. For different job schedulers, the code is flagged where changes may be needed. Running these simulations locally, while possible, will be extremely time consuming. 

The `simulations` directory contains subdirectories corresponding to each set of numerical experiments from the paper. Each set of experiments uses the following workflow: 

1. Navigate to the desired subdirectory. 

2. Create an `output` directory to contain simulation results and an `iotrash` directory to contain error files by running `mkdir output` and `mkdir iotrash`.  

3. Run `bash run_manage_sim.sh`. This will run 100 Monte Carlo replicates for each parameter combination in the experiment. The total number of replicates for each experiment is given by 100 times the number of parameter combinations (see the `run_manage_sim.sh` script).  

4. Wait for all jobs to complete. 

5. Run `Rscript load_results.R` to compile the results in a single file called `{subdirectory name}.rds`. 

When all simulations have finished and results have been compiled, run `Rscript figures.R` to create the figures corresponding to this set of experiments. 

### Analysis of publicly available data

The results of these analysis are shown in Table 3 of the main text. 

* FLCHAIN: Data are available through the `survival` package. Run `code/public_data_analysis/flchain.R`. 

* GBSG: Data are available through the `survival` package. Run `code/public_data_analysis/gbsg.R`.

* METABRIC: Data are available through the `DeepSurv` package but are provided here for convenience as `code/public_data_analysis/metabric.rds`. Run `code/public_data_analysis/metabric.R`. 

* NWTCO: Data are available through the `survival` package. Run `code/public_data_analysis/nwtco.R`.

* SUPPORT: Data are available on the Vanderbilt Biostatistics website but are provided here for convenience as `code/public_data_analysis/support.rds`. Run `code/public_data_analysis/support.R`. 

### Analysis of STEP data

Run `code/step_analysis/step_analysis.R` to perform the analysis and produce Figures 3 and 4 from the main text. 

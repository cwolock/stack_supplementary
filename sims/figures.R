#!/usr/local/bin/Rscript

## nice plots
suppressMessages(library(ggplot2))
# library("cowplot")
# library("grid")
# library("gridExtra")
# library("ggpubr")
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))

# colors 
stackG_cols <- brewer.pal(n = 9, "Blues")[c(3,5,8)]
stackL_cols <- brewer.pal(n = 9, "Greens")[c(3,5,7)]
survSL_cols <- brewer.pal(n = 9, "YlOrRd")[5]
cox_cols <- brewer.pal(n = 12, "Paired")[11]

###############################
### prospective left truncation
###############################
wd <- "C:/Users/cwolo/Dropbox/UW/DISSERTATION/stack_supplementary/sims/prospective_lefttruncation"
dat <- readRDS(paste0(wd, "/sim_062222_ltrunc.rds"))
y_max <- max(log(dat$MSE_uni))
y_min <- min(log(dat$MSE_uni))
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stacking (all times grid)",
                                   stackG_medium = "Global stacking (0.025 grid)",
                                   stackG_coarse = "Global stacking (0.1 grid)",
                                   stackL_fine = "Local stacking (all times grid)",
                                   stackL_medium = "Local stacking (0.025 grid)",
                                   stackL_coarse = "Local stacking (0.1 grid)",
                                   coxph = "Cox"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))
pal <- c(gridSurv_cols, stackSurv_cols, cox_cols)
dat %>%
  ggplot(aes(x = factor(n_train), y = log(MSE_uni), fill = Estimator)) +
  geom_boxplot() +
  facet_wrap(~dgp) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ylab("log (MISE)") +
  xlab("Sample size") +
  ylim(y_min, y_max)+
  scale_fill_manual(values = pal)
ggsave('prospective_leftruncation.png', device='png', width=12,
       height=6, units='in', dpi=800)

TC_rates <- integrate_df %>% group_by(dgp) %>%
  summarize(1-mean(censoring_rate),
            mean(truncation_rate))

#############################
### prospective no truncation
#############################
wd <- "C:/Users/cwolo/Dropbox/UW/DISSERTATION/stack_supplementary/sims/prospective_notruncation"
dat <- readRDS(paste0(wd, "/testing.rds"))
y_max <- max(log(dat$MSE_uni))
y_min <- min(log(dat$MSE_uni))
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stacking (all times grid)",
                                   stackG_medium = "Global stacking (0.025 grid)",
                                   stackG_coarse = "Global stacking (0.1 grid)",
                                   stackL_fine = "Local stacking (all times grid)",
                                   stackL_medium = "Local stacking (0.025 grid)",
                                   stackL_coarse = "Local stacking (0.1 grid)",
                                   survSL = "survSuperLearner",
                                   coxph = "Cox"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

pal <- c(stackG_cols, stackL_cols, survSL_cols, cox_cols)
dat %>%
  ggplot(aes(x = factor(n_train), y = log(MSE_uni), fill = Estimator)) +
  geom_boxplot() +
  facet_wrap(~dgp) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ylab("log (MISE)") +
  xlab("Sample size") +
  ylim(y_min, y_max)+
  scale_fill_manual(values = pal)
ggsave('prospective_notruncation.png', device='png', width=12,
       height=6, units='in', dpi=800)

TC_rates <- dat %>% group_by(dgp) %>%
  summarize(1-mean(censoring_rate),
            mean(truncation_rate))

#################
### retrospective
#################
wd <- "C:/Users/cwolo/Dropbox/UW/DISSERTATION/stack_supplementary/sims/retrospective"
dat <- readRDS(paste0(wd, "/testing.rds"))
y_max <- max(log(integrate_df$MSE_uni))
y_min <- min(log(integrate_df$MSE_uni))
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stacking (all times grid)",
                                   stackG_medium = "Global stacking (0.025 grid)",
                                   stackG_coarse = "Global stacking (0.1 grid)",
                                   stackL_fine = "Local stacking (all times grid)",
                                   stackL_medium = "Local stacking (0.025 grid)",
                                   stackL_coarse = "Local stacking (0.1 grid)"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))
pal <- c(gridSurv_cols, stackSurv_cols)
integrate_df %>%
  ggplot(aes(x = factor(n_train), y = log(MSE_uni), fill = Estimator)) +
  geom_boxplot() +
  facet_wrap(~dgp) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ylab("log (MISE)") +
  xlab("Sample size") +
  ylim(y_min, y_max)+
  scale_fill_manual(values = pal)
ggsave('retrospective.png', device='png', width=12,
       height=6, units='in', dpi=800)

TC_rates <- integrate_df %>% group_by(dgp) %>%
  summarize(1-mean(censoring_rate),
            mean(truncation_rate))

###################
### form comparison
###################
dat <- readRDS(paste0(wd, "/sim_072522_ltrunc_formcompare.rds"))
dat <- readRDS(paste0(wd, "/sim_100322_xgb_formcompare.rds"))
integrate_df <- dat
# time integrated uniformly
y_max <- max(log(integrate_df$MSE_uni))
y_min <- min(log(integrate_df$MSE_uni))
integrate_df <- integrate_df %>% filter(estimator == "gridSurv_exp_strat" | 
                                          estimator == "gridSurv_PI_strat")
integrate_df <- integrate_df %>%
  mutate(estimator = recode_factor(estimator,
                                   gridSurv_exp_strat = "Exponential",
                                   gridSurv_PI_strat = "Product integral"))
integrate_df$Estimator <- integrate_df$estimator
integrate_df <- integrate_df %>%
  mutate(dgp = recode_factor(dgp,
                             early = "Right skew",
                             late = "Left skew"))
integrate_df %>%
  ggplot(aes(x = factor(n_train), y = log(MSE_uni), fill = Estimator)) +
  geom_boxplot() +
  facet_wrap(~dgp) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ylab("log (MISE)") +
  xlab("Sample size") +
  ylim(y_min, y_max) +
  scale_fill_manual(values = brewer.pal(n = 8, "Dark2")[c(2,3)])
ggsave('LTRC_formcompare.png', device='png', width=12,
       height=6, units='in', dpi=800)

times <- integrate_df %>% filter(n_train == 1000) %>%
  group_by(estimator) %>%
  summarize(mean(runtime),
            sd(runtime),
            median(runtime),
            IQR(runtime))

TC_rates <- integrate_df %>% group_by(dgp) %>%
  summarize(1-mean(censoring_rate),
            mean(truncation_rate))

incompat <- integrate_df %>% group_by(estimator, n_train) %>%
  summarize(mean(prop_incompat))


### AWS timing
dat1 <- readRDS("C:/Users/cwolo/Dropbox/UW/DISSERTATION/conditional_surv/code/AWS/beta_fast/sim_081122_notrunc_timing_xgb_1.rds")
dat2 <- readRDS("C:/Users/cwolo/Dropbox/UW/DISSERTATION/conditional_surv/code/AWS/beta_fast/sim_100322_notrunc_timing_xgb_2.rds")
dat<- rbind(dat1, dat2)

times <- dat %>%
  group_by(estimator) %>%
  summarize(mean(runtime),
            sd(runtime),
            median(runtime),
            IQR(runtime))
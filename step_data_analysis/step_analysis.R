# set up
library(survML)
library(survival)
library(tidyverse)
library(SuperLearner)
library(squash)
set.seed(123)

# load data (replace with appropriate path based on your environment)
blind <- read.csv("C:/Users/cwolo/Dropbox/UW/DISSERTATION/conditional_surv/step_analysis/master_mitt_males_17OCT2007.csv")
blind <- read.csv("step_data.csv")

# transform Ad5
# following previous analysis of treating <= 18 (undetectable) as 18
blind <- blind %>% mutate(l_BASEAD5 = log(BASEAD5))

# subset into treated and placebo, getting rid of 43 without circumcision info
blind <-  blind %>% filter(!is.na(F_CRCM_UPDATE))
blind_treat <- blind %>% filter(F_TREAT == 1 & !is.na(F_CRCM_UPDATE))
blind_plac <- blind %>% filter(F_TREAT == 0  & !is.na(F_CRCM_UPDATE))

##############
# descriptives
##############

# loss to follow-up
sum(blind_treat$F_EVINF == 1 & blind_treat$VACC1_EVINF <= 365)/nrow(blind_treat)
sum(blind_treat$F_EVINF == 1 & blind_treat$VACC1_EVINF <= 730)/nrow(blind_treat)
sum(blind_plac$F_EVINF == 1 & blind_plac$VACC1_EVINF <= 365)/nrow(blind_plac)
sum(blind_plac$F_EVINF == 1 & blind_plac$VACC1_EVINF <= 730)/nrow(blind_plac)

# benchmark covariate values
benchmark_ad5 <- 3:8
newX <- as.data.frame(expand.grid(F_CRCM_UPDATE = c(0,1), l_BASEAD5 = benchmark_ad5))

# follow-up times at which to estimate
benchmark_times_blind <- seq(0, 730, by = 5)

############################
### global stacking analysis
############################

# set up super learner library
tune <- list(ntrees = c(250, 500, 1000), max_depth = c(1,2),
             shrinkage = c(0.01), minobspernode = 1)
xgb_grid <- create.SL.xgboost(tune = tune)
SL.library <- c("SL.mean", "SL.gam", "SL.earth",
                "SL.ranger", "SL.glm.interaction", xgb_grid$names)

# set up B grid (to approximate product integral)
approx_grid1 <- sort(unique(blind_treat$VACC1_EVINF))
approx_grid2 <- sort(unique(blind_plac$VACC1_EVINF))
approx_grid1 <- approx_grid1[approx_grid1 != 0]
approx_grid2 <- approx_grid2[approx_grid2 != 0]

# compute fits
blind_treat_survML  <- survML::stackG(time = blind_treat$VACC1_EVINF,
                                      event = blind_treat$F_EVINF,
                                      X = blind_treat[,c("F_CRCM_UPDATE", "l_BASEAD5")],
                                      newX = newX,
                                      newtimes = benchmark_times_blind,
                                      time_grid_approx = approx_grid1,
                                      bin_size = 0.025,
                                      time_basis = "continuous",
                                      surv_form = "exp",
                                      SL_control = list(SL.library = SL.library,
                                                        V = 5))$S_T_preds
blind_plac_survML  <- survML::stackG(time = blind_plac$VACC1_EVINF,
                                     event = blind_plac$F_EVINF,
                                     X = blind_plac[,c("F_CRCM_UPDATE", "l_BASEAD5")],
                                     newX = newX,
                                     newtimes = benchmark_times_blind,
                                     time_grid_approx = approx_grid2,
                                     bin_size = 0.025,
                                     time_basis = "continuous",
                                     surv_form = "exp",
                                     SL_control = list(SL.library = SL.library,
                                                       V = 5))$S_T_preds
################
### cox analysis
################
blind_treat_cox_fit  <- survival::coxph(survival::Surv(VACC1_EVINF, F_EVINF) ~ F_CRCM_UPDATE*l_BASEAD5,
                                        data = blind_treat)
blind_treat_cox <- t(summary(survival::survfit(blind_treat_cox_fit,
                                               newdata = newX,
                                               se.fit = FALSE,
                                               conf.int = FALSE),
                             times=benchmark_times_blind)$surv)
blind_plac_cox_fit <- survival::coxph(survival::Surv(VACC1_EVINF, F_EVINF) ~ F_CRCM_UPDATE*l_BASEAD5,
                                      data = blind_plac)
blind_plac_cox <- t(summary(survival::survfit(blind_plac_cox_fit,
                                              newdata = newX,
                                              se.fit = FALSE,
                                              conf.int = FALSE),
                            times=benchmark_times_blind)$surv)

####################################
### compile results and make figures
####################################
fitted <- list(blind_treat_survML = blind_treat_survML,
               blind_plac_survML = blind_plac_survML,
               blind_treat_cox = blind_treat_cox,
               blind_plac_cox = blind_plac_cox,
               times = benchmark_times_blind,
               newX = newX)


# compute differences in survival probabilities at landmark times of 1 and 2 years
year_indices <- which(fitted$times %in% c(365, 730))

surv_diff_cox <- c((1-fitted$blind_treat_cox[,year_indices]) - (1-fitted$blind_plac_cox[,year_indices]))
surv_diff_survML <- c((1-fitted$blind_treat_survML[,year_indices]) - (1-fitted$blind_plac_survML[,year_indices]))
end <- data.frame(surv_diff = c(surv_diff_cox,
                                surv_diff_survML),
                  time = rep(c(rep("One year", nrow(fitted$newX)), rep("Two years", nrow(fitted$newX))), 2),
                  rbind(fitted$newX, fitted$newX),
                  method = c(rep("cox", 2*nrow(fitted$newX)), rep("survML", 2*nrow(fitted$newX))))

end$F_CRCM_UPDATE = factor(end$F_CRCM_UPDATE, levels = c(0,1), labels = c("No", "Yes"))
colnames(end) <- c("Excess infection risk in vaccine arm",
                   "Time",
                   "Circumcision",
                   "Baseline log(Ad5) titer",
                   "Method")

end <- end %>% mutate(Method = ifelse(Method == "cox", "Cox", "Global survival stacking"))

# make Figure 4
end %>% ggplot(aes(x = `Baseline log(Ad5) titer`, y = `Excess infection risk in vaccine arm`)) +
  geom_line(aes(color = `Circumcision`, linetype = `Method`), alpha = 0.7, size = 1) +
  geom_abline(slope = 0, intercept = 0, color = "black") +
  facet_wrap(~ `Time`) +
  scale_color_manual(values = greyscale(2, start = 0.4, end = 0)) +
  theme_bw() +
  theme(text = element_text(size = 16),
        strip.background = element_blank(),) +
  ylab("Risk difference (vaccine - placebo)")
ggsave("step_analysis.png", # save wherever you like
       device='png', width=12,
       height=4, units='in', dpi=300)


# set up estimated representative survival curve
example_indices <- c(1,2,5,6,9,10)
example_surv_plac <- data.frame(fitted$newX[example_indices,], fitted$blind_plac_survML[example_indices,])
example_surv_treat <- data.frame(fitted$newX[example_indices,], fitted$blind_treat_survML[example_indices,])
colnames(example_surv_plac)[3:ncol(example_surv_plac)] <- as.character(fitted$times)
colnames(example_surv_treat)[3:ncol(example_surv_treat)] <- as.character(fitted$times)

example_surv_plac <- example_surv_plac %>% pivot_longer(cols = -c(1,2),
                                                        names_to = "time",
                                                        values_to = "surv")

example_surv_treat <- example_surv_treat %>% pivot_longer(cols = -c(1,2),
                                                          names_to = "time",
                                                          values_to = "surv")

example_surv_plac$treatment <- "placebo"
example_surv_treat$treatment <- "treatment"

example_surv <- rbind(example_surv_plac, example_surv_treat)

example_surv$`Baseline log(Ad5) titer` <- as.factor(example_surv$l_BASEAD5)
example_surv$Circumcision <- factor(example_surv$F_CRCM_UPDATE,
                                    labels = c("No", "Yes"))
example_surv$time <- as.numeric(example_surv$time)
example_surv$treatment <- factor(example_surv$treatment, labels = c("Placebo", "Vaccine"))

# make Figure 3
example_surv %>% ggplot(aes(x = time, y = surv)) +
  geom_line(aes(color = Circumcision, linetype = `Baseline log(Ad5) titer`), alpha = 0.7, size = 1) +
  geom_abline(slope = 0, intercept = 0, color = "black") +
  scale_linetype_manual(values = c(1,2,3)) +
  facet_wrap(~ treatment) +
  scale_color_manual(values = greyscale(2, start = 0.4, end = 0)) +
  theme_bw() +
  theme(text = element_text(size = 16),
        strip.background = element_blank(),) +
  ylab("1 - estimated probability of HIV-1 diagnosis") +
  xlab("Time (days)")
ggsave("step_curves.png",
       device='png', width=12,
       height=6, units='in', dpi=300)

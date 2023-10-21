library(tidyverse)
library(squash)
library(ggpubr)
library(cowplot)

small_estimators <- c("Global stacking (ours)",
                      "Local stacking",
                      "survSL",
                      "Linear Cox",
                      "Gen. additive Cox",
                      "LTRC forests")

wd <- "C:/Users/cwolo/Dropbox/UW/DISSERTATION/stack_supplementary/sims/scratch/"
setwd(wd)

################################
### Partition size in terms of n
################################
dat <- readRDS("rates.rds")
dat <- dat %>% mutate(C_partition_size = factor(rate,
                                                levels = c("1/4", "1/3", "1/2", "2/3", "3/4", "1"),
                                                labels = c(expression(n^{1/4}),
                                                           expression(n^{1/3}),
                                                           expression(n^{1/2}),
                                                           expression(n^{2/3}),
                                                           expression(n^{3/4}),
                                                           expression(n))),
                      B_partition_size = factor(approx_rate,
                                                levels = c("0.25", "0.333", "0.5", "0.667", "0.75", "1"),
                                                labels = c(expression(n^{1/4}),
                                                           expression(n^{1/3}),
                                                           expression(n^{1/2}),
                                                           expression(n^{2/3}),
                                                           expression(n^{3/4}),
                                                           expression(n))),
                      dgp = recode_factor(dgp,
                                          rightskew = "Right skew",
                                          leftskew = "Left skew")) %>%
  mutate(B_partition_size = fct_rev(B_partition_size))

summ <- dat %>% group_by(C_partition_size, B_partition_size, dgp, n_train) %>%
  summarize(MISE = log(mean(MSE_uni))) %>%
  mutate(n_train = factor(n_train,
                          levels = c(250, 500, 750, 1000),
                          labels = c("n = 250", "n = 500", "n = 750", "n = 1000"))) %>%
  group_by(dgp) %>%
  mutate(MISE = (MISE - min(MISE))/(max(MISE) - min(MISE)))

# pal <- greyscale(length(unique(dat$C_partition_size)), start = 1, end = 0.4)
# dat %>% ggplot(aes(x = factor(n_train), y = log(MSE_uni))) +
#   geom_boxplot(aes(fill = factor(C_partition_size)),
#                outlier.shape = NA) +
#   # linetype = C_partition_size,
#   # group = interaction(B_partition_size, C_partition_size))) +
#   facet_grid(B_partition_size ~ dgp) +
#   scale_fill_manual(values = pal) +
#   theme_bw()

full_plot <- summ %>% ggplot(aes(x = C_partition_size, y = B_partition_size, fill = MISE)) +
  geom_tile() +
  facet_grid(factor(n_train) ~ dgp) +
  theme_bw() +
  scale_fill_distiller(type = "seq",
                       direction = 1,
                       palette = "Greys") +
  ylab("Approximation grid size") +
  xlab("Regression grid size") +
  scale_x_discrete(labels = ~parse(text = .x)) +
  scale_y_discrete(labels = ~parse(text = .x)) +
  theme(text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        axis.text.x=element_text(size=8, hjust=0.5,vjust=0.2))



# dat <- dat %>% filter(rate != "1/4")
# pal <- greyscale(length(unique(dat$n_train)), start = 1, end = 0.4)
#
# p1 <- dat %>% mutate(type = "integrated") %>%
#   ggplot(aes(x = `Partition size`,
#              y = log(MSE_uni),
#              position = factor(n_train),
#              fill = factor(n_train))) +
#   geom_boxplot(outlier.shape = NA) +
#   facet_grid(type~dgp) +
#   theme_bw() +
#   ylab("log (MISE)") +
#   xlab("Number of cutpoints in terms of sample size") +
#   scale_fill_manual(values = pal, name = "Sample size") +
#   scale_x_discrete(labels = ~parse(text = .x)) +
#   theme(text = element_text(size = 10),
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         legend.position = "none",
#         axis.text.x=element_text(size=8, hjust=0.5,vjust=0.2))
#
# p2 <- dat %>% mutate(type = "50th percentile") %>%
#   ggplot(aes(x =  `Partition size`,
#              y = log(landmark_MSE_50),
#              position = factor(n_train),
#              fill = factor(n_train))) +
#   geom_boxplot(outlier.shape = NA) +
#   facet_grid(type~dgp) +
#   theme_bw() +
#   ylab("log (MSE)") +
#   xlab("Number of cutpoints in terms of sample size") +
#   scale_fill_manual(values = pal, name = "Sample size") +
#   scale_x_discrete(labels = ~parse(text = .x)) +
#   theme(text = element_text(size = 10),
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         legend.position = "none",
#         axis.text.x=element_text(size=8, hjust=0.5,vjust=0.2))
#
# p3 <- dat %>% mutate(type = "75th percentile") %>%
#   ggplot(aes(x =  `Partition size`,
#              y = log(landmark_MSE_75),
#              position = factor(n_train),
#              fill = factor(n_train))) +
#   geom_boxplot(outlier.shape = NA) +
#   facet_grid(type~dgp) +
#   theme_bw() +
#   ylab("log (MSE)") +
#   xlab("Number of cutpoints in terms of sample size") +
#   scale_fill_manual(values = pal, name = "Sample size") +
#   scale_x_discrete(labels = ~parse(text = .x)) +
#   theme(text = element_text(size = 10),
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         legend.position = "none",
#         axis.text.x=element_text(size=8, hjust=0.5,vjust=0.2))
#
# p4 <- dat%>% mutate(type = "90th pecentile") %>%
#   ggplot(aes(x =  `Partition size`,
#              y = log(landmark_MSE_90),
#              position = factor(n_train),
#              fill = factor(n_train))) +
#   geom_boxplot(outlier.shape = NA) +
#   facet_grid(type~dgp) +
#   theme_bw() +
#   ylab("log (MSE)") +
#   xlab("Number of cutpoints in terms of sample size") +
#   scale_fill_manual(values = pal, name = "Sample size") +
#   scale_x_discrete(labels = ~parse(text = .x)) +
#   theme(text = element_text(size = 10),
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         legend.position = "none",
#         axis.text.x=element_text(size=8, hjust=0.5,vjust=0.2))
#
# y_max <- max(layer_scales(p1)$y$range$range[2],
#              layer_scales(p2)$y$range$range[2],
#              layer_scales(p3)$y$range$range[2],
#              layer_scales(p4)$y$range$range[2])
# y_min <- min(layer_scales(p1)$y$range$range[1],
#              layer_scales(p2)$y$range$range[1],
#              layer_scales(p3)$y$range$range[1],
#              layer_scales(p4)$y$range$range[1])
#
# four_panel_plot <- ggarrange(
#   p1 + ylim(y_min, y_max) +theme(legend.position = "none",
#                                  axis.title.y = element_text(size = 14),
#                                  axis.title.x = element_blank(),
#                                  axis.text.x = element_blank(),
#                                  axis.ticks.x = element_blank(),
#                                  strip.text = element_text(size = 12),
#                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
#   p2 + ylim(y_min, y_max) +theme(legend.position = "none",
#                                  axis.title.y = element_text(size = 14),
#                                  strip.text.x = element_blank(),
#                                  strip.text.y = element_text(size = 12),
#                                  axis.title.x = element_blank(),
#                                  axis.text.x = element_blank(),
#                                  axis.ticks.x = element_blank(),
#                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
#   p3 +ylim(y_min, y_max) + theme(legend.position = "none",
#                                  axis.title.y = element_text(size = 14),
#                                  strip.text.x = element_blank(),
#                                  strip.text.y = element_text(size = 12),
#                                  axis.title.x = element_blank(),
#                                  axis.text.x = element_blank(),
#                                  axis.ticks.x = element_blank(),
#                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
#   p4 + ylim(y_min, y_max) +theme(legend.position = "none",
#                                  axis.title.y = element_text(size = 14),
#                                  axis.title.x = element_text(size = 14),
#                                  strip.text.x = element_blank(),
#                                  strip.text.y = element_text(size = 12),
#                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
#   nrow = 4,
#   ncol = 1
# )
# legend<- get_legend(
#   p1 +
#     guides(fill = guide_legend(nrow = 1)) +
#     theme(legend.direction = "horizontal",
#           legend.position = "bottom",
#           legend.key.size = unit(0.5, 'cm'),
#           legend.text = element_text(size = 12),
#           legend.title = element_text(size = 12))
# )
#
#
# full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
#                        rel_heights = c(1, 0.05))
ggsave(filename = "rates.png",
       plot = full_plot,
       device="png", width=10,
       height=12, units="in", dpi=300)

###################################################
### Scenario 1 (prospective, non-PH, no truncation)
###################################################
dat <- readRDS("scenario_1.rds")
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stacking (fine)",
                                   stackG_medium = "Global stacking (medium)",
                                   stackG_coarse = "Global stacking (coarse)",
                                   stackL_fine = "Local stacking (fine)",
                                   stackL_medium = "Local stacking (medium)",
                                   stackL_coarse = "Local stacking (coarse)",
                                   survSL = "survSL",
                                   LTRCforests = "LTRC forests",
                                   coxph = "Linear Cox",
                                   gam = "Gen. additive Cox"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

dat_small <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stacking (fine)",
                                   stackG_medium = "Global stacking (ours)",
                                   stackG_coarse = "Global stacking (coarse)",
                                   stackL_fine = "Local stacking (fine)",
                                   stackL_medium = "Local stacking",
                                   stackL_coarse = "Local stacking (coarse)",
                                   survSL = "survSL",
                                   LTRCforests = "LTRC forests",
                                   coxph = "Linear Cox",
                                   gam = "Gen. additive Cox"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew")) %>%
  filter(Estimator %in% small_estimators)

pal <- greyscale(length(unique(dat$n_train)), start = 1, end = 0.4)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = Estimator,
             y = log(MSE_uni),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_50),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p1_small <- dat_small %>% mutate(type = "integrated") %>%
  ggplot(aes(x = Estimator,
             y = log(MSE_uni),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p2_small <- dat_small %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_50),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_75),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_90),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])
y_max_small <- max(layer_scales(p1_small)$y$range$range[2],
                   layer_scales(p2_small)$y$range$range[2])
y_min_small <- min(layer_scales(p1_small)$y$range$range[1],
                   layer_scales(p2_small)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 strip.text = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 +ylim(y_min, y_max) + theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  nrow = 4,
  ncol = 1
)
two_panel_plot <- ggarrange(
  p1_small + ylim(y_min_small, y_max_small) +theme(legend.position = "none",
                                                   axis.title.y = element_text(size = 14),
                                                   axis.title.x = element_blank(),
                                                   axis.text.x = element_blank(),
                                                   axis.ticks.x = element_blank(),
                                                   strip.text = element_text(size = 12),
                                                   plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2_small +  ylim(y_min_small, y_max_small) +theme(legend.position = "none",
                                                    axis.title.y = element_text(size = 14),
                                                    axis.title.x = element_text(size = 14),
                                                    strip.text.x = element_blank(),
                                                    strip.text.y = element_text(size = 12),
                                                    axis.text.x = element_text(size = 13),
                                                    plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  nrow = 2,
  ncol = 1
)
legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
)
legend_small<- get_legend(
  p1_small +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
)


full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.05))
ggsave(filename = "scenario_1.png",
       plot = full_plot,
       device="png", width=10,
       height=12, units="in", dpi=300)

main_plot <- plot_grid(two_panel_plot, legend_small, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.05))
ggsave(filename = "scenario_1_2panel.png",
       plot = main_plot,
       device="png", width=10,
       height=7, units="in", dpi=300)

################################################
### Scenario 2 (prospective, non-PH, truncation)
################################################
dat <- readRDS("scenario_2.rds")
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stacking (fine)",
                                   stackG_medium = "Global stacking (medium)",
                                   stackG_coarse = "Global stacking (coarse)",
                                   stackL_fine = "Local stacking (fine)",
                                   stackL_medium = "Local stacking (medium)",
                                   stackL_coarse = "Local stacking (coarse)",
                                   LTRCforests = "LTRC forests",
                                   coxph = "Linear Cox"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

dat_small <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stacking (fine)",
                                   stackG_medium = "Global stacking (ours)",
                                   stackG_coarse = "Global stacking (coarse)",
                                   stackL_fine = "Local stacking (fine)",
                                   stackL_medium = "Local stacking",
                                   stackL_coarse = "Local stacking (coarse)",
                                   LTRCforests = "LTRC forests",
                                   coxph = "Linear Cox"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew")) %>%
  filter(Estimator %in% small_estimators)

pal <- greyscale(length(unique(dat$n_train)), start = 1, end = 0.4)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = Estimator,
             y = log(MSE_uni),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_50),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p1_small <- dat_small %>% mutate(type = "integrated") %>%
  ggplot(aes(x = Estimator,
             y = log(MSE_uni),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p2_small <- dat_small %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_50),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_75),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_90),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])
y_max_small <- max(layer_scales(p1_small)$y$range$range[2],
                   layer_scales(p2_small)$y$range$range[2])
y_min_small <- min(layer_scales(p1_small)$y$range$range[1],
                   layer_scales(p2_small)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 strip.text = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 +ylim(y_min, y_max) + theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  nrow = 4,
  ncol = 1
)
two_panel_plot <- ggarrange(
  p1_small + ylim(y_min_small, y_max_small) +theme(legend.position = "none",
                                                   axis.title.y = element_text(size = 14),
                                                   axis.title.x = element_blank(),
                                                   axis.text.x = element_blank(),
                                                   axis.ticks.x = element_blank(),
                                                   strip.text = element_text(size = 12),
                                                   plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2_small +  ylim(y_min_small, y_max_small) +theme(legend.position = "none",
                                                    axis.title.y = element_text(size = 14),
                                                    axis.title.x = element_text(size = 14),
                                                    strip.text.x = element_blank(),
                                                    strip.text.y = element_text(size = 12),
                                                    axis.text.x = element_text(size = 13),
                                                    plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  nrow = 2,
  ncol = 1
)
legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
)
legend_small<- get_legend(
  p1_small +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
)

full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.05))
ggsave(filename = "scenario_2.png",
       plot = full_plot,
       device="png", width=10,
       height=12, units="in", dpi=300)

main_plot <- plot_grid(two_panel_plot, legend_small, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.05))
ggsave(filename = "scenario_2_2panel.png",
       plot = main_plot,
       device="png", width=10,
       height=7, units="in", dpi=300)

##############################
### Scenario 3 (retrospective)
##############################
dat <- readRDS("scenario_3.rds")
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stack (fine)",
                                   stackG_medium = "Global stack (medium)",
                                   stackG_coarse = "Global stack (coarse)",
                                   stackL_fine = "Local stack (fine)",
                                   stackL_medium = "Local stack (medium)",
                                   stackL_coarse = "Local stack (coarse)"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

pal <- greyscale(length(unique(dat$n_train)), start = 1, end = 0.4)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = Estimator,
             y = log(MSE_uni),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_50),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_75),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal,name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_90),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal,name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 strip.text = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 +ylim(y_min, y_max) + theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  nrow = 4,
  ncol = 1
)

legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
)
full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.05))
ggsave(filename = "scenario_3.png",
       plot = full_plot,
       device="png", width=10,
       height=12, units="in", dpi=300)

#####################################
### Scenario 4 (proportional hazards)
#####################################
dat <- readRDS("scenario_4.rds")
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stacking (fine)",
                                   stackG_medium = "Global stacking (medium)",
                                   stackG_coarse = "Global stacking (coarse)",
                                   stackL_fine = "Local stacking (fine)",
                                   stackL_medium = "Local stacking (medium)",
                                   stackL_coarse = "Local stacking (coarse)",
                                   LTRCforests = "LTRC forests",
                                   coxph = "Linear Cox"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

pal <- greyscale(length(unique(dat$n_train)), start = 1, end = 0.4)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = Estimator,
             y = log(MSE_uni),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_50),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_75),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal,name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_90),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal,name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 strip.text = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 +ylim(y_min, y_max) + theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  nrow = 4,
  ncol = 1
)

legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
)
full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.05))
ggsave(filename = "scenario_4.png",
       plot = full_plot,
       device="png", width=10,
       height=12, units="in", dpi=300)

##############################
### Scenario 5 (discrete time)
##############################
dat <- readRDS("scenario_5.rds")
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stacking",
                                   stackL_fine = "Local stacking",
                                   LTRCforests = "LTRC forests",
                                   coxph = "Linear Cox"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

nbins <- c(10,20, 50)
pal <- greyscale(length(unique(dat$n_train)), start = 1, end = 0.4)

for (nbin_j in nbins){
  p1 <- dat %>% mutate(type = "integrated") %>%
    ggplot(aes(x = Estimator,
               y = log(MSE_uni),
               position = factor(n_train),
               fill = factor(n_train))) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(type~dgp) +
    theme_bw() +
    ylab("log (MISE)") +
    xlab("Estimator") +
    scale_fill_manual(values = pal, name = "Sample size") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
    theme(text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.position = "none")

  p2 <- dat %>% mutate(type = "50th percentile") %>%
    ggplot(aes(x = Estimator,
               y = log(landmark_MSE_50),
               position = factor(n_train),
               fill = factor(n_train))) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(type~dgp) +
    theme_bw() +
    ylab("log (MSE)") +
    xlab("Estimator") +
    scale_fill_manual(values = pal, name = "Sample size") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
    theme(text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.position = "none")

  p3 <- dat %>% mutate(type = "75th percentile") %>%
    ggplot(aes(x = Estimator,
               y = log(landmark_MSE_75),
               position = factor(n_train),
               fill = factor(n_train))) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(type~dgp) +
    theme_bw() +
    ylab("log (MSE)") +
    xlab("Estimator") +
    scale_fill_manual(values = pal,name = "Sample size") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
    theme(text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.position = "none")

  p4 <- dat%>% mutate(type = "90th pecentile") %>%
    ggplot(aes(x = Estimator,
               y = log(landmark_MSE_90),
               position = factor(n_train),
               fill = factor(n_train))) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(type~dgp) +
    theme_bw() +
    ylab("log (MSE)") +
    xlab("Estimator") +
    scale_fill_manual(values = pal,name = "Sample size") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
    theme(text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.position = "none")

  y_max <- max(layer_scales(p1)$y$range$range[2],
               layer_scales(p2)$y$range$range[2],
               layer_scales(p3)$y$range$range[2],
               layer_scales(p4)$y$range$range[2])
  y_min <- min(layer_scales(p1)$y$range$range[1],
               layer_scales(p2)$y$range$range[1],
               layer_scales(p3)$y$range$range[1],
               layer_scales(p4)$y$range$range[1])

  four_panel_plot <- ggarrange(
    p1 + ylim(y_min, y_max) +theme(legend.position = "none",
                                   axis.title.y = element_text(size = 14),
                                   axis.title.x = element_blank(),
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   strip.text = element_text(size = 12),
                                   plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
    p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                   axis.title.y = element_text(size = 14),
                                   strip.text.x = element_blank(),
                                   strip.text.y = element_text(size = 12),
                                   axis.title.x = element_blank(),
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
    p3 +ylim(y_min, y_max) + theme(legend.position = "none",
                                   axis.title.y = element_text(size = 14),
                                   strip.text.x = element_blank(),
                                   strip.text.y = element_text(size = 12),
                                   axis.title.x = element_blank(),
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
    p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                   axis.title.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   strip.text.x = element_blank(),
                                   strip.text.y = element_text(size = 12),
                                   plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
    nrow = 4,
    ncol = 1
  )

  legend<- get_legend(
    p1 +
      guides(fill = guide_legend(nrow = 1)) +
      theme(legend.direction = "horizontal",
            legend.position = "bottom",
            legend.key.size = unit(0.5, 'cm'),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))
  )
  full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                         rel_heights = c(1, 0.05))
  ggsave(filename = paste0("scenario_5_", nbin_j, ".png"),
         plot = full_plot,
         device="png",
         width=10,
         height=12,
         units="in",
         dpi=300)
}
###################
### form comparison
###################
dat <- readRDS("form_comparison.rds")
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_PI = "Product integral",
                                   stackG_exp = "Exponential"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

pal <- greyscale(length(unique(dat$n_train)), start = 1, end = 0.4)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = Estimator,
             y = log(MSE_uni),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_50),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_75),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal,name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_90),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal,name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 strip.text = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 +ylim(y_min, y_max) + theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  nrow = 4,
  ncol = 1
)

legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
)
full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.05))
ggsave(filename = "form_comparison.png",
       plot = full_plot,
       device="png", width=10,
       height=12, units="in", dpi=300)

incompat <- dat %>% group_by(estimator, n_train) %>%
  summarize(mean(prop_incompat))


###################
### grid comparison
###################
dat <- readRDS(paste0(wd, "grid_comparison.rds"))

dat <- dat %>%
  mutate(Estimator = estimator,#recode_factor(estimator,
         #  stackG_fine_Y = "Shared grid (fine)",
         # stackG_medium_Y = "Shared grid (medium)",
         #stackG_coarse_Y = "Shared grid (coarse)",
         #stackG_fine_W = "Different grid (fine)",
         # stackG_medium_W = "Different grid (medium)",
         # stackG_coarse_W = "Different grid (coarse)"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

pal <- greyscale(length(unique(dat$n_train)), start = 1, end = 0.4)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = Estimator,
             y = log(MSE_uni),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_50),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_75),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_90),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 strip.text = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 +ylim(y_min, y_max) + theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  nrow = 4,
  ncol = 1
)

legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
)

full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.05))
ggsave(filename = "grid_comparison.png",
       plot = full_plot,
       device="png", width=10,
       height=12, units="in", dpi=300)

##########################
### approx grid comparison
##########################
dat <- readRDS(paste0(wd, "approx_grid_comparison.rds"))

dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_50 = "50 cutpoint grid",
                                   stackG_100 = "100 cutpoint grid",
                                   stackG_250 = "250 cutpoint grid",
                                   stackG_all = "All times grid"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

pal <- greyscale(length(unique(dat$n_train)), start = 1, end = 0.4)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = Estimator,
             y = log(MSE_uni))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_50))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_75))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_90))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 strip.text = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 +ylim(y_min, y_max) + theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  nrow = 4,
  ncol = 1
)

legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
)

full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.05))
ggsave(filename = "approx_grid_comparison.png",
       plot = full_plot,
       device="png", width=10,
       height=12, units="in", dpi=300)

#########################################
### grid stabilization (reviewer comment)
#########################################
wd <- "C:/Users/cwolo/Dropbox/UW/DISSERTATION/conditional_surv/scratch/revision_sims/jcgs_revision/grid_stabilization/"
setwd(wd)

# obs event grid instead of follow-up time grid
dat <- readRDS(paste0(wd, "grid_stabilization_092323.rds"))

summ <- dat %>% group_by(dgp, grid) %>%
  summarize(mean(diffs))

p1 <- dat %>%
  ggplot(aes(x = grid, y = diffs)) +
  geom_point() +
  facet_wrap(~ dgp)

dat <- dat %>%
  mutate(Estimator = estimator,
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

pal <- greyscale(length(unique(dat$n_train)), start = 1, end = 0.4)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = Estimator,
             y = log(MSE_uni),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Estimator") +
  # scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_50),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_75),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = Estimator,
             y = log(landmark_MSE_90),
             position = factor(n_train),
             fill = factor(n_train))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Estimator") +
  scale_fill_manual(values = pal, name = "Sample size") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")


y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 strip.text = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 +ylim(y_min, y_max) + theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 14),
                                 axis.title.x = element_text(size = 14),
                                 strip.text.x = element_blank(),
                                 strip.text.y = element_text(size = 12),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  nrow = 4,
  ncol = 1
)
legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
)

full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.05))
ggsave(filename = 'C:/Users/cwolo/Dropbox/UW/DISSERTATION/conditional_surv/notes/figures/jcgs_revision/prospective_notrunc_091123.png',
       plot = full_plot,
       device='png', width=10,
       height=12, units='in', dpi=300)



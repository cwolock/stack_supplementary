#!/usr/local/bin/Rscript

library("cowplot")
library("tidyverse")
library("grid")
library("gridExtra")
library("ggpubr")
library("egg")
library("RColorBrewer")

# colors
stackG_cols <- brewer.pal(n = 9, "Blues")[c(3,5,8)]
stackL_cols <- brewer.pal(n = 9, "Greens")[c(3,5,7)]
survSL_cols <- brewer.pal(n = 9, "YlOrRd")[5]
cox_cols <- brewer.pal(n = 12, "Paired")[11]

#############################
### prospective no truncation 
#############################
dat <- readRDS(paste0(wd, "prospective_notruncation"))
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
p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = factor(n_train), y = log(MSE_uni), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_50), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_75), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_90), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

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
                                 axis.title.y = element_text(size = 10),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 10),
                                 strip.text.x = element_blank(),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 +ylim(y_min, y_max) + theme(legend.position = "none",
                                 axis.title.y = element_text(size = 10),
                                 strip.text.x = element_blank(),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 10),
                                 strip.text.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.3, 0), "cm")),
  nrow = 4, 
  ncol = 1
)
legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 3)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8))
)
full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.1))
ggsave('prospective_notruncation.png', device='png', width=6,
       height=7.5, units='in', dpi=800)

##########################
### prospective truncation
##########################
dat <- readRDS(paste0(wd, "prospective_truncation.rds"))
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

pal <- c(stackG_cols, stackL_cols, cox_cols)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = factor(n_train), y = log(MSE_uni), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_50), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_75), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_90), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) + theme(legend.position = "none",
                                  axis.title.y = element_text(size = 10),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) + theme(legend.position = "none",
                                  axis.title.y = element_text(size = 10),
                                  strip.text.x = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 + ylim(y_min, y_max) + theme(legend.position = "none",
                                  axis.title.y = element_text(size = 10),
                                  strip.text.x = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) + theme(legend.position = "none",
                                  axis.title.y = element_text(size = 10),
                                  strip.text.x = element_blank(),
                                  plot.margin = unit(c(0, 0, 0.3, 0), "cm")),
  nrow = 4, 
  ncol = 1
)
legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 3)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8))
)
full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.1))
ggsave('prospective_truncation.png', device='png', width=6,
       height=7.5, units='in', dpi=800)

#################
### retrospective
#################
dat <- readRDS(paste0(wd, "/retrospective.rds"))

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
pal <- c(stackG_cols, stackL_cols)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = factor(n_train), y = log(MSE_uni), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_50), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_75), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_90), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) + 
    theme(legend.position = "none",
          axis.title.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) + 
    theme(legend.position = "none",
          axis.title.y = element_text(size = 10),
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 + ylim(y_min, y_max) + 
    theme(legend.position = "none",
          axis.title.y = element_text(size = 10),
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) + 
    theme(legend.position = "none",
          axis.title.y = element_text(size = 10),
          strip.text.x = element_blank(),
          plot.margin = unit(c(0, 0, 0.3, 0), "cm")),
  nrow = 4, 
  ncol = 1
)
legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 3)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8))
)
full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.1))
ggsave('retrospective.png', device='png', width=6,
       height=7.5, units='in', dpi=800)

########################
### proportional hazards
########################
dat <- readRDS(paste0(wd, "proportional_hazards.rds"))

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

pal <- c(stackG_cols, stackL_cols, cox_cols)

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = factor(n_train), y = log(MSE_uni), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MISE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p2 <- dat %>% mutate(type = "50th percentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_50), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p3 <- dat %>% mutate(type = "75th percentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_75), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

p4 <- dat%>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_90), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  scale_fill_manual(values = pal) + 
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") 

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) + theme(legend.position = "none",
                                  axis.title.y = element_text(size = 10),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 10),
                                 strip.text.x = element_blank(),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 10),
                                 strip.text.x = element_blank(),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) +theme(legend.position = "none",
                                 axis.title.y = element_text(size = 10),
                                 strip.text.x = element_blank(),
                                 plot.margin = unit(c(0, 0, 0.3, 0), "cm")),
  nrow = 4, 
  ncol = 1
)
legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 3)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8))
)
full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.1))
ggsave('proportional_hazards.png', device='png', width=6,
       height=7.5, units='in', dpi=800)


#################
### discrete time
#################
dat <- readRDS(paste0(wd, "discrete_time.rds"))
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_fine = "Global stacking (all times grid)",
                                   stackL_fine = "Local stacking (all times grid)",
                                   coxph = "Cox"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

pal <- c(stackG_cols[1], stackL_cols[1], cox_cols)
nbins <- c(10, 20, 50)

for (nbin_j in nbins){
  p1 <- dat %>% mutate(type = "integrated") %>%
    filter(nbin == nbin_j) %>%
    ggplot(aes(x = factor(n_train), y = log(MSE_uni), fill = Estimator)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(type~dgp) +
    theme_bw() +
    ylab("log (MISE)") +
    xlab("Sample size") +
    scale_fill_manual(values = pal) + 
    theme(text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = "outside") 
  
  p2 <- dat %>% mutate(type = "50th percentile") %>%
    filter(nbin == nbin_j) %>%
    ggplot(aes(x = factor(n_train), y = log(landmark_MSE_50), fill = Estimator)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(type~dgp) +
    theme_bw() +
    ylab("log (MSE)") +
    xlab("Sample size") +
    scale_fill_manual(values = pal) + 
    theme(text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = "outside") 
  
  p3 <- dat %>% mutate(type = "75th percentile") %>%
    filter(nbin == nbin_j) %>%
    ggplot(aes(x = factor(n_train), y = log(landmark_MSE_75), fill = Estimator)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(type~dgp) +
    theme_bw() +
    ylab("log (MSE)") +
    xlab("Sample size") +
    scale_fill_manual(values = pal) + 
    theme(text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = "outside") 
  
  p4 <- dat%>% mutate(type = "90th pecentile") %>%
    filter(nbin == nbin_j) %>%
    ggplot(aes(x = factor(n_train), y = log(landmark_MSE_90), fill = Estimator)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(type~dgp) +
    theme_bw() +
    ylab("log (MSE)") +
    xlab("Sample size") +
    scale_fill_manual(values = pal) + 
    theme(text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = "outside")
  
  y_max <- max(layer_scales(p1)$y$range$range[2],
               layer_scales(p2)$y$range$range[2],
               layer_scales(p3)$y$range$range[2],
               layer_scales(p4)$y$range$range[2])
  y_min <- min(layer_scales(p1)$y$range$range[1],
               layer_scales(p2)$y$range$range[1],
               layer_scales(p3)$y$range$range[1],
               layer_scales(p4)$y$range$range[1])
  
  four_panel_plot <- ggarrange(
    p1 +ylim(y_min, y_max) +  theme(legend.position = "none",
                                    axis.title.y = element_text(size = 10),
                                    axis.title.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
    p2 + ylim(y_min, y_max) + theme(legend.position = "none",
                                    axis.title.y = element_text(size = 10),
                                    strip.text.x = element_blank(),
                                    axis.title.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
    p3 + ylim(y_min, y_max) + theme(legend.position = "none",
                                    axis.title.y = element_text(size = 10),
                                    strip.text.x = element_blank(),
                                    axis.title.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
    p4 +ylim(y_min, y_max) +  theme(legend.position = "none",
                                    axis.title.y = element_text(size = 10),
                                    strip.text.x = element_blank(),
                                    plot.margin = unit(c(0, 0, 0.3, 0), "cm")),
    nrow = 4, 
    ncol = 1
  )
  legend<- get_legend(
    p1 +
      guides(fill = guide_legend(nrow = 1)) +
      theme(legend.direction = "horizontal",
            legend.position = "bottom",
            legend.key.size = unit(0.5, 'cm'),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8))
  )
  full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                         rel_heights = c(1, 0.1))
  ggsave(paste0('discrete_time_', nbin_j, '.png'), device='png', width=6,
         height=7.5, units='in', dpi=800)
}

###################
### form comparison
###################
dat <- readRDS(paste0(wd, "/form_comparison.rds"))
# time integrated uniformly
y_max <- max(log(dat$MSE_uni))
y_min <- min(log(dat$MSE_uni))
dat <- dat %>%
  mutate(Estimator = recode_factor(estimator,
                                   stackG_PI = "Product integral",
                                   stackG_exp = "Exponential"),
         dgp = recode_factor(dgp,
                             rightskew = "Right skew",
                             leftskew = "Left skew"))

p1 <- dat %>% mutate(type = "integrated") %>%
  ggplot(aes(x = factor(n_train), y = log(MSE_uni), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") + 
  scale_fill_manual(values = brewer.pal(n = 8, "Dark2")[c(2,3)])

p2 <- dat %>% mutate(type = "50th pecentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_50), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") + 
  scale_fill_manual(values = brewer.pal(n = 8, "Dark2")[c(2,3)])

p3 <- dat %>% mutate(type = "75th pecentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_50), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") + 
  scale_fill_manual(values = brewer.pal(n = 8, "Dark2")[c(2,3)])

p4 <- dat %>% mutate(type = "90th pecentile") %>%
  ggplot(aes(x = factor(n_train), y = log(landmark_MSE_50), fill = Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(type~dgp) +
  theme_bw() +
  ylab("log (MSE)") +
  xlab("Sample size") +
  theme(text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside") + 
  scale_fill_manual(values = brewer.pal(n = 8, "Dark2")[c(2,3)])

y_max <- max(layer_scales(p1)$y$range$range[2],
             layer_scales(p2)$y$range$range[2],
             layer_scales(p3)$y$range$range[2],
             layer_scales(p4)$y$range$range[2])
y_min <- min(layer_scales(p1)$y$range$range[1],
             layer_scales(p2)$y$range$range[1],
             layer_scales(p3)$y$range$range[1],
             layer_scales(p4)$y$range$range[1])

four_panel_plot <- ggarrange(
  p1 + ylim(y_min, y_max) + theme(legend.position = "none",
                                  axis.title.y = element_text(size = 10),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p2 + ylim(y_min, y_max) + theme(legend.position = "none",
                                  axis.title.y = element_text(size = 10),
                                  strip.text.x = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p3 + ylim(y_min, y_max) + theme(legend.position = "none",
                                  axis.title.y = element_text(size = 10),
                                  strip.text.x = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  plot.margin = unit(c(0, 0, 0.1, 0), "cm")),
  p4 + ylim(y_min, y_max) + theme(legend.position = "none",
                                  axis.title.y = element_text(size = 10),
                                  strip.text.x = element_blank(),
                                  plot.margin = unit(c(0, 0, 0.3, 0), "cm")),
  nrow = 4, 
  ncol = 1
)

legend<- get_legend(
  p1 +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8))
)
full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                       rel_heights = c(1, 0.1))
ggsave('form_comparison.png', device='png', width=6,
       height=7.5, units='in', dpi=800)
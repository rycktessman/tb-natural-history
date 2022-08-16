#run this script after running IMIS_combine

library(psych)
library(ppcor) #for partial rank correlations
library(tidyverse)
library(data.table)
library(ggpubr)
library(cowplot)
library(scales)
library(ggridges)

setwd("C:/Users/Tess/OneDrive - Johns Hopkins/TB/Natural History Modeling")
path_out <- "output/IMIS Nov2021 v4/"
source("code/calib_functions2.R")

mult_expand <- 1 #whether a_tx is a free parameter or not
RR_free <- 0 #whether a_r_s and a_p_s vary from a_r_m and a_p_m
spont_progress <- 0 #whether those who have spontaneously resolved can progress back to smear- symptom- TB
spont_prog <- 0.15 #annual probability of returning from resolved (if spont_progress==1)
smear_hist_calib <- 0 #whether to include historical targets on bacillary status over time
country <- "Cambodia"

#param names and labels for graphing
param_names <- c("Smear Progression", "Symptom Progression", 
                 "Smear Regression", "Symptom Regression",
                 "Smear Progression Rel. Risk", "Smear Regression Rel. Risk",
                 "Symptom Progression Rel. Risk", "Symptom Regression Rel. Risk",
                 "Treatment (Smear-)", "Treatment (Smear+)",
                 "TB Mortality (Smear-)", "TB Mortality (Smear+)",
                 "Spontaneous Resolution")
names(param_names) <- c("p_m", "p_s",
                        "r_m", "r_s",
                        "a_p_m", "a_r_m",
                        "a_p_s", "a_r_s",
                        "c_tx", "c_tx_m",
                        "m_tb", "m_tb_m",
                        "c_sp")

#targets
load("data/params_targets_cambodia.Rda")

#target names 
names_prev <- pull_targets("prev", targets_all, country)[[2]] 
names_prev[[4]] <- "Prevalence:Notifications"
names_hist_pos <- pull_targets("hist_pos", targets_all, country)[[2]]
names_hist_neg <- pull_targets("hist_neg", targets_all, country)[[2]]
names <- c(names_prev, names_hist_pos, names_hist_neg)
#multipliers when graphing each target
mults <- c(100, 100, 100, 1, 1000, 100, 100, 100, 100, 100)
names(mults) <- names(names)

#colors for different sensitivity analyses
colors <- c(hue_pal()(4), "orange") 
names(colors) <- c("Main analysis",  
                   "No 10-year historical mortality targets",
                   "Independent progression & regression rel. risks",
                   "15% annual progression from spontaneous resolution",
                   "Adding historical smear status targets")

#read additional files in for sensitivity analysis 
out_post <- read.csv(paste0("output/IMIS Nov2021 v4 cambodia/", "out_IMIS_combined", ".csv")) %>% 
  mutate(type=names(colors)[[1]])
out_post_no10 <- read.csv(paste0("output/IMIS Nov2021 v4 cambodia no10/", "out_IMIS_combined", ".csv")) %>% 
  mutate(type=names(colors)[[2]])
out_post_RRfree <- read.csv(paste0("output/IMIS Nov2021 v4 cambodia RRfree/", "out_IMIS_combined", ".csv")) %>% 
  mutate(type=names(colors)[[3]])
out_post_spontprog <- read.csv(paste0("output/IMIS Nov2021 v4 cambodia spontprog/", "out_IMIS_combined", ".csv")) %>% 
  mutate(type=names(colors)[[4]])
out_post_smearhist <- read.csv(paste0("output/IMIS Nov2021 v4 cambodia smearhist/", "out_IMIS_combined", ".csv")) %>% 
  mutate(type=names(colors)[[5]])

#combine all analyses together
out_post_all <- bind_rows(out_post, out_post_no10)
out_post_all <- bind_rows(out_post_all, out_post_RRfree)
out_post_all <- bind_rows(out_post_all, out_post_spontprog)
out_post_all <- bind_rows(out_post_all, out_post_smearhist)

#POSTERIOR DISTRIBUTIONS BY SENSITIVITY ANALYSIS
out_params <- out_post_all %>% mutate(m_tb_m=m_tb*a_m, c_tx_m=c_tx*a_tx) 
out_params <- out_params %>% select(names(param_names), type)
out_params <- out_params %>% 
  mutate(a_p_s=if_else(type==names(colors)[[3]], a_p_s, a_p_m),
         a_r_s=if_else(type==names(colors)[[3]], a_r_s, a_r_m))
out_params <- out_params %>% mutate(type=factor(type, levels=names(colors)))
param_plots <- list()
for(i in names(param_names)) {
  plot <- ggplot() +
    geom_density_ridges(data=out_params, aes_string(x=i, y="type", fill="type"), 
                        alpha=0.5, color="black") + 
    scale_x_continuous(expand=expansion(mult=0.01)) +
    scale_y_discrete(limits = rev(levels(out_params$type)), expand=c(0,0)) +
    scale_fill_manual(values=colors[levels(out_params$type)]) +
    labs(x="", y="", fill="") +
    ggtitle(param_names[[i]]) + theme_bw() + 
    theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(size=9, face="bold"),
          axis.text.x=element_text(size=8), axis.title=element_blank(),
          axis.text.y=element_blank(),
          plot.margin = margin(5.5, 10, 5.5, 5.5))
  if(i=="m_tb") {
    plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                      labels = scales::percent_format(accuracy=1),
                                      breaks=c(0,0.02,0.04,0.06,0.08,0.1))
  }
  if(i=="c_tx") {
    plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                      labels = scales::percent_format(accuracy=1),
                                      breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.25))
  }
  if(i=="a_p_m"|i=="a_p_s") {
    plot <- plot + scale_x_continuous(breaks=c(0, 2.5, 5, 7.5, 10))
  }
  if(!(i %in% c("a_r_m", "a_r_s", "a_p_m", "a_p_s", "m_tb", "c_tx"))) {
    plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                      labels = scales::percent_format(accuracy=1))
  }
  param_plots[[i]] <- plot
}
param_plots[[14]] <- as_ggplot(get_legend(plot + 
                                            theme(legend.position="right",
                                                  legend.margin=margin(l = 4, unit='cm'),
                                                  legend.box.background=element_rect(fill="transparent", 
                                                                                     color=NA),
                                                  legend.background=element_rect(fill="transparent", 
                                                                                 color=NA))))
fig <- ggarrange(plotlist=param_plots, align="hv", ncol=3, nrow=5, common.legend=F)
fig2 <- annotate_figure(fig, left=text_grob("Density", size=9, rot=90), 
                        bottom=text_grob("Monthly Values", size=9))
ggsave(fig2, filename=paste0(path_out, "param_post_SA_cambodia_ridges.jpg"), dpi=500, height=12, width=8)


#POSTERIOR MODEL OUTPUT DISTRIBUTIONS AND TARGETS BY SENSITIVITY ANALYSIS
out_targets <- out_post_all %>% select(names(names), type)
out_targets <- out_targets %>% mutate(type=factor(type, levels=names(colors)))
target_plots <- list()
for(i in names(names)) {
  plot <- ggplot() +
    geom_density_ridges(data=out_targets, aes_string(x=i, y="type", fill="type"), 
                        alpha=0.5, color="black", scale=1.25, rel_min_height=-0.01) + 
    geom_vline(xintercept=targets_all[[i]]) +
    geom_vline(xintercept=targets_all_lb[[i]], linetype="dashed") +
    geom_vline(xintercept=targets_all_ub[[i]], linetype="dashed") +
    scale_x_continuous(expand=expansion(mult=0.02)) +
    scale_y_discrete(limits=rev(levels(out_targets$type)), expand=c(0,0)) +
    scale_fill_manual(values=colors[levels(out_targets$type)]) +
    labs(x="", y="", fill="") +
    ggtitle(names[[i]]) + theme_bw() + 
    theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(size=9, face="bold"),
          axis.text.x=element_text(size=8), axis.title=element_blank(),
          axis.text.y=element_blank(),
          plot.margin = margin(5.5, 10, 5.5, 5.5))
  if(i!="pnr") {
    plot <- plot + scale_x_continuous(expand=expansion(mult=0.02), labels = scales::percent_format(accuracy=1)) 
  }
  target_plots[[i]] <- plot
}
target_plots[[11]] <- as_ggplot(get_legend(plot + 
                                            theme(legend.position="right",
                                                  legend.margin=margin(l = 4, unit='cm'),
                                                  legend.box.background=element_rect(fill="transparent", 
                                                                                     color=NA),
                                                  legend.background=element_rect(fill="transparent", 
                                                                                 color=NA))))
fig <- ggarrange(plotlist=target_plots, align="hv", ncol=3, nrow=4, common.legend=F)
ggsave(fig, filename=paste0(path_out, "targets_post_SA_cambodia_ridges.jpg"), dpi=500, height=10, width=8.5)

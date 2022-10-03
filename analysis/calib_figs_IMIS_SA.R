#run this script after running IMIS_combine

library(psych)
library(ppcor) #for partial rank correlations
library(tidyverse)
library(data.table)
library(ggpubr)
library(cowplot)
library(scales)
library(ggridges)
library(RColorBrewer)

setwd("~/GitHub/tb-natural-history/")
path_out <- "output/main/"
source("code/calib_functions.R")
country <- "Philippines"

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
load("data/mort_to_prev_IHME.Rda")
load(paste0("data/targets_", tolower(country), ".Rda"))

#alt targets used in sensitivity analyses
targets_all[["deaths_tb_ihme"]] <- mort_to_prev_mu[[country]]/100
targets_all_lb[["deaths_tb_ihme"]] <- mort_to_prev_lb[[country]]/100
targets_all_ub[["deaths_tb_ihme"]] <- mort_to_prev_ub[[country]]/100
targets_all[["prop_m_notif_alt"]] <- 0.5
targets_all_lb[["prop_m_notif_alt"]] <- 0.4
targets_all_ub[["prop_m_notif_alt"]] <- 0.6

#target names 
names_prev <- pull_targets("prev", targets_all, country)[[2]] 
if(country %in% c("Cambodia", "Philippines")) {
  names_prev[["pnr_m_all"]] <- "Prevalence:Notifications"
  names(targets_all)[names(targets_all)=="pnr_m_all"] <- "pnr"
  names(targets_all_lb)[names(targets_all_lb)=="pnr_m_all"] <- "pnr"
  names(targets_all_ub)[names(targets_all_ub)=="pnr_m_all"] <- "pnr"
} else {
  names(targets_all)[names(targets_all)=="pnr_all"] <- "pnr"
  names(targets_all_lb)[names(targets_all_lb)=="pnr_all"] <- "pnr"
  names(targets_all_ub)[names(targets_all_ub)=="pnr_all"] <- "pnr"
}
names(names_prev)[names_prev=="Prevalence:Notifications"] <- "pnr"
names_hist_pos <- pull_targets("hist_pos", targets_all, country)[[2]]
names_hist_neg <- pull_targets("hist_neg", targets_all, country)[[2]]
names <- c(names_prev, names_hist_pos, names_hist_neg)
#multipliers when graphing each target
mults <- c(100, 100, 100, 100, 1, 1000, 100, 100, 100, 100, 100)
names(mults) <- names(names)

#colors for different sensitivity analyses
colors <- brewer.pal(n=7, name="Dark2")
names(colors) <- c("Main analysis",  
                   "15% annual progression from spontaneous resolution",
                   "Constrained progression & regression rel. risks",
                   "Adding historical smear status targets",
                   "IHME mortality targets",
                   "50% smear-positive notifications target",
                   "No 10-year historical mortality targets")
analyses <- c("base", "spontprog", "rrconstrain", "smearhist",
              "ihmedeaths", "smearnotif50", "no10")
names(analyses) <- names(colors)

#read additional files in for sensitivity analysis 
out_post_all <- data.frame()
for(i in names(analyses)) {
  out_post <- read.csv(paste0("output/", tolower(country), "_", analyses[[i]], "/out_IMIS_combined.csv")) %>%
    mutate(type=i, lab=analyses[[i]])
  out_post_all <- bind_rows(out_post, out_post_all)
}
if(country %in% c("Cambodia", "Philippines")) {
  out_post_all <- out_post_all %>% mutate(pnr=pnr_m_all)
} else {
  out_post_all <- out_post_all %>% mutate(pnr=pnr_all)
}


#POSTERIOR DISTRIBUTIONS BY SENSITIVITY ANALYSIS
out_params <- out_post_all %>% mutate(m_tb_m=m_tb*a_m, c_tx_m=c_tx*a_tx) 
out_params <- out_params %>% select(names(param_names), type, lab)
out_params <- out_params %>% 
  mutate(a_p_s=if_else(lab=="rrconstrain", a_p_m, a_p_s),
         a_r_s=if_else(lab=="rrconstrain", a_r_m, a_r_s))
out_params <- out_params %>% mutate(type=factor(type, levels=names(colors)))

param_plots <- list()
for(i in names(param_names)) {
  #manually impute bandwidth for params where prior is hard to see
  if(i %in% c("p_m", "r_m")) {
    bw <- 0.005
  }
  if(i=="c_tx") {
    bw <- 0.005
  }
  if(i=="p_s") {
    bw <- 0.005
  }
  if(i=="r_s") {
    bw <-0.0085
  }
  if(i=="c_sp") {
    bw <- 0.015
  }
  if(i %in% c("a_p_m", "a_p_s")) {
    bw <- 0.3
  }
  if(i %in% c("a_r_m", "a_r_s")) {
    bw <- 0.025
  }
  if(i=="c_tx_m") {
    bw <- 0.01
  }
  if(i %in% c("m_tb", "m_tb_m")) {
    bw <- 0.001
  }
  plot <- ggplot() +
    geom_density_ridges(data=out_params, aes_string(x=i, y="type", fill="type"), 
                        alpha=0.5, color="black", bandwidth=bw) + 
    scale_x_continuous(expand=expansion(mult=0.01),
                       labels = scales::percent_format(accuracy=1)) +
    scale_y_discrete(limits = rev(levels(out_params$type)), expand=c(0,0)) +
    scale_fill_manual(values=colors[levels(out_params$type)]) +
    labs(x="", y="", fill="") +
    ggtitle(param_names[[i]]) + theme_bw() + 
    theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(size=9, face="bold"),
          axis.text.x=element_text(size=8), axis.title=element_blank(),
          axis.text.y=element_blank(),
          plot.margin = margin(5.5, 10, 5.5, 5.5))
  if(i %in% c("a_r_m","a_r_s")) {
    plot <- plot + scale_x_continuous(expand=expansion(mult=0.01))
  }
  if(i %in% c("a_p_m","a_p_s")) {
    plot <- plot + scale_x_continuous(expand=expansion(mult=0.01),
                                      limits=c(0, 10))
  }
  if(FALSE) {
    if(i=="m_tb") {
      plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                        labels = scales::percent_format(accuracy=1),
                                        breaks=c(0,0.02,0.04,0.06,0.08,0.1),
                                        limits=c(x_min, x_max))
    }
    if(i=="m_tb_m") {
      plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                        labels = scales::percent_format(accuracy=1),
                                        breaks=c(0,0.02,0.04,0.06,0.08,0.1),
                                        limits=c(x_min, x_max))
    }
    if(i=="c_tx") {
      plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                        labels = scales::percent_format(accuracy=1),
                                        breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.25),
                                        limits=c(x_min, x_max))
    }
    if(i=="c_tx_m") {
      plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                        labels = scales::percent_format(accuracy=1),
                                        breaks=c(0, 0.25, 0.5, 0.75, 1),
                                        limits=c(x_min, x_max))
    }
    if(i=="p_m"|i=="r_m") {
      plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                        labels = scales::percent_format(accuracy=1),
                                        breaks=c(0, 0.05, 0.1, 0.15, 0.2),
                                        limits=c(x_min, x_max))
    }
    if(i=="a_p_m"|i=="a_p_s") {
      plot <- plot + scale_x_continuous(breaks=c(0, 2.5, 5, 7.5, 10),
                                        limits=c(x_min, x_max))
    }
    if(!(i %in% c("a_p_m", "a_p_s", "m_tb", "m_tb_m", "c_tx", "c_tx_m", 
                  "p_m", "r_m"))) {
      plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                        labels = scales::percent_format(accuracy=1),
                                        limits=c(x_min, x_max))
    }
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
ggsave(fig2, filename=paste0(path_out, "param_post_SA_", tolower(country), "_ridges.jpg"), 
       dpi=500, height=12, width=8)


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
  if(i=="deaths_tb") {
    plot <- plot + geom_vline(xintercept=targets_all[["deaths_tb_ihme"]], color="darkgrey") +
      geom_vline(xintercept=targets_all_lb[["deaths_tb_ihme"]], color="darkgrey", linetype="dashed") +
      geom_vline(xintercept=targets_all_ub[["deaths_tb_ihme"]], color="darkgrey", linetype="dashed")
  }
  if(i=="prop_m_notif") {
    plot <- plot + geom_vline(xintercept=targets_all[["prop_m_notif_alt"]], color="darkgrey") +
      geom_vline(xintercept=targets_all_lb[["prop_m_notif_alt"]], color="darkgrey", linetype="dashed") +
      geom_vline(xintercept=targets_all_ub[["prop_m_notif_alt"]], color="darkgrey", linetype="dashed")
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
ggsave(fig, filename=paste0(path_out, "targets_post_SA_", tolower(country), "_ridges.jpg"), 
       dpi=500, height=12, width=8.5)

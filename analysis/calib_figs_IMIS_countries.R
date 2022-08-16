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
country <- "Philippines"

#param names and labels for graphing
param_names <- c("Smear Progression", "Symptom Progression", 
                 "Smear Regression", "Symptom Regression",
                 "Progression Rel. Risk", "Regression Rel. Risk",
                 "Treatment (Smear-)", "Treatment (Smear+)",
                 "TB Mortality (Smear-)", "TB Mortality (Smear+)",
                 "Spontaneous Resolution")
names(param_names) <- c("p_m", "p_s",
                        "r_m", "r_s",
                        "a_p_m", "a_r_m",
                        "c_tx", "c_tx_m",
                        "m_tb", "m_tb_m",
                        "c_sp")

#colors for countries
countries <- c("Philippines", "Vietnam", 
               "Bangladesh", "Nepal", "Cambodia")
colors_c <- c(hue_pal()(5), "darkgrey")
names(colors_c) <- c(countries, "Prior")

#load targets
target_means <- list()
target_lbs <- list()
target_ubs <- list()
for(i in countries) {
    print(i)
    if(i=="Philippines") {
        load("data/params_targets.Rda")
    } else {
        load(paste0("data/params_targets_", tolower(i), ".Rda"))
    }
    target_means[[i]] <- targets_all
    target_lbs[[i]] <- targets_all_lb
    target_ubs[[i]] <- targets_all_ub
}
target_means <- bind_rows(target_means, .id="country")
target_means <- target_means %>% 
    mutate(country=factor(country, levels=sort(countries, decreasing=T)), country_num=as.numeric(country),
           country_num_end=country_num+1)
target_means <- target_means %>%
    mutate(pnr=if_else(country %in% c("Philippines", "Cambodia"), pnr_m_all, pnr_all))
target_lbs <- bind_rows(target_lbs, .id="country")
target_lbs <- target_lbs %>% 
    mutate(country=factor(country, levels=sort(countries, decreasing=T)), country_num=as.numeric(country),
           country_num_end=country_num+1)
target_lbs <- target_lbs %>%
    mutate(pnr=if_else(country %in% c("Philippines", "Cambodia"), pnr_m_all, pnr_all))
target_ubs <- bind_rows(target_ubs, .id="country")
target_ubs <- target_ubs %>% 
    mutate(country=factor(country, levels=sort(countries, decreasing=T)), country_num=as.numeric(country),
           country_num_end=country_num+1)
target_ubs <- target_ubs %>%
    mutate(pnr=if_else(country %in% c("Philippines", "Cambodia"), pnr_m_all, pnr_all))

#target names 
names_prev <- pull_targets("prev", targets_all, "Philippines")[[2]] 
names_prev[[4]] <- "Prevalence:Notifications"
names(names_prev)[[4]] <- "pnr"
names_hist_pos <- pull_targets("hist_pos", targets_all, country)[[2]]
names_hist_neg <- pull_targets("hist_neg", targets_all, country)[[2]]
names <- c(names_prev, names_hist_pos, names_hist_neg)
#multipliers when graphing each target
mults <- c(100, 100, 100, 1, 1000, 100, 100, 100, 100, 100)
names(mults) <- names(names)

#options for sensitivity analyses
if(mult_expand==1) {
    params_calib_prev <- params_calib_prev_mult
    names_params_calib <- names_params_calib_mult
    priors_prev_lb <- priors_prev_mult_lb
    priors_prev_ub <- priors_prev_mult_ub
    params_fixed_hist <- params_fixed_hist_mult
}
if(RR_free==1) {
    priors_prev_lb <- priors_prev_lb_RRfree
    priors_prev_ub <- priors_prev_ub_RRfree
    params_calib_prev <- params_calib_prev_RRfree
    params_calib_hist <- params_calib_hist_RRfree
    names_params_calib <- names_params_calib_RRfree
    param_names <- param_names_RRfree
}
if(spont_progress==1) {
    params_fixed_prev[["p_c"]] <- 1-exp(log(1-spont_prog)/12) #monthly probability corresponding to annual probability of 3%
    params_fixed_hist[["p_c"]] <- params_fixed_prev[["p_c"]]
}
if(smear_hist_calib==1) {
    #sinding-larsen
    targets_all[["tb_smear_4yr_1"]] <- 61
    targets_all[["alive_4yr_1"]] <- 366
    #braeuning neisen
    targets_all[["tb_smear_4yr_2"]] <- 84
    targets_all[["alive_4yr_2"]] <- 166
    #griep
    targets_all[["tb_smear_4yr_3"]] <- 14
    targets_all[["alive_4yr_3"]] <- 57
    names_hist_pos <- c(names_hist_pos, "tb_smear_4yr"="Still Smear+ at 4 Yrs (% Alive)")
    names <- c(names_prev, names_hist_pos, names_hist_neg)
    mults <- c(mults, 100)
    names(mults) <- names(names)
}

#read files in (priors are the same across countries)
out_prior <- read.csv(paste0(path_out, "out_prior_combined", ".csv")) %>% 
    mutate(type="Prior", country="Prior")

#posteriors and performance for all countries
out_post_all <- list()
stats_rounds_all <- list()
for(i in countries) {
    print(i)
    out_post <- read.csv(paste0("output/IMIS Nov2021 v4 ", tolower(i), "/", "out_IMIS_combined", ".csv")) %>% 
        mutate(type="Posterior", country=i)
    stats_rounds <- read.csv(paste0("output/IMIS Nov2021 v4 ", tolower(i), "/", "stats_rounds_combined", ".csv")) %>% 
        mutate(country=i)
    out_post_all[[i]] <- out_post
    stats_rounds_all[[i]] <- stats_rounds
}
out_post_all <- bind_rows(out_post_all)
stats_rounds_all <- bind_rows(stats_rounds_all)

#FIGURE 2: PRIOR AND POSTERIOR DISTRIBUTIONS BY COUNTRY
out_params <- out_post_all %>% select(names(params_calib_prev), country)
out_params <- bind_rows(out_params, out_prior %>% select(names(params_calib_prev), country))
out_params <- out_params %>% mutate(m_tb_m=m_tb*a_m, c_tx_m=c_tx*a_tx) %>% select(-c(a_m, a_tx))
out_params <- out_params %>% mutate(country=factor(country, levels=c(sort(countries), "Prior")))
param_plots <- list()
for(i in names(param_names)) {
    #manually impute bandwidth for params where prior is hard to see
    if(i %in% c("p_m", "r_m")) {
        bw <- 0.01
    }
    if(i=="c_tx") {
        bw <- 0.01
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
    if(i %in% c("a_p_m")) {
        bw <- 0.6
    }
    if(i %in% c("a_r_m")) {
        bw <- 0.025
    }
    if(i=="c_tx_m") {
        bw <- 0.05
    }
    if(i %in% c("m_tb", "m_tb_m")) {
        bw <- 0.005
    }
    plot <- ggplot() +
        geom_density_ridges(data=out_params, aes_string(x=i, y="country", fill="country"), 
                            alpha=0.5, color="black", bandwidth=bw) + 
        scale_x_continuous(expand=expansion(mult=0.01)) +
        scale_y_discrete(limits = rev(levels(out_params$country)), expand=c(0,0)) +
        scale_fill_manual(values=colors_c[levels(out_params$country)]) +
        labs(x="", y="", fill="") +
        ggtitle(param_names[[i]]) + theme_bw() + 
        theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(size=9, face="bold"),
              axis.text=element_text(size=8), axis.title=element_blank(),
              plot.margin = margin(5.5, 10, 5.5, 5.5))
    if(i=="m_tb_m") {
        plot <- plot + scale_x_continuous(limits=c(0, 0.1), expand=expansion(mult=c(0, 0.001)),
                                          labels=scales::percent_format(accuracy=1),
                                          breaks=c(0,0.02,0.04,0.06,0.08,0.1)) 
    }
    if(i=="m_tb") {
        plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                          labels = scales::percent_format(accuracy=1),
                                          breaks=c(0,0.02,0.04,0.06,0.08,0.1))
    }
    if(i=="c_tx_m") {
        plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                          labels = scales::percent_format(accuracy=1),
                                          breaks=c(0, 0.25, 0.5, 0.75, 1))
    }
    if(i=="c_tx") {
        plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                          labels = scales::percent_format(accuracy=1),
                                          breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.25))
    }
    if(i=="a_p_m") {
        plot <- plot + scale_x_continuous(breaks=c(0, 2.5, 5, 7.5, 10))
    }
    if(!(i %in% c("a_r_m", "a_r_s", "a_p_m", "a_p_s", "m_tb", "m_tb_m", "c_tx", "c_tx_m"))) {
        plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                          labels = scales::percent_format(accuracy=1))
    }
    param_plots[[i]] <- plot
}
fig <- ggarrange(plotlist=param_plots, align="hv", ncol=4, nrow=3, common.legend=F)
fig2 <- annotate_figure(fig, left=text_grob("Density", size=9, rot=90), 
                        bottom=text_grob("Monthly Values", size=9))
ggsave(fig2, filename=paste0(path_out, "param_post_countries_ridges.jpg"), dpi=500, height=8.5, width=10)

#version for slides with different dimensions and no country axis
param_plots <- list()
for(i in names(param_names)) {
    #manually impute bandwidth for params where prior is hard to see
    if(i %in% c("p_m", "r_m")) {
        bw <- 0.01
    }
    if(i=="c_tx") {
        bw <- 0.01
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
    if(i %in% c("a_p_m")) {
        bw <- 0.6
    }
    if(i %in% c("a_r_m")) {
        bw <- 0.025
    }
    if(i=="c_tx_m") {
        bw <- 0.05
    }
    if(i %in% c("m_tb", "m_tb_m")) {
        bw <- 0.005
    }
    plot <- ggplot() +
        geom_density_ridges(data=out_params, aes_string(x=i, y="country", fill="country"), 
                            alpha=0.5, color="black", bandwidth=bw) + 
        scale_x_continuous(expand=expansion(mult=0.01)) +
        scale_y_discrete(limits = rev(levels(out_params$country)), expand=c(0,0)) +
        scale_fill_manual(values=colors_c[levels(out_params$country)]) +
        labs(x="", y="", fill="") +
        ggtitle(param_names[[i]]) + theme_bw() + 
        theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(size=9, face="bold"),
              axis.text.x=element_text(size=8), axis.title=element_blank(),
              axis.text.y=element_blank(), axis.ticks.y=element_blank(),
              plot.margin = margin(5.5, 10, 5.5, 5.5))
    if(i=="m_tb_m") {
        plot <- plot + scale_x_continuous(limits=c(0, 0.1), expand=expansion(mult=c(0, 0.001)),
                                          labels=scales::percent_format(accuracy=1),
                                          breaks=c(0,0.02,0.04,0.06,0.08,0.1)) 
    }
    if(i=="m_tb") {
        plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                          labels = scales::percent_format(accuracy=1),
                                          breaks=c(0,0.02,0.04,0.06,0.08,0.1))
    }
    if(i=="c_tx_m") {
        plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                          labels = scales::percent_format(accuracy=1),
                                          breaks=c(0, 0.25, 0.5, 0.75, 1))
    }
    if(i=="c_tx") {
        plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                          labels = scales::percent_format(accuracy=1),
                                          breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.25))
    }
    if(i=="a_p_m") {
        plot <- plot + scale_x_continuous(breaks=c(0, 2.5, 5, 7.5, 10))
    }
    if(!(i %in% c("a_r_m", "a_r_s", "a_p_m", "a_p_s", "m_tb", "m_tb_m", "c_tx", "c_tx_m"))) {
        plot <- plot + scale_x_continuous(expand=expansion(mult=0.001), 
                                          labels = scales::percent_format(accuracy=1))
    }
    param_plots[[i]] <- plot
}
param_plots[[length(param_plots)+1]] <- get_legend(plot + theme(legend.position="right"))
fig <- plot_grid(plotlist=param_plots, align="hv", nrow=2)
fig2 <- annotate_figure(fig, left=text_grob("Density", size=9, rot=90), 
                        bottom=text_grob("Monthly Values", size=9))
ggsave(fig2, filename=paste0(path_out, "param_post_countries_ridges_slides.jpg"), 
       dpi=500, height=4.5, width=11)

#generate means and CIs for each country (for text)
out_post_sum <- out_params %>% 
    mutate(r_m_ms=r_m*a_r_m,
           r_s_ms=r_s*a_r_m)
out_post_sum <- pivot_longer(out_post_sum, 
                             cols=names(out_post_sum)[names(out_post_sum)!="country"],
                             names_to="param", values_to="value")
out_post_sum <- out_post_sum %>% group_by(country, param) %>% 
    summarise(mean=mean(value), 
              lb=quantile(value, probs=0.025),
              ub=quantile(value, probs=0.975),
              q5=quantile(value, probs=0.05),
              q10=quantile(value, probs=0.1),
              q90=quantile(value, probs=0.9),
              q95=quantile(value, probs=0.95)) 
out_post_sum <- out_post_sum %>% arrange(param, country)
write.csv(out_post_sum, file=paste0(path_out, "calib_summary.csv"), row.names=F)

#TARGET MEANS/CIs AND POSTERIOR OUTPUT BY COUNTRY
out_targets <- out_post_all %>% mutate(pnr=if_else(country %in% c("Philippines", "Cambodia"), pnr_m_all, pnr_all)) #combine PNR into 1 column
out_targets <- out_targets %>% select(names(names), country)
out_targets <- out_targets %>% mutate(country=factor(country, levels=sort(countries, decreasing=T)),
                                      country_num=as.numeric(country))

target_plots <- list()
for(i in names(names)) {
    plot <- ggplot() +
        geom_density_ridges(data=out_targets, aes_string(x=i, y="country_num", fill="country"), 
                            alpha=0.5, color="black", scale=1.25, rel_min_height=-0.01) + 
        geom_segment(data=target_means, aes_string(x=i, xend=i, y="country_num", yend="country_num_end"), color="black") +
        geom_segment(data=target_lbs, aes_string(x=i, xend=i, y="country_num", yend="country_num_end"), color="black", linetype="dashed") +
        geom_segment(data=target_ubs, aes_string(x=i, xend=i, y="country_num", yend="country_num_end"), color="black", linetype="dashed") +
        #include horizontal black lines that go all the way across - next line does this:
        geom_segment(data=target_means, aes(y=country_num, yend=country_num), 
                     x=min(target_lbs[[i]]) - (max(target_ubs[[i]]) - min(target_lbs[[i]]))*0.2,
                     xend=max(target_ubs[[i]]) + (max(target_ubs[[i]]) - min(target_lbs[[i]]))*0.2) +
        scale_x_continuous(expand=expansion(mult=0.02)) +
        scale_y_discrete(limits=levels(out_targets$country), expand=c(0,0)) +
        scale_fill_manual(values=colors_c[levels(out_targets$country)]) +
        labs(x="", y="", fill="") +
        ggtitle(names[[i]]) + theme_bw() + 
        theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(size=9, face="bold"),
              axis.text=element_text(size=8), axis.title=element_blank(),
              plot.margin = margin(5.5, 10, 5.5, 5.5))
    if(i!="pnr") {
        plot <- plot + scale_x_continuous(expand=expansion(mult=0.02), labels = scales::percent_format(accuracy=1)) 
    }
    target_plots[[i]] <- plot
}
fig <- ggarrange(plotlist=target_plots, align="hv", ncol=3, nrow=4, common.legend=F)
ggsave(fig, filename=paste0(path_out, "targets_post_countries_ridges.jpg"), dpi=500, height=10, width=8.5)

#version for slides with different dimensions and no country axis
target_plots <- list()
for(i in names(names)) {
    plot <- ggplot() +
        geom_density_ridges(data=out_targets, aes_string(x=i, y="country_num", fill="country"), 
                            alpha=0.5, color="black", scale=1.25, rel_min_height=-0.01) + 
        geom_segment(data=target_means, aes_string(x=i, xend=i, y="country_num", yend="country_num_end"), color="black") +
        geom_segment(data=target_lbs, aes_string(x=i, xend=i, y="country_num", yend="country_num_end"), color="black", linetype="dashed") +
        geom_segment(data=target_ubs, aes_string(x=i, xend=i, y="country_num", yend="country_num_end"), color="black", linetype="dashed") +
        #include horizontal black lines that go all the way across - next line does this:
        geom_segment(data=target_means, aes(y=country_num, yend=country_num), 
                     x=min(target_lbs[[i]]) - (max(target_ubs[[i]]) - min(target_lbs[[i]]))*0.2,
                     xend=max(target_ubs[[i]]) + (max(target_ubs[[i]]) - min(target_lbs[[i]]))*0.2) +
        scale_x_continuous(expand=expansion(mult=0.02)) +
        scale_y_discrete(limits=levels(out_targets$country), expand=c(0,0)) +
        scale_fill_manual(values=colors_c[levels(out_targets$country)]) +
        labs(x="", y="", fill="") +
        ggtitle(names[[i]]) + theme_bw() + 
        theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(size=9, face="bold"),
              axis.text.x=element_text(size=8), axis.title=element_blank(),
              axis.text.y=element_blank(), axis.ticks.y=element_blank(),
              plot.margin = margin(5.5, 10, 5.5, 5.5))
    if(i!="pnr") {
        plot <- plot + scale_x_continuous(expand=expansion(mult=0.02), labels = scales::percent_format(accuracy=1)) 
    }
    target_plots[[i]] <- plot
}
legend <- get_legend(plot + theme(legend.position="bottom"))
fig <- plot_grid(plotlist=target_plots, align="hv", ncol=5, nrow=2)
fig2 <- plot_grid(fig, legend, align="hv", nrow=2, rel_heights=c(0.9, 0.1))
ggsave(fig2, filename=paste0(path_out, "targets_post_countries_ridges_slides.jpg"), 
       dpi=500, height=6, width=12)


#CALIBRATION PERFORMANCE BY CHAIN FOR EACH COUNTRY
#get chains on same scale
out_post_temp <- out_post_all %>% group_by(country) %>%
    mutate(chain2=as.numeric(as.factor(chain)))
ggplot(out_post_temp, aes(x=factor(chain2), y=log_like)) + 
    geom_boxplot() + facet_wrap(~country, nrow=5) +
    labs(x="Chain", y="Log Likelihood") +
    scale_x_discrete(breaks=as.character(c(10, 20, 30, 40, 50))) +
    theme_bw() + theme(panel.grid=element_blank())
ggsave(paste0(path_out, "log_like_box_countries.jpg"), dpi=500, height=10, width=7)

#PROXY CALCULATE ESS BY COUNTRY
out_post_temp <- unique(out_post_all)
out_post_temp <- out_post_temp %>% group_by(country) %>%
    mutate(weight=like/sum(like, na.rm=T))
ess <- out_post_temp %>% group_by(country) %>%
    summarise(ess=1/sum(weight^2))
ess

stats_rounds_temp <- stats_rounds_all %>% 
    filter(round==12 & !is.na(V4))
ess2 <- stats_rounds_temp %>% group_by(country) %>%
    summarise(ess2=sum(V4))
ess2

ess3 <- out_post_temp %>% group_by(country) %>%
    summarise(ess3=n())
ess3

ess_all <- left_join(ess, ess2, by="country")
ess_all <- left_join(ess_all, ess3, by="country")
names(ess_all) <- c("country",
                    "inverse_summed_sqr_wts",
                    "sum_ess_rounds",
                    "unique_samples")
write.csv(ess_all, paste0(path_out, "ess_proxy_calcs.csv"), row.names=F)

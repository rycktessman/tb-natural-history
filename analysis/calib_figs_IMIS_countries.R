#run this script after running IMIS_combine

library(psych)
library(ppcor) #for partial rank correlations
library(tidyverse)
library(data.table)
library(ggpubr)
library(cowplot)
library(scales)
library(ggridges)

setwd("~/GitHub/tb-natural-history")
path_out <- "output/main/"
source("code/calib_functions.R")
load("data/params_all.Rda")

RR_free <- 1 #4 free RR parameters in this version
spont_progress <- 0 #whether those who have spontaneously resolved can progress back to smear- symptom- TB
spont_prog <- 0.15 #what probability to use if spont_progress is 1
smear_hist_calib <- 0 #whether to include historical targets on bacillary status over time
deaths_targets <- "base" #"base", or "ihme" or use ihme targets
no_10yr_hist <- 0 #whether to include 10 year historical survival as calibration targets
smear_notif_override <- NA #NA, or an alt estimate +/- 10% (uniformly distributed)
RR_regress_recip <- 0 #if 1, fit one over the regression relative risk(s) instead of the regression relative risk(s)

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
countries <- c("Philippines", 
               "Vietnam", 
               "Bangladesh", 
               "Nepal", 
               "Cambodia"
               )
colors_c <- c(hue_pal()(5), "darkgrey")
names(colors_c) <- c(countries, "Prior")

#load targets
target_means <- list()
target_lbs <- list()
target_ubs <- list()
if(deaths_targets=="ihme") {
    load("data/mort_to_prev_ihme.Rda")
}
for(i in countries) {
    print(i)
    load(paste0("data/targets_", tolower(i), ".Rda"))
    target_means[[i]] <- targets_all
    target_lbs[[i]] <- targets_all_lb
    target_ubs[[i]] <- targets_all_ub
    if(deaths_targets=="ihme") {
        target_means[[i]][["deaths_tb_ihme"]] <- mort_to_prev_mu[[i]]/100
        target_lbs[[i]][["deaths_tb_ihme"]] <- mort_to_prev_lb[[i]]/100
        target_ubs[[i]][["deaths_tb_ihme"]] <- mort_to_prev_ub[[i]]/100
    }
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
names_prev[["pnr_m_all"]] <- "Prevalence:Notifications"
names(names_prev)[names_prev=="Prevalence:Notifications"] <- "pnr"
names_hist_pos <- pull_targets("hist_pos", targets_all, country)[[2]]
names_hist_neg <- pull_targets("hist_neg", targets_all, country)[[2]]
names <- c(names_prev, names_hist_pos, names_hist_neg)
#multipliers when graphing each target
mults <- c(100, 100, 100, 100, 1, 1000, 100, 100, 100, 100, 100)
names(mults) <- names(names)

#options for sensitivity analyses
scenario_lab <- "_base"
if(RR_free==1) {
    priors_prev_lb <- c(priors_prev_lb, "a_p_s"=priors_prev_lb$a_p_m, "a_r_s"=priors_prev_lb$a_r_m)
    priors_prev_ub <- c(priors_prev_ub, "a_p_s"=priors_prev_ub$a_p_m, "a_r_s"=priors_prev_ub$a_r_m)
    names_params_calib <- c(names_params_calib, "a_p_s"="Progression Multiplier (symptoms)",
                            "a_r_s"="Regression Multiplier (symptoms)")
    names_params_calib[["a_r_m"]] <- "Regression Multiplier (smear)"
    names_params_calib[["a_p_m"]] <- "Progression Multiplier (smear)"
    param_names <- c(param_names, "a_p_s"="Progression Rel. Risk (symptom)",
                     "a_r_s"="Regression Rel. Risk (symptom)")
    param_names[["a_r_m"]] <- "Regression Rel. Risk (smear)"
    param_names[["a_p_m"]] <- "Progression Rel. Risk (smear)"
}
if(RR_free==0) {
    scenario_lab <- "_rrconstrain"
}
if(spont_progress==1) {
    params_fixed_prev[["p_c"]] <- 1-exp(log(1-spont_prog)/12) #monthly probability corresponding to annual probability of 3%
    params_fixed_hist[["p_c"]] <- params_fixed_prev[["p_c"]]
    scenario_lab <- "_spontprog"
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
    scenario_lab <- "_smearhist"
    target_means[["tb_smear_4yr"]] <- 159/589
    target_lbs[["tb_smear_4yr"]] <- qbeta(0.025, 159, 589-159)
    target_ubs[["tb_smear_4yr"]] <- qbeta(0.975, 159, 589-159)
    
}
if(deaths_targets=="ihme") {
    scenario_lab <- "_ihmedeaths"
}
if(no_10yr_hist==1) {
    scenario_lab <- "_no10"
}
if(!is.na(smear_notif_override)) {
    scenario_lab <- paste0("_smearnotif", as.character(round(smear_notif_override*100)))
}
if(RR_regress_recip==1) {
    scenario_lab <- "_recipprior"
}


#read files in (priors are the same across countries)
out_prior <- read.csv(paste0(path_out, "out_priors_combined.csv")) %>% mutate(type="Prior", country="Prior")

if(scenario_lab=="_rrconstrain") {
    out_prior <- out_prior %>% mutate(a_r_s=a_r_m,
                                      a_p_s=a_p_m)
}
if(scenario_lab=="_recipprior") {
    out_prior <- read.csv(paste0(path_out, "out_priors_combined_recipprior.csv")) %>%
        mutate(type="Prior", country="Prior")
    out_prior <- out_prior %>% mutate(a_r_m=1/a_r_m_recip,
                                      a_r_s=1/a_r_s_recip)
}

#posteriors and performance for all countries
out_post_all <- list()
stats_rounds_all <- list()
ess_all <- list()
ess_each <- list()
for(i in countries) {
    path_in <- paste0("output/", tolower(i), scenario_lab, "/")
    print(i)
    out_post <- read.csv(paste0(path_in, "out_IMIS_combined", ".csv")) %>% 
        mutate(type="Posterior", country=i)
    stats_rounds <- read.csv(paste0(path_in, "stats_rounds_combined", ".csv")) %>% 
        mutate(country=i)
    ess <- read.csv(paste0(path_in, "ess_chains.csv"))
    out_post_all[[i]] <- out_post
    stats_rounds_all[[i]] <- stats_rounds
    ess_each[[i]] <- ess
    ess_all[[i]] <- sum(ess)
}
out_post_all <- bind_rows(out_post_all)
stats_rounds_all <- bind_rows(stats_rounds_all)
ess_all <- bind_rows(ess_all)
ess_each <- bind_cols(ess_each)
names(ess_each) <- countries

#FIGURE 2: PRIOR AND POSTERIOR DISTRIBUTIONS BY COUNTRY
if(scenario_lab=="_recipprior") {
    out_post_all <- out_post_all %>% mutate(a_r_m=1/a_r_m_recip, a_r_s=1/a_r_s_recip)
}
out_params <- out_post_all %>% select(names(names_params_calib), country)
out_params <- bind_rows(out_params, out_prior %>% select(names(names_params_calib), country))
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
    if(i %in% c("a_p_m", "a_p_s")) {
        bw <- 0.6
    }
    if(i %in% c("a_r_m", "a_r_s")) {
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
fig <- plot_grid(plotlist=param_plots, align="hv", ncol=3)
fig2 <- annotate_figure(fig, left=text_grob("Density", size=9, rot=90), 
                        bottom=text_grob("Monthly Values", size=9))
ggsave(fig2, filename=paste0(path_out, "param_post_countries_ridges", scenario_lab, ".jpg"), 
       dpi=500, height=12, width=10)

#version for slides with different dimensions and no country axis
if(FALSE) {
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
}

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
write.csv(out_post_sum, file=paste0(path_out, "calib_summary", scenario_lab, ".csv"), row.names=F)

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
    if(deaths_targets=="ihme" & i=="deaths_tb") {
        plot <- plot + 
            geom_segment(data=target_means, aes(x=deaths_tb_ihme, xend=deaths_tb_ihme, y=country_num, yend=country_num_end), color="red") +
            geom_segment(data=target_lbs, aes(x=deaths_tb_ihme, xend=deaths_tb_ihme, y=country_num, yend=country_num_end), color="red", linetype="dashed") +
            geom_segment(data=target_ubs, aes(x=deaths_tb_ihme, xend=deaths_tb_ihme, y=country_num, yend=country_num_end), color="red", linetype="dashed")
    }
    target_plots[[i]] <- plot
}
fig <- ggarrange(plotlist=target_plots, align="hv", ncol=3, nrow=4, common.legend=F)
ggsave(fig, filename=paste0(path_out, "targets_post_countries_ridges", scenario_lab, ".jpg"), 
       dpi=500, height=10, width=8.5)

#version for slides with different dimensions and no country axis
if(FALSE) {
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
}


#CALIBRATION PERFORMANCE BY CHAIN FOR EACH COUNTRY
#get chains on same scale
out_post_temp <- out_post_all %>% group_by(country) %>%
    mutate(chain2=as.numeric(as.factor(chain)))
ggplot(out_post_temp, aes(x=factor(chain2), y=log_like)) + 
    geom_boxplot() + facet_wrap(~country, nrow=5) +
    labs(x="Chain", y="Log Likelihood") +
    scale_x_discrete(breaks=as.character(c(10, 20, 30, 40, 50))) +
    theme_bw() + theme(panel.grid=element_blank())
ggsave(paste0(path_out, "log_like_box_countries", scenario_lab, ".jpg"), dpi=500, height=10, width=7)

#ESS by country
ess_all
colMeans(ess_each)

#Unique samples by country
unique(out_post_all) %>% group_by(country) %>% summarise(samples=n())
chain_sets <- unique(out_post_all) %>% group_by(country, chain) %>% summarise(samples=n())
chain_sets %>% group_by(country) %>% summarise(samples=mean(samples))

#improvement in likelihood over rounds
ggplot(stats_rounds_all, aes(x=round, y=V1, color=as.character(chain))) +
    geom_line() + facet_wrap(~country) +
    labs(x="", y="", color="") + 
    theme_bw() + theme(legend.position="none",
                       panel.grid=element_blank())

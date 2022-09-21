#generate figures from combined microsim output from MARCC
#these are the figures that have multiple countries combined
#for other figures with just 1 country, use microsim_figs_all_states
library(tidyverse)
library(data.table)
library(matrixStats)
library(Hmisc) #for weighted quantile function
library(ggpubr)
library(cowplot)
library(scales)
library(patchwork)
library(ggridges)
library(MetBrewer)

setwd("~/GitHub/tb-natural-history")
source("code/model_functions.R")
source("code/calib_functions.R")
path_out <- "output/main"
scenario_lab <- "base" #base, rrconstrain, spontprog, smearhist, ihmedeaths, no10, smearnotif50

#% by TB state for steady state calcs
countries <- c("Philippines", "Vietnam", "Nepal", "Cambodia", "Bangladesh")
steady_wts <- list()
sizes <- list() #prev survey counts
for(i in countries) {
  print(i)
  load(paste0("data/targets_", tolower(i), ".Rda"))
  wts <- targets_all[1:4]
  steady_wts[[i]] <- wts
  sizes[[i]] <- prev_cases
}

#keep consistent colors w/ calibration figs
colors_c <- c(hue_pal()(5))
names(colors_c) <- c("Philippines", "Vietnam", 
                     "Bangladesh", "Nepal", "Cambodia")

#also assign colors for each state
colors_s <- met.brewer("Veronese", n=5, type="discrete")[c(1,2,4,5)]
state_names <- c("Smear- Subclinical",
                 "Smear+ Subclinical", 
                 "Smear- Symptomatic", 
                 "Smear+ Symptomatic")
names(colors_s) <- state_names

#load files for each country and each starting population
out_post_all <- list()
props_all <- list()
times_all_all <- list()
for(i in countries) {
  print(i)
  path <- paste0("output/", tolower(i), "_", scenario_lab, "/")
  out_post <- read.csv(paste0(path, "out_IMIS_combined.csv"), stringsAsFactors=F)
  out_post_all[[i]] <- out_post
  props_c <- list()
  times_all_c <- list()
  for(start_pop in 1:4) {
    props <- read.csv(paste0(path, "props_comb", start_pop, ".csv"), stringsAsFactors=F)
    load(paste0(path, "times_all_comb", start_pop, ".Rda"))
    props_c[[start_pop]] <- props
    times_all_c[[start_pop]] <- times_all_comb
  }
  props_c <- bind_rows(props_c, .id="start_pop")
  times_all_c <- bind_rows(times_all_c, .id="start_pop")

  props_all[[i]] <- props_c
  times_all_all[[i]] <- times_all_c
}

#########################################################################
#PART 1: INDIVIDUAL TRAJECTORIES SUMMARY STATS
#CIs HERE INCORPORATE PARAMETER UNCERTAINTY, NOT STOCHASTIC UNCERTAINTY
#USED TO FORM HEATMAP/TABLE (heatmap is in separate script)
#########################################################################

#each country and initial state separately
times_out_all <- list()
path_out <- "output/main/"
for(i in countries) {
  times_out <- list()
  symptom_available <- T
  for(j in 1:4) {
    print(j)
    times_pop <- c()
    times_all_comb <- times_all_all[[i]] %>% filter(start_pop==j)
    props <- props_all[[i]] %>% filter(start_pop==j)
    #TOTAL TIME (unconditional)
    #mean and CI of total time smear- symptom-
    mu <- weighted.mean(times_all_comb[["mean"]][,"X1"], weights=times_all_comb[["mean"]]$weight)
    lb <- wtd.quantile(x=times_all_comb[["mean"]]$X1, probs=0.025, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    ub <- wtd.quantile(x=times_all_comb[["mean"]]$X1, probs=0.975, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
    #mean and CI of total time smear+ symptom- 
    mu <- weighted.mean(times_all_comb[["mean"]][,"X2"], weights=times_all_comb[["mean"]]$weight)
    lb <- wtd.quantile(x=times_all_comb[["mean"]]$X2, probs=0.025, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    ub <- wtd.quantile(x=times_all_comb[["mean"]]$X2, probs=0.975, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
    #mean and CI of total time smear- symptom+
    mu <- weighted.mean(times_all_comb[["mean"]][,"X3"], weights=times_all_comb[["mean"]]$weight)
    lb <- wtd.quantile(x=times_all_comb[["mean"]]$X3, probs=0.025, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    ub <- wtd.quantile(x=times_all_comb[["mean"]]$X3, probs=0.975, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
    #mean and CI of total time smear+ symptom+
    mu <- weighted.mean(times_all_comb[["mean"]][,"X4"], weights=times_all_comb[["mean"]]$weight)
    lb <- wtd.quantile(x=times_all_comb[["mean"]]$X4, probs=0.025, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    ub <- wtd.quantile(x=times_all_comb[["mean"]]$X4, probs=0.975, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
    #mean and CI of total time smear+ (any)
    mu <- weighted.mean(times_all_comb[["mean"]][,"smear"], weights=times_all_comb[["mean"]]$weight)
    lb <- wtd.quantile(x=times_all_comb[["mean"]]$smear, probs=0.025, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    ub <- wtd.quantile(x=times_all_comb[["mean"]]$smear, probs=0.975, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
    if(symptom_available==T) {
      #mean and CI of total time symptom+ (any)
      mu <- weighted.mean(times_all_comb[["mean"]][,"symptom"], weights=times_all_comb[["mean"]]$weight)
      lb <- wtd.quantile(x=times_all_comb[["mean"]]$symptom, probs=0.025, weight=times_all_comb[["mean"]]$weight,
                         normwt=T)
      ub <- wtd.quantile(x=times_all_comb[["mean"]]$symptom, probs=0.975, weight=times_all_comb[["mean"]]$weight,
                         normwt=T)
      times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
    } else {
      times_pop <- c(times_pop, "NA")
    }
    #mean and CI of total time w/ TB (any)
    mu <- weighted.mean(times_all_comb[["mean"]][,"tb_any"], weights=times_all_comb[["mean"]]$weight)
    lb <- wtd.quantile(x=times_all_comb[["mean"]]$tb_any, probs=0.025, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    ub <- wtd.quantile(x=times_all_comb[["mean"]]$tb_any, probs=0.975, weight=times_all_comb[["mean"]]$weight,
                       normwt=T)
    times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
    
    #PROPORTION EVER REACHING A GIVEN STATE
    #mean and CI of proportion reaching smear- symptom-
    mu <- 100*weighted.mean(props$X1, weights=props$weight)
    lb <- 100*wtd.quantile(x=props$X1, probs=0.025, weight=props$weight, normwt=T)
    ub <- 100*wtd.quantile(x=props$X1, probs=0.975, weight=props$weight, normwt=T)
    times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
    #mean and CI of proportion reaching smear+ symptom- 
    mu <- 100*weighted.mean(props[,"X2"], weights=props$weight)
    lb <- 100*wtd.quantile(x=props$X2, probs=0.025, weight=props$weight, normwt=T)
    ub <- 100*wtd.quantile(x=props$X2, probs=0.975, weight=props$weight, normwt=T)
    times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
    #mean and CI of proportion reaching smear- symptom+
    mu <- 100*weighted.mean(props[,"X3"], weights=props$weight)
    lb <- 100*wtd.quantile(x=props$X3, probs=0.025, weight=props$weight, normwt=T)
    ub <- 100*wtd.quantile(x=props$X3, probs=0.975, weight=props$weight, normwt=T)
    times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
    #mean and CI of proportion reaching smear+ symptom+
    mu <- 100*weighted.mean(props[,"X4"], weights=props$weight)
    lb <- 100*wtd.quantile(x=props$X4, probs=0.025, weight=props$weight, normwt=T)
    ub <- 100*wtd.quantile(x=props$X4, probs=0.975, weight=props$weight, normwt=T)
    times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
    #mean and CI of proportion reaching smear+ (any)
    mu <- 100*weighted.mean(props[,"smear"], weights=props$weight)
    lb <- 100*wtd.quantile(x=props$smear, probs=0.025, weight=props$weight, normwt=T)
    ub <- 100*wtd.quantile(x=props$smear, probs=0.975, weight=props$weight, normwt=T)
    times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
    #mean and CI of proportion reaching symptom+ (any)
    if(symptom_available==T) {
      mu <- 100*weighted.mean(props[,"symptom"], weights=props$weight)
      lb <- 100*wtd.quantile(x=props$symptom, probs=0.025, weight=props$weight, normwt=T)
      ub <- 100*wtd.quantile(x=props$symptom, probs=0.975, weight=props$weight, normwt=T)
      times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
    } else {
      times_pop <- c(times_pop, "NA")
    }
    #mean and CI of proportion reaching spontaneous resolution
    mu <- 100*weighted.mean(props[,"X5"], weights=props$weight)
    lb <- 100*wtd.quantile(x=props$X5, probs=0.025, weight=props$weight, normwt=T)
    ub <- 100*wtd.quantile(x=props$X5, probs=0.975, weight=props$weight, normwt=T)
    times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
    #mean and CI of proportion reaching treatment
    mu <- 100*weighted.mean(props[,"X6"], weights=props$weight)
    lb <- 100*wtd.quantile(x=props$X6, probs=0.025, weight=props$weight, normwt=T)
    ub <- 100*wtd.quantile(x=props$X6, probs=0.975, weight=props$weight, normwt=T)
    times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
    #mean and CI of proportion reaching death (TB only)
    mu <- 100*weighted.mean(props[,"X7"], weights=props$weight)
    lb <- 100*wtd.quantile(x=props[,"X7"], probs=0.025, weight=props$weight, normwt=T)
    ub <- 100*wtd.quantile(x=props[,"X7"], probs=0.975, weight=props$weight, normwt=T)
    times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
    times_out[[j]] <- times_pop
  }
  times_out <- bind_cols(times_out)
  times_out_all[[i]] <- times_out
  write.csv(times_out, file=paste0(path_out, "times_out_", 
                                   tolower(i), "_", scenario_lab, ".csv"), row.names=F)
}

#steady state distribution (each country separately)
times_pop_steady_all <- list()
for(i in countries) {
  print(i)
  symptom_available <- T
  times_pop_steady <- list()
  
  times_all_steady <- bind_rows(lapply(1:4, function(x) times_all_all[[i]] %>% filter(start_pop==x) %>% pull(mean)),
                                .id="start_pop")
  times_all_steady <- times_all_steady %>% mutate(start_pop=as.numeric(start_pop))
  #exclude parameter sets if sims didn't finish for all 4 starting populations
  times_all_steady <- times_all_steady %>% group_by(id) %>% mutate(id_count=n())
  times_all_steady <- times_all_steady %>% filter(id_count==4)
  #apply state population weights
  times_all_steady <- times_all_steady %>% mutate(weight_steady=unname(unlist(steady_wts[[i]][start_pop])))
  #calculate averages across all states (weighted by state dist) for each param set
  times_all_steady <- times_all_steady %>% group_by(id) %>% select(-c(id_count, start_pop)) %>%
    summarise_all(~weighted.mean(., w=weight_steady))
  #TOTAL TIME (unconditional)
  #mean and CI of total time smear- symptom-
  mu <- weighted.mean(times_all_steady$X1, weights=times_all_steady$weight)
  lb <- wtd.quantile(x=times_all_steady$X1, probs=0.025, weight=times_all_steady$weight,
                     normwt=T)
  ub <- wtd.quantile(x=times_all_steady$X1, probs=0.975, weight=times_all_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear+ symptom- 
  mu <- weighted.mean(times_all_steady$X2, weights=times_all_steady$weight)
  lb <- wtd.quantile(x=times_all_steady$X2, probs=0.025, weight=times_all_steady$weight,
                     normwt=T)
  ub <- wtd.quantile(x=times_all_steady$X2, probs=0.975, weight=times_all_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear- symptom+
  mu <- weighted.mean(times_all_steady$X3, weights=times_all_steady$weight)
  lb <- wtd.quantile(x=times_all_steady$X3, probs=0.025, weight=times_all_steady$weight,
                     normwt=T)
  ub <- wtd.quantile(x=times_all_steady$X3, probs=0.975, weight=times_all_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear+ symptom+
  mu <- weighted.mean(times_all_steady$X4, weights=times_all_steady$weight)
  lb <- wtd.quantile(x=times_all_steady$X4, probs=0.025, weight=times_all_steady$weight,
                     normwt=T)
  ub <- wtd.quantile(x=times_all_steady$X4, probs=0.975, weight=times_all_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear+ (any)
  mu <- weighted.mean(times_all_steady$smear, weights=times_all_steady$weight)
  lb <- wtd.quantile(x=times_all_steady$smear, probs=0.025, weight=times_all_steady$weight,
                     normwt=T)
  ub <- wtd.quantile(x=times_all_steady$smear, probs=0.975, weight=times_all_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time symptom+ (any)
  if(symptom_available) {
    mu <- weighted.mean(times_all_steady$symptom, weights=times_all_steady$weight)
    lb <- wtd.quantile(x=times_all_steady$symptom, probs=0.025, weight=times_all_steady$weight,
                       normwt=T)
    ub <- wtd.quantile(x=times_all_steady$symptom, probs=0.975, weight=times_all_steady$weight,
                       normwt=T)
    times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  } else {
    times_pop_steady <- c(times_pop_steady, "NA")
  }
  #mean and CI of total time w/ TB (any)
  mu <- weighted.mean(times_all_steady$tb_any, weights=times_all_steady$weight)
  lb <- wtd.quantile(x=times_all_steady$tb_any, probs=0.025, weight=times_all_steady$weight,
                     normwt=T)
  ub <- wtd.quantile(x=times_all_steady$tb_any, probs=0.975, weight=times_all_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  
  #proportion ever reaching a state
  props_steady <- props_all[[i]]
  props_steady <- props_steady %>% mutate(start_pop=as.numeric(start_pop))
  #exclude parameter sets if sims didn't finish for all 4 starting populations
  props_steady <- props_steady %>% group_by(id) %>% mutate(id_count=n())
  props_steady <- props_steady %>% filter(id_count==4)
  #apply state population weights
  props_steady <- props_steady %>% mutate(weight_steady=unname(unlist(steady_wts[[i]][start_pop])))
  #calculate averages across all states (weighted by state dist) for each param set
  props_steady <- props_steady %>% group_by(id) %>% select(-c(id_count, start_pop)) %>%
    summarise_all(~weighted.mean(., w=weight_steady))
  #mean and CI of total time smear- symptom-
  mu <- 100*weighted.mean(props_steady$X1, weights=props_steady$weight)
  lb <- 100*wtd.quantile(x=props_steady$X1, probs=0.025, weight=props_steady$weight,
                     normwt=T)
  ub <- 100*wtd.quantile(x=props_steady$X1, probs=0.975, weight=props_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of total time smear+ symptom- 
  mu <- 100*weighted.mean(props_steady$X2, weights=props_steady$weight)
  lb <- 100*wtd.quantile(x=props_steady$X2, probs=0.025, weight=props_steady$weight,
                     normwt=T)
  ub <- 100*wtd.quantile(x=props_steady$X2, probs=0.975, weight=props_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of total time smear- symptom+
  mu <- 100*weighted.mean(props_steady$X3, weights=props_steady$weight)
  lb <- 100*wtd.quantile(x=props_steady$X3, probs=0.025, weight=props_steady$weight,
                     normwt=T)
  ub <- 100*wtd.quantile(x=props_steady$X3, probs=0.975, weight=props_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of total time smear+ symptom+
  mu <- 100*weighted.mean(props_steady$X4, weights=props_steady$weight)
  lb <- 100*wtd.quantile(x=props_steady$X4, probs=0.025, weight=props_steady$weight,
                     normwt=T)
  ub <- 100*wtd.quantile(x=props_steady$X4, probs=0.975, weight=props_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of total time smear+ (any)
  mu <- 100*weighted.mean(props_steady$smear, weights=props_steady$weight)
  lb <- 100*wtd.quantile(x=props_steady$smear, probs=0.025, weight=props_steady$weight,
                     normwt=T)
  ub <- 100*wtd.quantile(x=props_steady$smear, probs=0.975, weight=props_steady$weight,
                     normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of total time symptom+ (any)
  if(symptom_available) {
    mu <- 100*weighted.mean(props_steady$symptom, weights=props_steady$weight)
    lb <- 100*wtd.quantile(x=props_steady$symptom, probs=0.025, weight=props_steady$weight,
                         normwt=T)
    ub <- 100*wtd.quantile(x=props_steady$symptom, probs=0.975, weight=props_steady$weight,
                         normwt=T)
    times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  } else {
    times_pop_steady <- c(times_pop_steady, "NA")
  }
  #mean and CI of proportion reaching spontaneous resolution
  mu <- 100*weighted.mean(props_steady$X5, weights=props_steady$weight)
  lb <- 100*wtd.quantile(x=props_steady$X5, probs=0.025, weight=props_steady$weight,
                         normwt=T)
  ub <- 100*wtd.quantile(x=props_steady$X5, probs=0.975, weight=props_steady$weight,
                         normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching treatment
  mu <- 100*weighted.mean(props_steady$X6, weights=props_steady$weight)
  lb <- 100*wtd.quantile(x=props_steady$X6, probs=0.025, weight=props_steady$weight,
                         normwt=T)
  ub <- 100*wtd.quantile(x=props_steady$X6, probs=0.975, weight=props_steady$weight,
                         normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching death (TB only)
  mu <- 100*weighted.mean(props_steady$X7, weights=props_steady$weight)
  lb <- 100*wtd.quantile(x=props_steady$X7, probs=0.025, weight=props_steady$weight, normwt=T)
  ub <- 100*wtd.quantile(x=props_steady$X7, probs=0.975, weight=props_steady$weight, normwt=T)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  times_pop_steady_all[[i]] <- unlist(times_pop_steady)
}
times_pop_steady_all <- bind_cols(times_pop_steady_all)
write.csv(times_pop_steady_all, 
          file=paste0(path_out, "times_out_steady", "_", scenario_lab ,".csv"), row.names=F)

#countries pooled together
symptom_available <- T
times_pool_out <- list()
#for each state
for(j in 1:4) {
  print(j)
  times_pop <- c()
  #sample from times for each country and pool together
  times_pool <- sapply(countries, function(x)
    (times_all_all[[x]] %>% filter(start_pop==j))[["mean"]], simplify=F)
  times_pool <- sapply(countries, function(x)
    times_pool[[x]][sample(1:nrow(times_pool[[x]]), 50000, replace=T,
           prob=times_pool[[x]][["weight"]]), ], simplify=F)
  times_pool <- bind_rows(times_pool, .id="country")
  #sample from props for each country and pool together
  props_pool <- sapply(countries, function(x)
    props_all[[x]] %>% filter(start_pop==j), simplify=F)
  props_pool <- sapply(countries, function(x)
    props_pool[[x]][sample(1:nrow(props_pool[[x]]), 50000, replace=T,
                           prob=props_pool[[x]][["weight"]]), ], simplify=F)
  props_pool <- bind_rows(props_pool, .id="country")
  #TOTAL TIME (unconditional)
  #mean and CI of total time smear- symptom-
  mu <- mean(times_pool[,"X1"])
  lb <- quantile(x=times_pool$X1, probs=0.025)
  ub <- quantile(x=times_pool$X1, probs=0.975)
  times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear+ symptom- 
  mu <- mean(times_pool[,"X2"])
  lb <- quantile(x=times_pool$X2, probs=0.025)
  ub <- quantile(x=times_pool$X2, probs=0.975)
  times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear- symptom+
  mu <- mean(times_pool[,"X3"])
  lb <- quantile(x=times_pool$X3, probs=0.025)
  ub <- quantile(x=times_pool$X3, probs=0.975)
  times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear+ symptom+
  mu <- mean(times_pool[,"X4"])
  lb <- quantile(x=times_pool$X4, probs=0.025)
  ub <- quantile(x=times_pool$X4, probs=0.975)
  times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear+ (any)
  mu <- mean(times_pool[,"smear"])
  lb <- quantile(x=times_pool$smear, probs=0.025)
  ub <- quantile(x=times_pool$smear, probs=0.975)
  times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  if(symptom_available==T) {
    #mean and CI of total time symptom+ (any)
    mu <- mean(times_pool[,"symptom"], na.rm=T)
    lb <- quantile(x=times_pool$symptom, probs=0.025, na.rm=T)
    ub <- quantile(x=times_pool$symptom, probs=0.975, na.rm=T)
    times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  } else {
    times_pop <- c(times_pop, "NA")
  }
  #mean and CI of total time w/ TB (any)
  mu <- mean(times_pool[,"tb_any"])
  lb <- quantile(x=times_pool$tb_any, probs=0.025)
  ub <- quantile(x=times_pool$tb_any, probs=0.975)
  times_pop <- c(times_pop, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  
  #PROPORTION EVER REACHING A GIVEN STATE
  #mean and CI of proportion reaching smear- symptom-
  mu <- 100*mean(props_pool$X1)
  lb <- 100*quantile(x=props_pool$X1, probs=0.025)
  ub <- 100*quantile(x=props_pool$X1, probs=0.975)
  times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching smear+ symptom- 
  mu <- 100*mean(props_pool[,"X2"])
  lb <- 100*quantile(x=props_pool$X2, probs=0.025)
  ub <- 100*quantile(x=props_pool$X2, probs=0.975)
  times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching smear- symptom+
  mu <- 100*mean(props_pool[,"X3"])
  lb <- 100*quantile(x=props_pool$X3, probs=0.025)
  ub <- 100*quantile(x=props_pool$X3, probs=0.975)
  times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching smear+ symptom+
  mu <- 100*mean(props_pool[,"X4"])
  lb <- 100*quantile(x=props_pool$X4, probs=0.025)
  ub <- 100*quantile(x=props_pool$X4, probs=0.975)
  times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching smear+ (any)
  mu <- 100*mean(props_pool[,"smear"])
  lb <- 100*quantile(x=props_pool$smear, probs=0.025)
  ub <- 100*quantile(x=props_pool$smear, probs=0.975)
  times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching symptom+ (any)
  if(symptom_available==T) {
    mu <- 100*mean(props_pool[,"symptom"], na.rm=T)
    lb <- 100*quantile(x=props_pool$symptom, probs=0.025, na.rm=T)
    ub <- 100*quantile(x=props_pool$symptom, probs=0.975, na.rm=T)
    times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  } else {
    times_pop <- c(times_pop, "NA")
  }
  #mean and CI of proportion reaching spontaneous resolution
  mu <- 100*mean(props_pool[,"X5"])
  lb <- 100*quantile(x=props_pool$X5, probs=0.025)
  ub <- 100*quantile(x=props_pool$X5, probs=0.975)
  times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching treatment
  mu <- 100*mean(props_pool[,"X6"])
  lb <- 100*quantile(x=props_pool$X6, probs=0.025)
  ub <- 100*quantile(x=props_pool$X6, probs=0.975)
  times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching death (TB only)
  mu <- 100*mean(props_pool[,"X7"])
  lb <- 100*quantile(x=props_pool[,"X7"], probs=0.025)
  ub <- 100*quantile(x=props_pool[,"X7"], probs=0.975)
  times_pop <- c(times_pop, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  times_pool_out[[j]] <- times_pop
}
#repeat for steady state
if(TRUE) {
  times_pop_steady <- list()
  times_pool_steady <- sapply(countries, function(x)
    cbind("start_pop"=times_all_all[[x]][["start_pop"]], 
              times_all_all[[x]][["mean"]]), simplify=F)
  #exclude parameter sets if sims didn't finish for all 4 starting populations
  times_pool_steady <- sapply(countries, function(x)
    (times_pool_steady[[x]]) %>% group_by(id) %>% mutate(id_count=n()), simplify=F)
  times_pool_steady <- sapply(countries, function(x)
    (times_pool_steady[[x]]) %>% filter(id_count==4) %>%
      mutate(id=as.character(id)), simplify=F)
  #assign unique ID for each country row
  times_pool_steady <- sapply(countries, function(x)
    times_pool_steady[[x]] %>% ungroup() %>% group_by(start_pop) %>%
      mutate(id_unique=rep(1:n())),
    simplify=F)
  #sample 50k parameter sets for each country
  ids <- sapply(countries, function(x)
    sample(times_pool_steady[[x]] %>% filter(start_pop==1) %>% 
             pull("id_unique"), 50000, 
           replace=T,prob=times_pool_steady[[x]] %>% 
             filter(start_pop==1) %>% pull(weight)),
    simplify=F)
  row_nums <- sapply(countries, function(x)
    c(match(ids[[x]], times_pool_steady[[x]] %>% filter(start_pop==1) %>%
            pull("id_unique")),
      nrow(times_pool_steady[[x]])/4 + 
        match(ids[[x]], times_pool_steady[[x]] %>% filter(start_pop==1) %>%
              pull("id_unique")),
      2*nrow(times_pool_steady[[x]])/4 + 
      match(ids[[x]], times_pool_steady[[x]] %>% filter(start_pop==1) %>%
              pull("id_unique")),
      3*nrow(times_pool_steady[[x]])/4 + 
      match(ids[[x]], times_pool_steady[[x]] %>% filter(start_pop==1) %>%
              pull("id_unique"))),
    simplify=F)
  times_pool_steady <- sapply(countries, function(x)
    times_pool_steady[[x]][row_nums[[x]], ], simplify=F)
  #merge in pop weights
  times_pool_steady <- sapply(countries, function(x)
    times_pool_steady[[x]] %>% 
      mutate(start_pop=as.numeric(start_pop)),
    simplify=F)
  times_pool_steady <- sapply(countries, function(x)
    times_pool_steady[[x]] %>%
      mutate(weight_steady=unname(unlist(steady_wts[[x]][start_pop]))),
    simplify=F) 
  #bind rows together 
  times_pool_steady <- bind_rows(times_pool_steady, .id="country")
  
  #calculate averages across all states (weighted by state dist) for each param set
  symptom_steady <- times_pool_steady %>%
    group_by(country, id) %>%
    summarise(symptom=weighted.mean(symptom, w=weight_steady))
  times_pool_steady <- times_pool_steady %>% group_by(country, id) %>% 
    select(-c(id_count, start_pop, id_unique, weight, symptom)) %>%
    summarise_all(~weighted.mean(., w=weight_steady, na.rm=T))
  times_pool_steady <- cbind(times_pool_steady, "symptom"=(symptom_steady %>% pull(symptom)))
  #TOTAL TIME (unconditional)
  #mean and CI of total time smear- symptom-
  mu <- mean(times_pool_steady$X1)
  lb <- quantile(x=times_pool_steady$X1, probs=0.025)
  ub <- quantile(x=times_pool_steady$X1, probs=0.975)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear+ symptom- 
  mu <- mean(times_pool_steady$X2)
  lb <- quantile(x=times_pool_steady$X2, probs=0.025)
  ub <- quantile(x=times_pool_steady$X2, probs=0.975)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear- symptom+
  mu <- mean(times_pool_steady$X3)
  lb <- quantile(x=times_pool_steady$X3, probs=0.025)
  ub <- quantile(x=times_pool_steady$X3, probs=0.975)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear+ symptom+
  mu <- mean(times_pool_steady$X4)
  lb <- quantile(x=times_pool_steady$X4, probs=0.025)
  ub <- quantile(x=times_pool_steady$X4, probs=0.975)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time smear+ (any)
  mu <- mean(times_pool_steady$smear)
  lb <- quantile(x=times_pool_steady$smear, probs=0.025)
  ub <- quantile(x=times_pool_steady$smear, probs=0.975)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  #mean and CI of total time symptom+ (any)
  if(symptom_available) {
    mu <- mean(times_pool_steady$symptom, na.rm=T)
    lb <- quantile(x=times_pool_steady$symptom, probs=0.025, na.rm=T)
    ub <- quantile(x=times_pool_steady$symptom, probs=0.975, na.rm=T)
    times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  } else {
    times_pop_steady <- c(times_pop_steady, "NA")
  }
  #mean and CI of total time w/ TB (any)
  mu <- mean(times_pool_steady$tb_any)
  lb <- quantile(x=times_pool_steady$tb_any, probs=0.025)
  ub <- quantile(x=times_pool_steady$tb_any, probs=0.975)
  times_pop_steady <- c(times_pop_steady, paste0(format(round(mu, 1), nsmall=1), " [", format(round(lb, 1), nsmall=1), "-", format(round(ub, 1), nsmall=1), "]"))
  
  #proportion ever reaching a state
  #exclude parameter sets if sims didn't finish for all 4 starting populations
  props_pool_steady <- sapply(countries, function(x)
    props_all[[x]] %>% mutate(start_pop=as.numeric(start_pop)) %>%
      group_by(id) %>% mutate(id_count=n())  %>% filter(id_count==4),
    simplify=F)
  #assign unique ID for each country row
  props_pool_steady <- sapply(countries, function(x)
    props_pool_steady[[x]] %>% ungroup() %>% group_by(start_pop) %>%
      mutate(id_unique=rep(1:n())),
    simplify=F)
  #sample 50k parameter sets for each country
  ids <- sapply(countries, function(x)
    sample(props_pool_steady[[x]] %>% filter(start_pop==1) %>% 
             pull("id_unique"), 50000, 
           replace=T,prob=props_pool_steady[[x]] %>% 
             filter(start_pop==1) %>% pull(weight)),
    simplify=F)
  row_nums <- sapply(countries, function(x)
    c(match(ids[[x]], props_pool_steady[[x]] %>% filter(start_pop==1) %>%
              pull("id_unique")),
      nrow(props_pool_steady[[x]])/4 + 
        match(ids[[x]], props_pool_steady[[x]] %>% filter(start_pop==1) %>%
                pull("id_unique")),
      2*nrow(props_pool_steady[[x]])/4 + 
        match(ids[[x]], props_pool_steady[[x]] %>% filter(start_pop==1) %>%
                pull("id_unique")),
      3*nrow(props_pool_steady[[x]])/4 + 
        match(ids[[x]], props_pool_steady[[x]] %>% filter(start_pop==1) %>%
                pull("id_unique"))),
    simplify=F)
  props_pool_steady <- sapply(countries, function(x)
    props_pool_steady[[x]][row_nums[[x]], ], simplify=F)
  #apply state population weights
  props_pool_steady <- sapply(countries, function(x)
    props_pool_steady[[x]] %>% 
      mutate(start_pop=as.numeric(start_pop)),
    simplify=F)
  props_pool_steady <- sapply(countries, function(x)
    props_pool_steady[[x]] %>%
      mutate(weight_steady=unname(unlist(steady_wts[[x]][start_pop]))),
    simplify=F) 
  #bind rows
  props_pool_steady <- bind_rows(props_pool_steady, .id="country")
  #calculate averages across all states (weighted by state dist) for each param set
  symptom_steady <- props_pool_steady %>%
    group_by(country, id) %>%
    summarise(symptom=weighted.mean(symptom, w=weight_steady))
  props_pool_steady <- props_pool_steady %>% group_by(country, id) %>% 
    select(-c(id_count, start_pop, id_unique, symptom)) %>%
    summarise_all(~weighted.mean(., w=weight_steady))
  props_pool_steady <- cbind(props_pool_steady, "symptom"=(symptom_steady %>% pull(symptom)))
  #mean and CI of total time smear- symptom-
  mu <- 100*mean(props_pool_steady$X1, weights=props_pool_steady$weight)
  lb <- 100*quantile(x=props_pool_steady$X1, probs=0.025, weight=props_pool_steady$weight)
  ub <- 100*quantile(x=props_pool_steady$X1, probs=0.975, weight=props_pool_steady$weight)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of total time smear+ symptom- 
  mu <- 100*mean(props_pool_steady$X2, weights=props_pool_steady$weight)
  lb <- 100*quantile(x=props_pool_steady$X2, probs=0.025, weight=props_pool_steady$weight)
  ub <- 100*quantile(x=props_pool_steady$X2, probs=0.975, weight=props_pool_steady$weight)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of total time smear- symptom+
  mu <- 100*mean(props_pool_steady$X3, weights=props_pool_steady$weight)
  lb <- 100*quantile(x=props_pool_steady$X3, probs=0.025, weight=props_pool_steady$weight)
  ub <- 100*quantile(x=props_pool_steady$X3, probs=0.975, weight=props_pool_steady$weight)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of total time smear+ symptom+
  mu <- 100*mean(props_pool_steady$X4, weights=props_pool_steady$weight)
  lb <- 100*quantile(x=props_pool_steady$X4, probs=0.025, weight=props_pool_steady$weight)
  ub <- 100*quantile(x=props_pool_steady$X4, probs=0.975, weight=props_pool_steady$weight)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of total time smear+ (any)
  mu <- 100*mean(props_pool_steady$smear, weights=props_pool_steady$weight)
  lb <- 100*quantile(x=props_pool_steady$smear, probs=0.025, weight=props_pool_steady$weight)
  ub <- 100*quantile(x=props_pool_steady$smear, probs=0.975, weight=props_pool_steady$weight)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of total time symptom+ (any)
  symptom_available <- T
  if(symptom_available) {
    mu <- 100*mean(props_pool_steady$symptom, weights=props_pool_steady$weight, na.rm=T)
    lb <- 100*quantile(x=props_pool_steady$symptom, probs=0.025, weight=props_pool_steady$weight, na.rm=T)
    ub <- 100*quantile(x=props_pool_steady$symptom, probs=0.975, weight=props_pool_steady$weight, na.rm=T)
    times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  } else {
    times_pop_steady <- c(times_pop_steady, "NA")
  }
  #mean and CI of proportion reaching spontaneous resolution
  mu <- 100*mean(props_pool_steady$X5, weights=props_pool_steady$weight)
  lb <- 100*quantile(x=props_pool_steady$X5, probs=0.025, weight=props_pool_steady$weight)
  ub <- 100*quantile(x=props_pool_steady$X5, probs=0.975, weight=props_pool_steady$weight)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching treatment
  mu <- 100*mean(props_pool_steady$X6, weights=props_pool_steady$weight)
  lb <- 100*quantile(x=props_pool_steady$X6, probs=0.025, weight=props_pool_steady$weight)
  ub <- 100*quantile(x=props_pool_steady$X6, probs=0.975, weight=props_pool_steady$weight)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  #mean and CI of proportion reaching death (TB only)
  mu <- 100*mean(props_pool_steady$X7, weights=props_pool_steady$weight)
  lb <- 100*quantile(x=props_pool_steady$X7, probs=0.025, weight=props_pool_steady$weight)
  ub <- 100*quantile(x=props_pool_steady$X7, probs=0.975, weight=props_pool_steady$weight)
  times_pop_steady <- c(times_pop_steady, paste0(round(mu), "% [", round(lb), "-", round(ub), "%]"))
  times_pool_out[[5]] <- unlist(times_pop_steady)
}
times_pool_out <- bind_cols(times_pool_out)
write.csv(times_pool_out, file=paste0(path_out, "times_out_pooled_", scenario_lab, ".csv"), row.names=F)


#########################################################################
#PART 2: RELATIVE INFECTIONS CALCULATIONS
#########################################################################
#calculate the relative number of infections produced by avg. person starting from each state
#only include param IDs that ran across all 4 sets of simulations (each start_pop)
rel_inf_out_all <- list()
rel_inf_samples_all <- list()
prop_trans_sum_all <- list()
rel_inf_smear_all <- list()
rel_smear_samples_all <- list()
rel_inf_smear_sa_all <- list()
#define rel infectiousness parameters
rr_m_mu <- 0.35 #relative infectiousness of smear-, compared to smear+
rr_s_mu <- 0.7 #relative infectiousness of subclinical, compared to symptomatic
rr_m <- rbeta(n=100000, shape1=9.7, shape2=17.6) #fit to 35% [20-55%]
rr_s <- rgamma(n=100000, shape=34, scale=0.021) #use gamma given skew: 70% [50-100%]
rr_s[rr_s>1] <- 1
inf_mu <- 1
inf_m_mu <- inf_mu/rr_m_mu
inf_s_mu <- inf_mu/rr_s_mu
inf_ms_mu <- inf_m_mu/rr_s_mu
inf <- 1
inf_m <- inf/rr_m
inf_s <- inf/rr_s 
inf_ms <- inf_m/rr_s
rel_inf_params <- data.frame("inf"=rep(1, 100000), "inf_m"=inf_m, "inf_s"=inf_s, "inf_ms"=inf_ms)
for(i in countries) {
  wts <- unname(unlist(steady_wts[[i]]))
  size <- sizes[[i]]
  times_all <- times_all_all[[i]]
  ids_1 <- as.character(unique((times_all %>% filter(start_pop==1))[["mean"]] %>% pull(id)))
  ids_2 <- as.character(unique((times_all %>% filter(start_pop==2))[["mean"]] %>% pull(id)))
  ids_3 <- as.character(unique((times_all %>% filter(start_pop==3))[["mean"]] %>% pull(id)))
  ids_4 <- as.character(unique((times_all %>% filter(start_pop==4))[["mean"]] %>% pull(id)))
  ids_sub <- Reduce(intersect, list(ids_1, ids_2, ids_3, ids_4))
  
  #used as reference
  times_sub_1 <- (times_all %>% filter(start_pop==1))[["mean"]]
  times_sub_1 <- times_sub_1 %>% filter(id %in% ids_sub)
  
  #sample population weights to use in calculations
  props <- data.frame(t(rmultinom(n=length(ids_sub), size=size, prob=wts))/size)
  names(props) <- c("prop1", "prop2", "prop3", "prop4")
  props <- props %>% mutate(id=ids_sub)
  props_lb <- colQuantiles(as.matrix(props[,1:4]), probs=0.025)
  props_ub <- colQuantiles(as.matrix(props[,1:4]), probs=0.975)
  
  rel_inf_out <- list()
  rel_inf_samples <- list()
  prop_trans_out <- list()
  for(j in 1:4) {
    print(j)
    times_pop <- c()
    times_all_state <- (times_all %>% filter(start_pop==j))[["mean"]]
    #need to match IDs to calculate relative infections for each param set
    times_sub_state <- times_all_state %>% filter(id %in% ids_sub)
    times_sub_state <- left_join(times_sub_state, props, by="id")
    #match relative infectiousness samples too - i.e. rel_inf is consistent among all ppl within a simulation, but varies across simulatinos
    times_sub_state <- cbind(times_sub_state, rel_inf_params[1:nrow(times_sub_state),])
    
    #calculate steady state % of transmission
    times_sub_state <- times_sub_state %>% 
      mutate(prop_trans=inf*X1*prop1+inf_m*X2*prop2+inf_s*X3*prop3+inf_ms*X4*prop4)
    prop_trans_out[[j]] <- times_sub_state$prop_trans
    
    #only calculate relative infectiousness for pops 2-4 (1 is reference)
    rel_inf_all <- (times_sub_state$inf*times_sub_state$X1 + 
                      times_sub_state$inf_m*times_sub_state$X2 + 
                      times_sub_state$inf_s*times_sub_state$X3 + 
                      times_sub_state$inf_ms*times_sub_state$X4)/(
                        times_sub_state$inf*times_sub_1$X1 + 
                          times_sub_state$inf_m*times_sub_1$X2 + 
                          times_sub_state$inf_s*times_sub_1$X3 + 
                          times_sub_state$inf_ms*times_sub_1$X4)
    if(j!=1) {
      rel_inf_mu <- wtd.mean(x=(times_sub_state$inf*times_sub_state$X1 + 
                                  times_sub_state$inf_m*times_sub_state$X2 + 
                                  times_sub_state$inf_s*times_sub_state$X3 + 
                                  times_sub_state$inf_ms*times_sub_state$X4)/
                               (times_sub_state$inf*times_sub_1$X1 + 
                                  times_sub_state$inf_m*times_sub_1$X2 + 
                                  times_sub_state$inf_s*times_sub_1$X3 + 
                                  times_sub_state$inf_ms*times_sub_1$X4),
                             weight=times_sub_state$weight)
      
      rel_inf_lb <- wtd.quantile(x=(times_sub_state$inf*times_sub_state$X1 + 
                                      times_sub_state$inf_m*times_sub_state$X2 + 
                                      times_sub_state$inf_s*times_sub_state$X3 + 
                                      times_sub_state$inf_ms*times_sub_state$X4)/
                                   (times_sub_state$inf*times_sub_1$X1 + 
                                      times_sub_state$inf_m*times_sub_1$X2 + 
                                      times_sub_state$inf_s*times_sub_1$X3 + 
                                      times_sub_state$inf_ms*times_sub_1$X4),
                                 probs=0.025,
                                 weight=times_sub_state$weight,
                                 normwt=T)[[1]]
      
      rel_inf_ub <- wtd.quantile(x=(times_sub_state$inf*times_sub_state$X1 + 
                                      times_sub_state$inf_m*times_sub_state$X2 + 
                                      times_sub_state$inf_s*times_sub_state$X3 + 
                                      times_sub_state$inf_ms*times_sub_state$X4)/
                                   (times_sub_state$inf*times_sub_1$X1 + 
                                      times_sub_state$inf_m*times_sub_1$X2 + 
                                      times_sub_state$inf_s*times_sub_1$X3 + 
                                      times_sub_state$inf_ms*times_sub_1$X4),
                                 probs=0.975,
                                 weight=times_sub_state$weight,
                                 normwt=T)[[1]]
      
      rel_inf <- c("mu"=rel_inf_mu, "lb"=rel_inf_lb, "ub"=rel_inf_ub)
      rel_inf_out[[j]] <- rel_inf
      rel_inf_samples[[j]] <- rel_inf_all
    }
  }
  rel_inf_out <- bind_rows(rel_inf_out)
  rel_inf_out <- data.frame(rel_inf_out, "name"=c("Smear+ Subclinical", 
                                                  "Smear- Symptomatic", 
                                                  "Smear+ Symptomatic"))
  rel_inf_out <- rel_inf_out %>% mutate(name=factor(name, levels=c("Smear+ Subclinical", 
                                                                   "Smear- Symptomatic", 
                                                                   "Smear+ Symptomatic")))
  prop_trans_out <- bind_cols(prop_trans_out)
  names(prop_trans_out) <- names(props)[1:4]
  prop_trans_out <- prop_trans_out %>% mutate(id=ids_sub)
  prop_trans_out <- left_join(prop_trans_out, times_sub_state %>% select(id, weight), by="id")
  #normalize to 100%
  prop_trans_out <- prop_trans_out %>% mutate(norm_sum=prop1+prop2+prop3+prop4,
                                              prop1=prop1/norm_sum,
                                              prop2=prop2/norm_sum,
                                              prop3=prop3/norm_sum,
                                              prop4=prop4/norm_sum
  )
  #summary states
  prop_trans_sum <- data.frame("prop_trans_mu"=c(wtd.mean(prop_trans_out$prop1, weight=prop_trans_out$weight),
                                                 wtd.mean(prop_trans_out$prop2, weight=prop_trans_out$weight),
                                                 wtd.mean(prop_trans_out$prop3, weight=prop_trans_out$weight),
                                                 wtd.mean(prop_trans_out$prop4, weight=prop_trans_out$weight)),
                               "prop_trans_lb"=c(wtd.quantile(x=prop_trans_out$prop1, probs=0.025,
                                                              weight=prop_trans_out$weight,
                                                              normwt=T)[[1]],
                                                 wtd.quantile(x=prop_trans_out$prop2, probs=0.025,
                                                              weight=prop_trans_out$weight,
                                                              normwt=T)[[1]],
                                                 wtd.quantile(x=prop_trans_out$prop3, probs=0.025,
                                                              weight=prop_trans_out$weight,
                                                              normwt=T)[[1]],
                                                 wtd.quantile(x=prop_trans_out$prop4, probs=0.025,
                                                              weight=prop_trans_out$weight,
                                                              normwt=T)[[1]]),
                               "prop_trans_ub"=c(wtd.quantile(x=prop_trans_out$prop1, probs=0.975,
                                                              weight=prop_trans_out$weight,
                                                              normwt=T)[[1]],
                                                 wtd.quantile(x=prop_trans_out$prop2, probs=0.975,
                                                              weight=prop_trans_out$weight,
                                                              normwt=T)[[1]],
                                                 wtd.quantile(x=prop_trans_out$prop3, probs=0.975,
                                                              weight=prop_trans_out$weight,
                                                              normwt=T)[[1]],
                                                 wtd.quantile(x=prop_trans_out$prop4, probs=0.975,
                                                              weight=prop_trans_out$weight,
                                                              normwt=T)[[1]])
  )
  prop_trans_sum <- prop_trans_sum %>% mutate(pop_prop=wts,
                                              pop_prop_lb=props_lb,
                                              pop_prop_ub=props_ub)
  prop_trans_sum <- prop_trans_sum %>% 
    mutate(name=state_names)
  prop_trans_sum <- prop_trans_sum %>%
    mutate(lab=paste0(round(prop_trans_mu*100, 0), "% [", 
                      round(prop_trans_lb*100, 0), "-",
                      round(prop_trans_ub*100, 0), "%]"),
           lab2=paste0(round(prop_trans_mu*100, 0), "%"),
           pop_lab=paste0(round(pop_prop*100, 0), "% [", 
                          round(pop_prop_lb*100, 0), "-",
                          round(pop_prop_ub*100, 0), "%]"),
           pop_lab2=paste0(round(pop_prop*100, 0), "%")
    )
  prop_trans_sum <- prop_trans_sum %>%
    mutate(name=factor(name, levels=prop_trans_sum$name),
           lab=factor(lab, levels=prop_trans_sum$lab),
           pop_lab=factor(pop_lab, levels=unique(prop_trans_sum$pop_lab)))
  
  rel_inf_samples <- bind_cols(rel_inf_samples)
  names(rel_inf_samples) <- c("Smear+ Subclinical", 
                              "Smear- Symptomatic", 
                              "Smear+ Symptomatic")
  
  #calculate smear+ subclinical relative to smear+ symptomatic (to see if statistically significant)
  times_all_2 <- (times_all %>% filter(start_pop==2))[["mean"]]
  times_all_4 <- (times_all %>% filter(start_pop==4))[["mean"]]
  times_sub_2 <- times_all_2 %>% filter(id %in% ids_sub)
  times_sub_4 <- times_all_4 %>% filter(id %in% ids_sub)
  rel_smear_samples <- (times_sub_state$inf*times_sub_2$X1 + 
                          times_sub_state$inf_m*times_sub_2$X2 + 
                          times_sub_state$inf_s*times_sub_2$X3 + 
                          times_sub_state$inf_ms*times_sub_2$X4)/
    (times_sub_state$inf*times_sub_4$X1 + 
       times_sub_state$inf_m*times_sub_4$X2 + 
       times_sub_state$inf_s*times_sub_4$X3 + 
       times_sub_state$inf_ms*times_sub_4$X4)
  
  rel_inf_mu <- wtd.mean(x=(times_sub_state$inf*times_sub_2$X1 + 
                              times_sub_state$inf_m*times_sub_2$X2 + 
                              times_sub_state$inf_s*times_sub_2$X3 + 
                              times_sub_state$inf_ms*times_sub_2$X4)/
                           (times_sub_state$inf*times_sub_4$X1 + 
                              times_sub_state$inf_m*times_sub_4$X2 + 
                              times_sub_state$inf_s*times_sub_4$X3 + 
                              times_sub_state$inf_ms*times_sub_4$X4),
                         weight=times_sub_2$weight)
  
  rel_inf_lb <- wtd.quantile(x=(times_sub_state$inf*times_sub_2$X1 + 
                                  times_sub_state$inf_m*times_sub_2$X2 + 
                                  times_sub_state$inf_s*times_sub_2$X3 + 
                                  times_sub_state$inf_ms*times_sub_2$X4)/
                               (times_sub_state$inf*times_sub_4$X1 + 
                                  times_sub_state$inf_m*times_sub_4$X2 + 
                                  times_sub_state$inf_s*times_sub_4$X3 + 
                                  times_sub_state$inf_ms*times_sub_4$X4),
                             probs=0.025,
                             weight=times_sub_2$weight,
                             normwt=T)[[1]]
  
  rel_inf_ub <- wtd.quantile(x=(times_sub_state$inf*times_sub_2$X1 + 
                                  times_sub_state$inf_m*times_sub_2$X2 + 
                                  times_sub_state$inf_s*times_sub_2$X3 + 
                                  times_sub_state$inf_ms*times_sub_2$X4)/
                               (times_sub_state$inf*times_sub_4$X1 + 
                                  times_sub_state$inf_m*times_sub_4$X2 + 
                                  times_sub_state$inf_s*times_sub_4$X3 + 
                                  times_sub_state$inf_ms*times_sub_4$X4),
                             probs=0.975,
                             weight=times_sub_2$weight,
                             normwt=T)[[1]]
  
  rel_inf_smear <- data.frame("mu"=rel_inf_mu, "lb"=rel_inf_lb, "ub"=rel_inf_ub,
                              "name"="Smear+ Subclinical")
  
  #sensitivity analysis on (smear+) subclinical vs. (smear+) symptomatic
  rel_inf_smear_sa <- list()
  for(rr_s in (1:100)/100) {
    inf_s_mu <- inf_mu/rr_s 
    inf_ms_mu <- inf_m_mu/rr_s
    ids_2 <- unique(times_all_2$id)
    ids_4 <- unique(times_all_4$id)
    times_sub_2 <- times_all_2 %>% filter(id %in% ids_4)
    times_sub_4 <- times_all_4 %>% filter(id %in% ids_2)
    rel_inf_mu <- wtd.mean(x=(inf_mu*times_sub_2$X1 + 
                                inf_m_mu*times_sub_2$X2 + 
                                inf_s_mu*times_sub_2$X3 + 
                                inf_ms_mu*times_sub_2$X4)/
                             (inf_mu*times_sub_4$X1 + 
                                inf_m_mu*times_sub_4$X2 + 
                                inf_s_mu*times_sub_4$X3 + 
                                inf_ms_mu*times_sub_4$X4),
                           weight=times_sub_2$weight)
    
    rel_inf_lb <- wtd.quantile(x=(inf_mu*times_sub_2$X1 + 
                                    inf_m_mu*times_sub_2$X2 + 
                                    inf_s_mu*times_sub_2$X3 + 
                                    inf_ms_mu*times_sub_2$X4)/
                                 (inf_mu*times_sub_4$X1 + 
                                    inf_m_mu*times_sub_4$X2 + 
                                    inf_s_mu*times_sub_4$X3 + 
                                    inf_ms_mu*times_sub_4$X4),
                               probs=0.025,
                               weight=times_sub_2$weight,
                               normwt=T)[[1]]
    
    rel_inf_ub <- wtd.quantile(x=(inf_mu*times_sub_2$X1 + 
                                    inf_m_mu*times_sub_2$X2 + 
                                    inf_s_mu*times_sub_2$X3 + 
                                    inf_ms_mu*times_sub_2$X4)/
                                 (inf_mu*times_sub_4$X1 + 
                                    inf_m_mu*times_sub_4$X2 + 
                                    inf_s_mu*times_sub_4$X3 + 
                                    inf_ms_mu*times_sub_4$X4),
                               probs=0.975,
                               weight=times_sub_2$weight,
                               normwt=T)[[1]]
    
    rel_inf <- c("mu"=rel_inf_mu, "lb"=rel_inf_lb, "ub"=rel_inf_ub,
                 "rr_s"=rr_s)
    rel_inf_smear_sa[[as.character(rr_s)]] <- rel_inf
  }
  rel_inf_smear_sa <- bind_rows(rel_inf_smear_sa)
  rel_inf_smear_sa <- rel_inf_smear_sa %>% mutate(name="Smear+ Subclinical")
  
  rel_inf_out_all[[i]] <- rel_inf_out
  rel_inf_samples_all[[i]] <- rel_inf_samples
  rel_inf_smear_all[[i]] <- rel_inf_smear
  rel_smear_samples_all[[i]] <- rel_smear_samples
  prop_trans_sum_all[[i]] <- prop_trans_sum
  rel_inf_smear_sa_all[[i]] <- rel_inf_smear_sa
}
rel_inf_out_all <- bind_rows(rel_inf_out_all, .id="country")
rel_inf_samples_all <- bind_rows(rel_inf_samples_all, .id="country")
rel_inf_samples_all <- pivot_longer(rel_inf_samples_all, cols=2:4, 
                                    names_to="state", values_to="rel_inf")
rel_smear_samples_all <- stack(rel_smear_samples_all)
names(rel_smear_samples_all) <- c("rel_inf", "country")
rel_smear_samples_all <- rel_smear_samples_all %>% mutate(state="Smear+ Subclinical",
                                                          country=as.character(country))
rel_inf_smear_all <- bind_rows(rel_inf_smear_all, .id="country")
rel_inf_smear_all <- 
  rel_inf_smear_all %>% mutate(inf_m=0.7)
rel_inf_smear_sa_all <- bind_rows(rel_inf_smear_sa_all, .id="country")
prop_trans_sum_all <- bind_rows(prop_trans_sum_all, .id="country")

#min rr_s by country where rel_inf > 1
rel_inf_smear_sa_all %>% filter(mu>=1) %>% group_by(country) %>%
  summarise(rr_s=min(rr_s)) %>% select(country, rr_s)

#export data for appendix Table 2
#part 1: durations
rel_inf_out <- rel_inf_out_all %>% mutate(rel_to="Smear- Subclinical")
rel_inf_pool <- rel_inf_samples_all %>% group_by(state) %>% 
  summarise(mu=mean(rel_inf), lb=quantile(rel_inf, 0.025), ub=quantile(rel_inf, 0.975)) %>%
  mutate(country="Pooled", rel_to="Smear- Subclinical") %>% 
  rename("name"="state")
rel_inf_out <- bind_rows(rel_inf_out, rel_inf_pool)

rel_smear_out <- rel_inf_smear_all %>% mutate(rel_to="Smear+ Symptomatic") %>%
  select(-inf_m)
rel_smear_pool <- rel_smear_samples_all %>% group_by(state) %>% 
  summarise(mu=mean(rel_inf), lb=quantile(rel_inf, 0.025), ub=quantile(rel_inf, 0.975)) %>%
  mutate(country="Pooled", rel_to="Smear+ Symptomatic") %>%
  rename("name"="state")
rel_inf_out <- bind_rows(rel_inf_out, rel_smear_out)
rel_inf_out <- bind_rows(rel_inf_out, rel_smear_pool)

rel_inf_out <- rel_inf_out %>% 
  mutate(est=paste0(format(round(mu, 1), nsmall=1), " [", 
                    format(round(lb, 1), nsmall=1), "-", 
                    format(round(ub, 1), nsmall=1), "]")) %>% 
  select(-c(mu, lb, ub))
rel_inf_out <- pivot_wider(rel_inf_out, id_cols=c("name", "rel_to"),
                           names_from=country, values_from=est)
rel_inf_out <- rel_inf_out[, c("name", sort(countries), "Pooled", "rel_to")]
write.csv(rel_inf_out, file=paste0(path_out, "rel_inf_pp_table_", scenario_lab, ".csv"), row.names=F)

#part 2: proportions
props_out <- prop_trans_sum_all %>% select(country, name, lab, pop_lab)
props_out <- props_out %>% pivot_longer(props_out, cols=contains("lab"),
                                        names_to="estimate", values_to="value") #warning here is fine
props_out <- props_out %>% 
  mutate(estimate=if_else(estimate=="lab", "trans_prop", "pop_prop"))
props_out <- pivot_wider(props_out, id_cols=c("name", "estimate"),
                         names_from=country, values_from="value")
props_out <- props_out %>% arrange(estimate, name)
props_out <- props_out[, c("estimate", "name", sort(countries))]
write.csv(props_out, file=paste0(path_out, "/rel_inf_pop_table_", scenario_lab, ".csv"), row.names=F)


  
#graphs (main text figure 4)
if(scenario_lab=="base")  {   #don't need this figure for sensitivity analyses - just use the tables
  #labels for graphs
  rel_inf_out_long <- pivot_longer(rel_inf_out, cols=sort(countries),
                                   names_to="country", 
                                   values_to="rel_inf")
  rel_inf_out_long <- rel_inf_out_long %>%
    mutate(x_pos=as.numeric(str_split(rel_inf, " ", simplify=T)[,1]) +
             case_when(name=="Smear- Symptomatic"~1.2, 
                       rel_to=="Smear+ Symptomatic"~0.6, 
                       TRUE~0), #place to right of tall peaks
           y_pos=as.numeric(factor(country, levels=sort(countries, decreasing=T)))+
             case_when(name=="Smear- Symptomatic"~0.8,
                       name=="Smear+ Symptomatic"~0.4,
                       rel_to=="Smear+ Symptomatic"~0.75,
                       TRUE~0.6)) #make position higher for tall peaks
  fig1 <- ggplot(rel_inf_samples_all, aes(x=rel_inf, y=country, fill=state)) +
    geom_density_ridges(alpha=0, scale=1.5,
                        quantile_lines=T, quantile_fun=function(x,...)mean(x), bandwidth=0.05) +
    stat_density_ridges(scale=1.5, quantile_lines=T, quantiles=c(0.025, 0.975), alpha=0.5, bandwidth=0.05) +
    geom_vline(aes(xintercept=1), linetype="dashed") +
    geom_text(data=rel_inf_out_long %>% filter(rel_to=="Smear- Subclinical"), 
              aes(x=x_pos, y=y_pos, label=rel_inf, fill=NULL), size=3) +
    scale_y_discrete(limits = rev,
                     expand=expansion(mult=c(0, 0.05))) +
    scale_x_continuous(limits=c(0, 11), breaks=1:10,
                       expand=c(0,0)) +
    scale_fill_manual(values=colors_s) +
    ggtitle("A. Relative cumulative secondary infections") +
    labs(x="vs. smear-negative subclinical", y="", fill="") +
    theme_bw()  + theme(plot.title=element_text(size=11, face="bold"),
                        axis.title.y=element_blank(),
                        axis.title.x=element_text(size=10),
                        plot.margin = unit(c(0,0.2,0,0), "cm"),
                        legend.position="none")
  fig2 <- ggplot(rel_smear_samples_all, aes(x=rel_inf, y=country, fill=state)) +
    geom_density_ridges(alpha=0, scale=1.5,
                        quantile_lines=T, quantile_fun=function(x,...)mean(x)) +
    stat_density_ridges(scale=1.5, quantile_lines=T, quantiles=c(0.025, 0.975), alpha=0.5) +
    geom_vline(aes(xintercept=1), linetype="dashed") +
    geom_text(data=rel_inf_out_long %>% filter(rel_to=="Smear+ Symptomatic"), 
              aes(x=x_pos, y=y_pos, label=rel_inf, fill=NULL), size=3) +
    scale_y_discrete(limits = rev,
                     expand=expansion(mult=c(0, 0.05))) +
    scale_x_continuous(limits=c(0, 3), breaks=1:3,
                       expand=c(0,0)) +
    scale_fill_manual(values=colors_s) +
    ggtitle("") +
    labs(x="vs. smear-pos. symptomatic", y="", fill="") +
    theme_bw()  + theme(plot.title=element_text(size=11, face="bold"),
                        axis.title.y=element_blank(),
                        axis.title.x=element_text(size=10),
                        plot.margin = unit(c(0,0.1,0,0.2), "cm"),
                        legend.position="none",
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank())
  plot_top <- wrap_plots(fig1, fig2, ncol=2, nrow=1, 
                         widths=c(0.65, 0.35))
  
  #MIDDLE PANEL: sensitivity analysis on smear
  fig3 <- ggplot(rel_inf_smear_all, aes(x=inf_m, color=country)) +
    #geom_point(aes(y=mu), position=position_jitter(seed=123, width=0.02)) +
    #geom_errorbar(aes(y=mu, ymin=lb, ymax=ub), width=0.02,
    #               position=position_jitter(seed=123, width=0.02)) +
    geom_line(data=rel_inf_smear_sa_all, 
              aes(x=rr_s, y=mu, color=country), alpha=0.25, size=1) +
    geom_line(data=rel_inf_smear_sa_all, 
              aes(x=rr_s, y=mu, color=country, linetype=country), size=1) +
    geom_hline(aes(yintercept=1)) +
    geom_vline(aes(xintercept=0.7), linetype="dashed") +
    geom_text(x=0.71, y=1.7, 
              label="Base Estimate\n(70%)", 
              color="black", size=3, fontface=1, hjust=0) +
    scale_y_continuous(expand=expansion(mult=c(.05, .05))) +
    scale_x_continuous(expand=expansion(mult=c(0, 0)),
                       breaks=c(0.25, 0.5, 0.75, 1),
                       labels=scales::percent_format(accuracy=1)) +
    scale_color_manual(values=colors_c[sort(unique(rel_inf_smear_all$country))]) +
    scale_linetype_manual(values=c("dashed", "dotted", "longdash", "solid", "twodash")) +
    labs(x="Relative infectiousness of subclinical (vs. symptomatic) TB", 
         y="Relative cumulative\nsecondary infections",
         fill="", color="", linetype="") +
    theme_bw() + theme(panel.grid=element_blank()) + 
    ggtitle("B. Per-person contribution to transmission over 5 years,\nsmear-positive subclinical vs. smear-positive symptomatic") +
    theme(plot.title=element_text(size=11, face="bold"),
          axis.title=element_text(size=10),
          plot.margin = unit(c(0.5,0,0,0), "cm"),
          legend.key.width=unit(1, "cm"))
  #BOTTOM PANEL: stacked bar charts
  prop_pop <- prop_trans_sum_all %>% 
    select(country, name, starts_with("pop_lab"), starts_with("pop_prop"))
  prop_trans <- prop_trans_sum_all %>% 
    select(country, name, starts_with("lab"), starts_with("prop_trans")) 
  fig4 <- ggplot(prop_pop) + 
    geom_col(aes(x=country, y=pop_prop, fill=fct_rev(name)), 
             width=0.95, position="stack") +
    geom_text(aes(x=country, y=pop_prop, 
                  label=fct_rev(pop_lab2)), 
              position=position_stack(vjust=0.5), size=3, fontface=1, color="white") +
    scale_y_continuous(expand=expansion(mult=c(0, 0)),
                       labels=scales::percent_format(accuracy=1)) +
    scale_fill_manual(values=colors_s) + scale_color_manual(values=c("white", "black", "black", "black")) +
    labs(x="", y="", fill="Initial TB state", color="", group="") +
    ggtitle("C. Population contribution to\nTB prevalence") +
    guides(color="none") +
    theme_bw() + theme(panel.grid=element_blank()) +
    theme(plot.title=element_text(size=11, face="bold"),
          legend.margin=margin(t=0, unit="cm"),
          axis.title=element_blank(),
          plot.margin = unit(c(0.5,0,0,0), "cm"))
  
  fig5 <- ggplot(prop_trans) + 
    geom_col(aes(x=country, y=prop_trans_mu, fill=fct_rev(name)),
             width=0.95, position="stack") +
    geom_text(aes(x=country, y=prop_trans_mu, 
                  label=fct_rev(lab2)), 
              position=position_stack(vjust=0.5), size=3, fontface=1, color="white") +
    scale_y_continuous(expand=expansion(mult=c(0, 0)),
                       labels=scales::percent_format(accuracy=1)) +
    scale_fill_manual(values=colors_s) + scale_color_manual(values=c("white", "black", "black", "black")) +
    labs(x="", y="", fill="Initial TB state", color="", group="") +
    ggtitle("D. Population contribution to\ncumulative 5-year transmission") +
    guides(color="none") +
    theme_bw() + theme(panel.grid=element_blank()) +
    theme(plot.title=element_text(size=11, face="bold"),
          legend.margin=margin(t=0, unit="cm"),
          axis.title=element_blank(),
          plot.margin = unit(c(0.5,0.2,0,0), "cm"),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
  plot_bottom <- plot_grid(fig4 + theme(legend.position="none"),
                           fig5 + theme(legend.position="none"),
                           nrow=1, ncol=2, align="hv")
  #combine all panels
  plot <- plot_grid(get_legend(fig4+theme(legend.position="bottom")),
                    plot_top, fig3, plot_bottom, 
                    nrow=4, ncol=1, align="hv",
                    rel_heights=c(0.1, 0.7, 0.5, 0.5))
  ggsave(plot, filename=paste0(path_out, "/rel_inf.jpg"), 
         dpi=500, height=11, width=7.5)
}


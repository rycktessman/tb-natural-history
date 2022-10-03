#generate output for heatmap and for appendix tables for sens analyses
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
setwd("C:/Users/Tess/OneDrive - Johns Hopkins/TB/Natural History Modeling")
source("code/model_v3.R")
source("code/calib_functions2.R")

#% by TB state for steady state calcs
country <- "Nepal"
steady_wts <- list()
load(paste0("data/params_targets_", tolower(country), ".Rda"))
steady_wts <- c(1-(targets_all[["prop_m_all"]] + targets_all[["prop_s_all"]] - targets_all[["prop_ms"]]),
         targets_all[["prop_m_all"]] - targets_all[["prop_ms"]],
         targets_all[["prop_s_all"]] - targets_all[["prop_ms"]],
         targets_all[["prop_ms"]])

countries <- c("Philippines", "Vietnam", "Nepal", "Cambodia", "Bangladesh")
sizes <- c(356, 219, 225, 314, 278)
names(sizes) <- countries
size <- sizes[[country]]

state_names <- c("Smear- Subclinical",
                 "Smear+ Subclinical", 
                 "Smear- Symptomatic", 
                 "Smear+ Symptomatic")

#load files for each sensitivity analysis and each starting population
#run analyses for each sensitivity analysis

#repeat for each sensitivity analysis - JUST CHANGE THE PATH
path <- paste0("output/IMIS Nov2021 v4 ", tolower(country), " smearhist")
out_post <- read.csv(paste0(path, "/out_IMIS_combined.csv"), stringsAsFactors=F)
props_c <- list()
times_all_c <- list()
for(start_pop in 1:4) {
  props <- read.csv(paste0(path, "/props_comb", start_pop, ".csv"), stringsAsFactors=F)
  load(paste0(path, "/times_all_comb", start_pop, ".Rda"))
  props_c[[start_pop]] <- props
  times_all_c[[start_pop]] <- times_all_comb
}
props_c <- bind_rows(props_c, .id="start_pop")
times_all_c <- bind_rows(times_all_c, .id="start_pop")

#PART 1A: INDIVIDUAL TRAJECTORIES SUMMARY STATS
#USED TO FORM HEATMAP/TABLE (heatmap is in separate script)

#each initial state separately
times_out_all <- list()
symptom_available <- T
times_out <- list()
for(j in 1:4) {
  print(j)
  times_pop <- c()
  times_all_comb <- times_all_c %>% filter(start_pop==j)
  props <- props_c %>% filter(start_pop==j)
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
write.csv(times_out, file=paste0(path, "/times_out.csv"), 
          row.names=F)

#steady state distribution
if(TRUE) {
  times_pop_steady <- list()
  
  times_all_steady <- bind_rows(lapply(1:4, function(x) times_all_c %>% filter(start_pop==x) %>% pull(mean)),
                                .id="start_pop")
  times_all_steady <- times_all_steady %>% mutate(start_pop=as.numeric(start_pop))
  #exclude parameter sets if sims didn't finish for all 4 starting populations
  times_all_steady <- times_all_steady %>% group_by(id) %>% mutate(id_count=n())
  times_all_steady <- times_all_steady %>% filter(id_count==4)
  #apply state population weights
  times_all_steady <- times_all_steady %>% mutate(weight_steady=steady_wts[start_pop])
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
  props_steady <- props_c
  props_steady <- props_steady %>% mutate(start_pop=as.numeric(start_pop))
  #exclude parameter sets if sims didn't finish for all 4 starting populations
  props_steady <- props_steady %>% group_by(id) %>% mutate(id_count=n())
  props_steady <- props_steady %>% filter(id_count==4)
  #apply state population weights
  props_steady <- props_steady %>% mutate(weight_steady=steady_wts[start_pop])
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
}
times_pop_steady <- data.frame(unlist(times_pop_steady))
write.csv(times_pop_steady, 
          file=paste0(path, "/times_out_steady.csv"), row.names=F)

#PART 1B: RELATIVE INFECTIONS CALCULATIONS
#calculate the relative number of infections produced by avg. person starting from each state
#only include param IDs that ran across all 4 sets of simulations (each start_pop)
if(TRUE) {
  #define rel infectiousness parameters
  rr_m <- 0.35 #relative infectiousness of smear-, compared to smear+
  rr_s <- 0.7 #relative infectiousness of subclinical, compared to symptomatic
  inf <- 1
  inf_m <- inf/rr_m
  inf_s <- inf/rr_s 
  inf_ms <- inf_m/rr_s
  wts <- steady_wts
  times_all <- times_all_c
  ids_1 <- as.character(unique((times_all %>% filter(start_pop==1))[["mean"]] %>% pull(id)))
  ids_2 <- as.character(unique((times_all %>% filter(start_pop==2))[["mean"]] %>% pull(id)))
  ids_3 <- as.character(unique((times_all %>% filter(start_pop==3))[["mean"]] %>% pull(id)))
  ids_4 <- as.character(unique((times_all %>% filter(start_pop==4))[["mean"]] %>% pull(id)))
  ids_sub <- Reduce(intersect, list(ids_1, ids_2, ids_3, ids_4))
  
  #used as reference
  times_sub_1 <- (times_all %>% filter(start_pop==1))[["mean"]]
  times_sub_1 <- times_sub_1 %>% filter(id %in% ids_sub)
  
  #sample population weights to use in calculations
  props <- data.frame(t(rmultinom(n=length(ids_sub), size=size, prob=wts))/size) #356=cases in prev survey
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
    
    #calculate steady state % of transmission
    times_sub_state <- times_sub_state %>% 
      mutate(prop_trans=inf*X1*prop1+inf_m*X2*prop2+inf_s*X3*prop3+inf_ms*X4*prop4)
    prop_trans_out[[j]] <- times_sub_state$prop_trans
    
    #only calculate relative infectiousness for pops 2-4 (1 is reference)
    rel_inf_all <- (inf*times_sub_state$X1 + 
                  inf_m*times_sub_state$X2 + 
                  inf_s*times_sub_state$X3 + 
                  inf_ms*times_sub_state$X4)/(inf*times_sub_1$X1 + 
                                                inf_m*times_sub_1$X2 + 
                                                inf_s*times_sub_1$X3 + 
                                                inf_ms*times_sub_1$X4)
    if(j!=1) {
      rel_inf_mu <- wtd.mean(x=(inf*times_sub_state$X1 + 
                                  inf_m*times_sub_state$X2 + 
                                  inf_s*times_sub_state$X3 + 
                                  inf_ms*times_sub_state$X4)/
                               (inf*times_sub_1$X1 + 
                                  inf_m*times_sub_1$X2 + 
                                  inf_s*times_sub_1$X3 + 
                                  inf_ms*times_sub_1$X4),
                             weight=times_sub_state$weight)
      
      rel_inf_lb <- wtd.quantile(x=(inf*times_sub_state$X1 + 
                                      inf_m*times_sub_state$X2 + 
                                      inf_s*times_sub_state$X3 + 
                                      inf_ms*times_sub_state$X4)/
                                   (inf*times_sub_1$X1 + 
                                      inf_m*times_sub_1$X2 + 
                                      inf_s*times_sub_1$X3 + 
                                      inf_ms*times_sub_1$X4),
                                 probs=0.025,
                                 weight=times_sub_state$weight,
                                 normwt=T)[[1]]
      
      rel_inf_ub <- wtd.quantile(x=(inf*times_sub_state$X1 + 
                                      inf_m*times_sub_state$X2 + 
                                      inf_s*times_sub_state$X3 + 
                                      inf_ms*times_sub_state$X4)/
                                   (inf*times_sub_1$X1 + 
                                      inf_m*times_sub_1$X2 + 
                                      inf_s*times_sub_1$X3 + 
                                      inf_ms*times_sub_1$X4),
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
  rel_smear_samples <- (inf*times_sub_2$X1 + 
                          inf_m*times_sub_2$X2 + 
                          inf_s*times_sub_2$X3 + 
                          inf_ms*times_sub_2$X4)/
    (inf*times_sub_4$X1 + 
       inf_m*times_sub_4$X2 + 
       inf_s*times_sub_4$X3 + 
       inf_ms*times_sub_4$X4)
  
  rel_inf_mu <- wtd.mean(x=(inf*times_sub_2$X1 + 
                              inf_m*times_sub_2$X2 + 
                              inf_s*times_sub_2$X3 + 
                              inf_ms*times_sub_2$X4)/
                           (inf*times_sub_4$X1 + 
                              inf_m*times_sub_4$X2 + 
                              inf_s*times_sub_4$X3 + 
                              inf_ms*times_sub_4$X4),
                         weight=times_sub_2$weight)
  
  rel_inf_lb <- wtd.quantile(x=(inf*times_sub_2$X1 + 
                                  inf_m*times_sub_2$X2 + 
                                  inf_s*times_sub_2$X3 + 
                                  inf_ms*times_sub_2$X4)/
                               (inf*times_sub_4$X1 + 
                                  inf_m*times_sub_4$X2 + 
                                  inf_s*times_sub_4$X3 + 
                                  inf_ms*times_sub_4$X4),
                             probs=0.025,
                             weight=times_sub_2$weight,
                             normwt=T)[[1]]
  
  rel_inf_ub <- wtd.quantile(x=(inf*times_sub_2$X1 + 
                                  inf_m*times_sub_2$X2 + 
                                  inf_s*times_sub_2$X3 + 
                                  inf_ms*times_sub_2$X4)/
                               (inf*times_sub_4$X1 + 
                                  inf_m*times_sub_4$X2 + 
                                  inf_s*times_sub_4$X3 + 
                                  inf_ms*times_sub_4$X4),
                             probs=0.975,
                             weight=times_sub_2$weight,
                             normwt=T)[[1]]
  
  rel_inf_smear <- data.frame("mu"=rel_inf_mu, "lb"=rel_inf_lb, "ub"=rel_inf_ub,
                              "name"="Smear+ Subclinical")
  
  #sensitivity analysis on (smear+) subclinical vs. (smear+) symptomatic
  rel_inf_smear_sa <- list()
  for(rr_s in (1:100)/100) {
    inf_m <- inf/rr_m
    inf_s <- inf/rr_s 
    inf_ms <- inf_m/rr_s
    ids_2 <- unique(times_all_2$id)
    ids_4 <- unique(times_all_4$id)
    times_sub_2 <- times_all_2 %>% filter(id %in% ids_4)
    times_sub_4 <- times_all_4 %>% filter(id %in% ids_2)
    rel_inf_mu <- wtd.mean(x=(inf*times_sub_2$X1 + 
                                inf_m*times_sub_2$X2 + 
                                inf_s*times_sub_2$X3 + 
                                inf_ms*times_sub_2$X4)/
                             (inf*times_sub_4$X1 + 
                                inf_m*times_sub_4$X2 + 
                                inf_s*times_sub_4$X3 + 
                                inf_ms*times_sub_4$X4),
                           weight=times_sub_2$weight)
    
    rel_inf_lb <- wtd.quantile(x=(inf*times_sub_2$X1 + 
                                    inf_m*times_sub_2$X2 + 
                                    inf_s*times_sub_2$X3 + 
                                    inf_ms*times_sub_2$X4)/
                                 (inf*times_sub_4$X1 + 
                                    inf_m*times_sub_4$X2 + 
                                    inf_s*times_sub_4$X3 + 
                                    inf_ms*times_sub_4$X4),
                               probs=0.025,
                               weight=times_sub_2$weight,
                               normwt=T)[[1]]
    
    rel_inf_ub <- wtd.quantile(x=(inf*times_sub_2$X1 + 
                                    inf_m*times_sub_2$X2 + 
                                    inf_s*times_sub_2$X3 + 
                                    inf_ms*times_sub_2$X4)/
                                 (inf*times_sub_4$X1 + 
                                    inf_m*times_sub_4$X2 + 
                                    inf_s*times_sub_4$X3 + 
                                    inf_ms*times_sub_4$X4),
                               probs=0.975,
                               weight=times_sub_2$weight,
                               normwt=T)[[1]]
    
    rel_inf <- c("mu"=rel_inf_mu, "lb"=rel_inf_lb, "ub"=rel_inf_ub,
                 "rr_s"=rr_s)
    rel_inf_smear_sa[[as.character(rr_s)]] <- rel_inf
  }
  rel_inf_smear_sa <- bind_rows(rel_inf_smear_sa)
  rel_inf_smear_sa <- rel_inf_smear_sa %>% mutate(name="Smear+ Subclinical")
  
}

#export data for appendix Table 2
#part 1: durations
rel_inf_out <- rel_inf_out %>% mutate(rel_to="Smear- Subclinical")
rel_inf_smear <- rel_inf_smear %>% mutate(rel_to="Smear+ Symptomatic")
rel_inf_out <- bind_rows(rel_inf_out, rel_inf_smear)

rel_inf_out <- rel_inf_out %>% 
  mutate(est=paste0(format(round(mu, 1), nsmall=1), " [", 
                    format(round(lb, 1), nsmall=1), "-", 
                    format(round(ub, 1), nsmall=1), "]")) %>% 
  select(-c(mu, lb, ub))
write.csv(rel_inf_out, file=paste0(path, "/rel_inf_pp_table.csv"), row.names=F)

#part 2: proportions
props_out <- prop_trans_sum %>% select(name, lab, pop_lab)
write.csv(props_out, file=paste0(path, "/rel_inf_pop_table.csv"), row.names=F)


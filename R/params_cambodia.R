setwd("C:/Users/Tess/OneDrive - Johns Hopkins/TB/Natural History Modeling")
load("data/params_targets.Rda")
library(MASS)
library(tidyverse)
library(dampack)
library(readxl)

targets_prev <- read_xlsx("data/country_targets.xlsx", range="A1:E1000")[,c(1,5)]
names(targets_prev) <- c("param", "value")

#define sorted rank function used to induce correlations
# Function inputs:
# 1. X: matrix of uncorrelated parameter draws (n samples by 4 parameters, in this case) to induce correlation
# 2. normals: matrix of multinormal samples (same dimensions) used to induce correlation on X
sorted_rank <- function(X, normals) {
  col <- ncol(normals)
  row <- nrow(normals)
  Xsorted <- matrix(0, nrow=row, ncol=col)
  Nrank <- rep(0, row)
  Xstar <- matrix(0, nrow=row, ncol=col)
  for (j in 1:col) {
    Xsorted[,j] <- sort(X[,j])
    Nrank[order(normals[,j])] <- seq(1, row, by=1)
    Xstar[,j] <- Xsorted[Nrank, j]
  }
  return(Xstar)
}

#replace present-day targets and present-day all-cause mortality, which are country specific

#######################
#PRESENT-DAY TARGETS###
#######################

#prevalence smear and symptom breakdowns
prev_cases <- 314 #282 w/ no history of tx, 308 not currently on tx
m_cases <- 103 #89 w/ no history of tx, 99 not currently on tx
s_cases <- 93 #info on current/tx history not given
ms_cases <- 45 #info on current/tx history not given
targets_all[["prop_m_all"]] <- m_cases/prev_cases
targets_all[["prop_s_all"]] <- s_cases/prev_cases
targets_all[["prop_ms"]] <- ms_cases/prev_cases
targets_all_lb[["prop_m_all"]] <- qbeta(p=0.025, shape1=m_cases, shape2=prev_cases-m_cases)
targets_all_lb[["prop_s_all"]] <- qbeta(p=0.025, shape1=s_cases, shape2=prev_cases-s_cases)
targets_all_lb[["prop_ms"]] <- qbeta(p=0.025, shape1=ms_cases, shape2=prev_cases-ms_cases)
targets_all_ub[["prop_m_all"]] <- qbeta(p=0.975, shape1=m_cases, shape2=prev_cases-m_cases)
targets_all_ub[["prop_s_all"]] <- qbeta(p=0.975, shape1=s_cases, shape2=prev_cases-s_cases)
targets_all_ub[["prop_ms"]] <- qbeta(p=0.975, shape1=ms_cases, shape2=prev_cases-ms_cases)

#smear-positive prevalence to notification ratio
notif_sp <- 15812 #new smear-positive pulmonary notifications, 2011 (new_sp). assume no children are sp - fairly consistent w/ limited Cambodia data. 
prev_sp <- 271*(1-0.04)/100000 #smear-positive prev from exec summary, adjusting for 4% currently on tx
prev_sp_lb <- 212*(1-0.04)/100000 #lb smear-positive prev from exec summary, adjusting for 4% currently on tx
prev_sp_ub <- 348*(1-0.04)/100000 #ub smear-positive prev from exec summary, adjusting for 4% currently on tx
prev_sp_samples <- rbinom(n=100000, size=18000, prob=prev_sp)/18000 #sample size chosen to match prev_lb and prev_ub
notif_adj <- rbeta(n=100000, shape1=57, shape2=16) #adjustment factor for those on tx that don't get reported - based on % tx in public sector from prev survey
pop_adult <- 9753000 #adult pop in 2011 from WPP
pnr_m_samples <- prev_sp_samples/((notif_sp/notif_adj)/pop_adult)
targets_all[["pnr_m_all"]] <- mean(pnr_m_samples)
targets_all_lb[["pnr_m_all"]] <- quantile(pnr_m_samples, 0.025)
targets_all_ub[["pnr_m_all"]] <- quantile(pnr_m_samples, 0.975)
#calculate quantities used in parameterizing the likelihood for the PNR target - using dampack gamma_params
sd_gamma <- (targets_all_ub[["pnr_m_all"]]-targets_all_lb[["pnr_m_all"]])/3.9 
pnr_gamma_shape <- gamma_params(targets_all[["pnr_m_all"]], sd_gamma)$shape 
pnr_gamma_scale <- gamma_params(targets_all[["pnr_m_all"]], sd_gamma)$scale
pnr_params <- list("pnr_gamma_shape"=pnr_gamma_shape, "pnr_gamma_scale"=pnr_gamma_scale)

#percent of notifications that are smear-positive
#note: proportion ends up being so high already that we don't try to remove kids from this - probably a very small % of pulmonary TB diagnoses anyway
notif_sn <- 7686 #new smear-negative pulmonary notifications, 2011 (new_sn) - these are clinical diagnoses
all_cxr_s <- 710 #total ppl screened w/ +CXR and +symptom screen
case_cxr_s <- 88 #bac-confirmed TB cases w/ +CXR and +symptom screen
prop_clindx_tb_lb <- case_cxr_s/all_cxr_s
prop_clindx_tb_mean <- 0.254 #Cambodia-specific, PPV from cough>2wks in table 3: from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018502#s3
prop_clindx_tb_ub <- 0.37 #ub from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018502#s3
prop_clindx_tb_samples <- rbeta(n=100000, shape1=25/3, shape2=(100-25)/3) #to fit mean and CI above
notif_sn_adj <- notif_sn*prop_clindx_tb_samples
prop_m_notif <- notif_sp/(notif_sn_adj+notif_sp)
#add random noise to achieve 5% widening of the 2.5 and 97.5th CIs
prop_m_notif <- prop_m_notif + rnorm(1000000, mean=0, sd=0.05) 
#truncate at 0% and 100% (this doesn't actually affect any samples)
prop_m_notif[prop_m_notif>1] <- 1
prop_m_notif[prop_m_notif<0] <- 0
#convert to prob of each probability
prop_m_notif_tmp <- data.frame(table(round(prop_m_notif*100))/length(prop_m_notif))
prop_m_notif_samples <- prop_m_notif_tmp$Freq
names(prop_m_notif_samples) <- prop_m_notif_tmp$Var1
#smooth
prop_m_notif_smooth <- rep(0, 101)
names(prop_m_notif_smooth) <- as.character(0:100)
prop_m_notif_smooth[names(prop_m_notif_samples)] <- prop_m_notif_samples
minprop <- min(as.integer(names(prop_m_notif_samples)))
maxprop <- max(as.integer(names(prop_m_notif_samples)))
prop_m_notif_smooth[prop_m_notif_smooth==0 & 
                      as.integer(names(prop_m_notif_smooth))>minprop &
                      as.integer(names(prop_m_notif_smooth))<maxprop] <- 0.000001
prop_m_notif_smooth[as.integer(names(prop_m_notif_smooth))<minprop] <- 
  0.000001/(1+minprop-as.integer(names(prop_m_notif_smooth[as.integer(names(prop_m_notif_smooth))<minprop])))
prop_m_notif_smooth[as.integer(names(prop_m_notif_smooth))>maxprop] <- 
  0.000001/(1+as.integer(names(prop_m_notif_smooth[as.integer(names(prop_m_notif_smooth))>maxprop]))-maxprop)
prop_m_notif_smooth[["100"]] <- pmin(prop_m_notif_smooth[["100"]], prop_m_notif_smooth[["99"]])
#update targets
targets_all[["prop_m_notif"]] <- mean(prop_m_notif)
targets_all_lb[["prop_m_notif"]] <- quantile(prop_m_notif, 0.025)
targets_all_ub[["prop_m_notif"]] <- quantile(prop_m_notif, 0.975)

#untreated mortality to prevalence ratio
#treated CFR - use lognormal distribution given skew
cfr_samples <- rlnorm(n=100000, meanlog=log(targets_prev %>% filter(param=="treat_tb_cfr") %>% pull(value)), 
                      sdlog=0.5)
#total estimated deaths
deaths_samples <- rnorm(n=100000, 
                        mean=targets_prev %>% filter(param=="deaths_all_count") %>% pull(value), 
                        sd=(targets_prev %>% filter(param=="deaths_all_count_ub") %>% pull(value) - 
                              targets_prev %>% filter(param=="deaths_all_count_lb") %>% pull(value))/4) #total estimated cases (from WHO, so consistent with deaths estimates)
#estimated incident cases
case_samples <- rgamma(n=100000, 
                       shape=gamma_params(targets_prev %>% filter(param=="tb_inc_all") %>% pull(value), 
                                          (targets_prev %>% filter(param=="tb_inc_all_ub") %>% pull(value) -
                                             targets_prev %>% filter(param=="tb_inc_all_lb") %>% pull(value))/
                                            (4))$shape, 
                       scale=gamma_params(targets_prev %>% filter(param=="tb_inc_all") %>% pull(value), 
                                          (targets_prev %>% filter(param=="tb_inc_all_ub") %>% pull(value) -
                                             targets_prev %>% filter(param=="tb_inc_all_lb") %>% pull(value))/
                                            (4))$scale)
#estimate all ages prev (excluding those currently on tx)
prev_samples <- rnorm(n=100000, mean=targets_prev %>% filter(param=="prev_all_nt") %>% pull(value),
                      sd=(targets_prev %>% filter(param=="prev_all_nt_ub") %>% pull(value) -
                            targets_prev %>% filter(param=="prev_all_nt_lb") %>% pull(value))/(2*1.96))
#calculate numbers treated that didn't fail/LTFU and numbers treated that did fail/LTRU
cases_failLTFU <- targets_prev %>% filter(param=="failure_tx_count") %>% pull(value) +
  targets_prev %>% filter(param=="ltfu_tx_count") %>% pull(value)
cases_tx_nofailLTFU <- targets_prev %>% filter(param=="tx_count_all") %>% pull(value) - 
  cases_failLTFU
#estimate of % treated not captured in notifications that should still be removed from mortality estimate
prop_tx_no_notif_samples <- 1 - notif_adj
#combine all components to calculate targets
deaths_untx_samples <- deaths_samples - (1/(1-prop_tx_no_notif_samples))*
  (cfr_samples*cases_tx_nofailLTFU + (deaths_samples/case_samples)*cases_failLTFU)
deaths_untx_per_case_samples <- deaths_untx_samples/
  (prev_samples*targets_prev %>% filter(param=="pop") %>% pull(value)/1000)
deaths_untx_per_case <- mean(deaths_untx_per_case_samples)
deaths_untx_per_case_lb <- quantile(deaths_untx_per_case_samples, 0.025)[[1]]
deaths_untx_per_case_ub <- quantile(deaths_untx_per_case_samples, 0.975)[[1]]
#remove negative deaths
deaths_untx_per_case_samples <- deaths_untx_per_case_samples[deaths_untx_per_case_samples>=0]
#convert to prob of each number
mort_samples_tmp <- data.frame(table(round(deaths_untx_per_case_samples*1000))/length(deaths_untx_per_case_samples))
mort_samples <- mort_samples_tmp$Freq
names(mort_samples) <- mort_samples_tmp$Var1
#update targets
targets_all[["deaths_tb"]] <- deaths_untx_per_case
targets_all_lb[["deaths_tb"]] <- deaths_untx_per_case_lb
targets_all_ub[["deaths_tb"]] <- deaths_untx_per_case_ub

#present-day all-cause mortality
m_ac_annual <- 0.00655 #all-cause/non-TB mortality probability (applied to symptom- TB too - annual avg. mortality probability among 15-64 years in Cambodia)
m_ac <- 1-exp(-1*-log(1-m_ac_annual)/12)
#update params
params$m_ac <- m_ac
params_fixed_prev$m_ac <- m_ac

#save params and targets to RDA file for faster loading in the future/on MARCC
save(n, verbose, cyc_len, params, names_params_calib, names_params_calib_mult,
     params_calib_prev, params_fixed_prev, params_calib_prev_mult,
     params_calib_hist, params_fixed_hist, params_fixed_hist_mult, 
     params_depend, 
     priors_prev_lb, priors_prev_ub, priors_prev_mult_lb, priors_prev_mult_ub,
     priors_hist_lb, priors_hist_ub, mort_samples, prop_m_notif_smooth,
     targets_all, targets_all_lb, targets_all_ub, pnr_params,
     file="data/params_targets_cambodia.Rda")
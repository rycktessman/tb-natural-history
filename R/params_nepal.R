setwd("C:/Users/Tess/OneDrive - Johns Hopkins/TB/Natural History Modeling")
load("data/params_targets.Rda")
library(MASS)
library(tidyverse)
library(dampack)
library(readxl)

targets_prev <- read_xlsx("data/country_targets.xlsx", range="A1:F1000")[,c(1,6)]
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
prev_cases <- 225 #3 currently on tx, 22 history of treatment within past 2 years
m_cases <- 68 #info on current/tx history not given
s_cases <- 60 #info on current/tx history not given
ms_cases <- 0.375*m_cases #% of all smear+ that were symptomatic (slightly different but close to # of smear+ cases)
targets_all[["prop_m_all"]] <- m_cases/prev_cases
targets_all[["prop_s_all"]] <- s_cases/prev_cases
targets_all[["prop_ms"]] <- ms_cases/prev_cases
targets_all_lb[["prop_m_all"]] <- qbeta(p=0.025, shape1=m_cases, shape2=prev_cases-m_cases)
targets_all_lb[["prop_s_all"]] <- qbeta(p=0.025, shape1=s_cases, shape2=prev_cases-s_cases)
targets_all_lb[["prop_ms"]] <- qbeta(p=0.025, shape1=ms_cases, shape2=prev_cases-ms_cases)
targets_all_ub[["prop_m_all"]] <- qbeta(p=0.975, shape1=m_cases, shape2=prev_cases-m_cases)
targets_all_ub[["prop_s_all"]] <- qbeta(p=0.975, shape1=s_cases, shape2=prev_cases-s_cases)
targets_all_ub[["prop_ms"]] <- qbeta(p=0.975, shape1=ms_cases, shape2=prev_cases-ms_cases)

#prevalence to notification ratio
#[treatment-naive] prevalence from survey
prev <- 416.3*(1-15/225)/100000 #estimated prev from survey, adjusting for 15/225 on tx in last 2 yrs
prev_lb <- 314.1*(1-15/225)/100000 #estimated prev from survey, adjusting for 15/225 on tx in last 2 yrs
prev_ub <- 518.5*(1-15/225)/100000 #estimated prev from survey, adjusting for 15/225 on tx in last 2 yrs
#notifications from WHO data
notif_clindx <- 4286 #new pulmonary clinically diagnosed notifications (excludes relapse & EP), 2018 (new_clindx)
notif_labconf <- 15770 #new pulmonary lab-confirmed notifications (excludes relapse & EP), 2018 (new_labconf)
notif_child <- 880+507 #child notifications: newrel_m014+newrel_f014 - assume negligible amount of relapses among children
#population 
pop_adult <- 19555000 #adult pop in 2018 from WPP
#parameterize prevalence and notifications distribution
prev_samples <- rbinom(n=100000, size=15000, prob=prev)/15000 #sample size chosen to match lb and ub
prop_clindx_tb_lb <- 59/1621  #% of CXR+ symptom+ from the survey that were TB+
prop_clindx_tb_mean <- 85/360 #% of presumptive smear- TB cases that were Xpert+ in Nepal: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6937953/
prop_clindx_tb_samples <- rbeta(n=100000, shape1=85/35, shape2=(360-85)/35) #matches mean and lower bound above
prop_child_remove <- rbeta(n=100000, shape1=1, shape2=1) #wide uncertainty here - child notifications are very small
notif_samples <- notif_labconf + notif_clindx*prop_clindx_tb_samples - notif_child*prop_child_remove
pnr_samples <- prev_samples/(notif_samples/pop_adult)

#update targets - replace pnr_m_all with pnr_all
targets_all <- targets_all[names(targets_all)!="pnr_m_all"]
targets_all_lb <- targets_all_lb[names(targets_all_lb)!="pnr_m_all"]
targets_all_ub <- targets_all_ub[names(targets_all_ub)!="pnr_m_all"]
targets_all[["pnr_all"]] <- mean(pnr_samples)
targets_all_lb[["pnr_all"]] <- quantile(pnr_samples, 0.025)
targets_all_ub[["pnr_all"]] <- quantile(pnr_samples, 0.975)
#calculate quantities used in parameterizing the likelihood for the PNR target - using dampack gamma_params
sd_gamma <- (targets_all_ub[["pnr_all"]]-targets_all_lb[["pnr_all"]])/3.9 
pnr_gamma_shape <- gamma_params(targets_all[["pnr_all"]], sd_gamma)$shape 
pnr_gamma_scale <- gamma_params(targets_all[["pnr_all"]], sd_gamma)$scale
pnr_params <- list("pnr_gamma_shape"=pnr_gamma_shape, "pnr_gamma_scale"=pnr_gamma_scale)

#untreated mortality to prevalence ratio
#treated CFR - use lognormal distribution given skew
cfr_samples <- rlnorm(n=100000, meanlog=log(targets_prev %>% filter(param=="treat_tb_cfr") %>% pull(value)), 
                      sdlog=0.4)
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
prop_tx_no_notif_samples <- rbeta(n=100000, shape1=10, shape2=48) #48 out of 58 ppl currently on TB treatment said govt facility was first choice for tx
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


#percent of notifications that are smear-positive
notif_xpert_lb <- 4743*31885/33199  #newinc_rdx  WHO notifications spreadsheet, 2015
notif_xpert_ub <- 10199
notif_xpert_samples <- runif(n=100000, min=notif_xpert_lb, max=notif_xpert_ub) #assume somewhere between no increase and doubling of Xpert usage from 2015-18
pxpertpos_testtreat_samples <- rbeta(n=100000, shape1=5.5, shape2=0.69) #to match mean 90%, LB 57%, UB 100%
psmearpos_xpertpos_samples <- rbeta(n=100000, shape1=10, shape2=12.27) #to match mean 44%, LB 26%, UB 68%

pnosmear_clindx_samples <- rbeta(n=100000, shape1=10, shape2=90) #publications and pre-2012 data indicate most TB gets smear-tested, but added some uncertainty here
psmearpos_TBnosmearclindx_samples <- rbeta(n=100000, shape1=53/6, shape2=47/6) #lb=30% smear+ in prev survey, ub= 75% mean (notif_labconf - notif_xpert_smearneg / notif_TB)

notif_xpert_smearneg <- notif_xpert_samples*pxpertpos_testtreat_samples*(1-psmearpos_xpertpos_samples)
notif_smearpos_tested_samples <- notif_labconf - 
  notif_xpert_samples*pxpertpos_testtreat_samples*(1-psmearpos_xpertpos_samples)
notif_TB_tested_samples <- notif_labconf + 
  notif_clindx*prop_clindx_tb_samples*(1-pnosmear_clindx_samples)
psmearpos_TBnosmearclindx_samples_alt1 <- notif_smearpos_tested_samples/notif_TB_tested_samples #alternative upper bound
rho <- cor(psmearpos_TBnosmearclindx_samples_alt1, prop_clindx_tb_samples)
means <- rep(0, 2)
cov_pd <- diag(nrow=2, ncol=2)
cov_pd[cov_pd==0] <- rho
norms <- mvrnorm(1000000, means, cov_pd)
samples_sorted <- sorted_rank(cbind(prop_clindx_tb_samples, psmearpos_TBnosmearclindx_samples),  norms)
pTB_clindx_sorted <- samples_sorted[,1]
psmearpos_TBnosmearclindx_sorted <- samples_sorted[,2]

notif_clindx_smearpos <- notif_clindx*pnosmear_clindx_samples*pTB_clindx_sorted*psmearpos_TBnosmearclindx_sorted
notif_smearpos <- notif_labconf - notif_xpert_smearneg + notif_clindx_smearpos
notif_TB <- notif_labconf + notif_clindx*prop_clindx_tb_samples

prop_m_notif <- notif_smearpos/notif_TB
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
#update targets
targets_all[["prop_m_notif"]] <- mean(prop_m_notif)
targets_all_lb[["prop_m_notif"]] <- quantile(prop_m_notif, 0.025)
targets_all_ub[["prop_m_notif"]] <- quantile(prop_m_notif, 0.975)

#present-day all-cause mortality
m_ac_annual <- 0.00483 #all-cause/non-TB mortality probability (applied to symptom- TB too - annual avg. mortality probability among 15-64 years in Cambodia)
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
     file="data/params_targets_nepal.Rda")
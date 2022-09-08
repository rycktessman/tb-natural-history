#generate RDA file with target info

setwd("~/GitHub/tb-natural-history")
library(MASS)
library(tidyverse)
library(dampack)
library(readxl)

targets_all <- c()
targets_all_lb <- c()
targets_all_ub <- c()

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

#######################
#PRESENT-DAY TARGETS###
#######################

#prevalence smear and symptom breakdowns
prev_cases <- 116+25+25+53 #from doc shared by Hoa and Dr. Hai
m_cases <- 25+25
s_cases <- 53+25
ms_cases <- 25
m_only_cases <- m_cases - ms_cases
s_only_cases <- s_cases - ms_cases
none_cases <- prev_cases - (m_only_cases + s_only_cases + ms_cases)
targets_all[["prop"]] <- none_cases/prev_cases
targets_all[["prop_m"]] <- m_only_cases/prev_cases
targets_all[["prop_s"]] <- s_only_cases/prev_cases
targets_all[["prop_ms"]] <- ms_cases/prev_cases
targets_all_lb[["prop"]] <- qbeta(p=0.025, shape1=none_cases, shape2=prev_cases-none_cases)
targets_all_lb[["prop_m"]] <- qbeta(p=0.025, shape1=m_only_cases, shape2=prev_cases-m_only_cases)
targets_all_lb[["prop_s"]] <- qbeta(p=0.025, shape1=s_only_cases, shape2=prev_cases-s_only_cases)
targets_all_lb[["prop_ms"]] <- qbeta(p=0.025, shape1=ms_cases, shape2=prev_cases-ms_cases)
targets_all_ub[["prop"]] <- qbeta(p=0.975, shape1=none_cases, shape2=prev_cases-none_cases)
targets_all_ub[["prop_m"]] <- qbeta(p=0.975, shape1=m_only_cases, shape2=prev_cases-m_only_cases)
targets_all_ub[["prop_s"]] <- qbeta(p=0.975, shape1=s_only_cases, shape2=prev_cases-s_only_cases)
targets_all_ub[["prop_ms"]] <- qbeta(p=0.975, shape1=ms_cases, shape2=prev_cases-ms_cases)

#prevalence to notification ratio
prev <- 322*(1-0.054)/100000 #from PLOS ONE paper, adjusting for 5.4% on tx in last 2 yrs
prev_lb <- 260*(1-0.054)/100000 #from PLOS ONE paper, adjusting for 5.4% on tx in last 2 yrs
prev_ub <- 399*(1-0.054)/100000 #from PLOS ONE paper, adjusting for 5.4% on tx in last 2 yrs
notif_clindx <- (23894+23117)/2 #new pulmonary clinically diagnosed notifications (excludes relapse & EP), avg 2017 and 2018 (new_clindx)
notif_labconf <- (50726+49720)/2#new pulmonary lab-confirmed notifications (excludes relapse & EP), avg 2017 and 2018 (new_labconf)
notif_child <- (966+854+750+802)/2 #child notifications: newrel_m014+newrel_f014/c_newinc - assume negligible amount of relapses among children
pop_adult <- (72749000+73405000)/2 #avg of adult pop in 17 and 18 from WPP
#parameterize prevalence and notifications distribution
prev_samples <- rbinom(n=100000, size=24000, prob=prev)/24000
prop_clindx_tb_samples <- rbeta(n=100000, shape1=103/20, shape2=(299-103)/20) #from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7157739/, adding more uncertainty
prop_child_remove <- rbeta(n=100000, shape1=1, shape2=1) #wide uncertainty here - child notifications are very small
notif_samples <- notif_labconf + notif_clindx*prop_clindx_tb_samples - notif_child*prop_child_remove
pnr_samples <- prev_samples/(notif_samples/pop_adult)
#add to targets 
targets_all[["pnr_all"]] <- mean(pnr_samples)
targets_all_lb[["pnr_all"]] <- quantile(pnr_samples, 0.025)
targets_all_ub[["pnr_all"]] <- quantile(pnr_samples, 0.975)
#calculate quantities used in parameterizing the likelihood for the PNR target - using dampack gamma_params
sd_gamma <- (targets_all_ub[["pnr_all"]]-targets_all_lb[["pnr_all"]])/3.9 
pnr_gamma_shape <- gamma_params(targets_all[["pnr_all"]], sd_gamma)$shape 
pnr_gamma_scale <- gamma_params(targets_all[["pnr_all"]], sd_gamma)$scale
pnr_params <- list("pnr_gamma_shape"=pnr_gamma_shape, "pnr_gamma_scale"=pnr_gamma_scale)

#mortality target: parameterize distributions around each component
#1. treated CFR - use lognormal distribution given skew
cfr_samples <- rlnorm(n=100000, meanlog=log(0.0244), sdlog=0.5) #from treatment outcomes reported to WHO
#2. total deaths estimate = overall CFR numerator: deaths - use normal distribution (matches WHO CI)
deaths_samples <- rnorm(n=100000, mean=13500, sd=(18000-9150)/(2*1.96)) #WHO-estimated total TB deaths
#3. overall CFR denominator: incident cases - use gamma distribution
inc_samples <- rgamma(n=100000, shape=gamma_params(176000,(254000-112000)/(2*1.96))$shape, 
                      scale=gamma_params(176000,(254000-112000)/(2*1.96))$scale)
#4. all ages untreated prevalence - estimated, parameterize normal distribution
pop <- 95073000
pop_prop_15plus <- 73077000/pop #% of pop that is 15+
notif_ratio_15plus <- (101192-99392)/101192 #ratio of TB notifications that are < 15 to 15+
ontx_prop <- 0.023 #% of prevalent TB cases that were on tx
prev <- 322/100000 #prev. mean and CI reported in survey among 15+
prev_lb <- 260/100000
prev_ub <- 399/100000
prev_adj <- prev*(1-ontx_prop)*(pop_prop_15plus+(1-pop_prop_15plus)*notif_ratio_15plus) #adjust for treatment history & estimate rel.# cases among < 15
prev_adj_lb <- prev_lb*(1-ontx_prop)*(pop_prop_15plus+(1-pop_prop_15plus)*notif_ratio_15plus)
prev_adj_ub <- prev_ub*(1-ontx_prop)*(pop_prop_15plus+(1-pop_prop_15plus)*notif_ratio_15plus)
prev_samples <- rnorm(n=100000, mean=prev_adj, sd=(prev_adj_ub-prev_adj_lb)/(2*1.96))
cases_failLTFU <- 2527 + 596 #numbers treated that didn't fail/LTFU and numbers treated that did fail/LTRU
cases_tx_nofailLTFU <- 100876 - cases_failLTFU
#estimate of % treated not captured in notifications that should still be removed from mortality estimate
prop_tx_no_notif_samples <- rbeta(n=100000, shape1=12, shape2=100-12)
#combine all components to calculate targets
deaths_untx_samples <- deaths_samples - (1/(1-prop_tx_no_notif_samples))*(cfr_samples*cases_tx_nofailLTFU + 
                                                                            (deaths_samples/inc_samples)*cases_failLTFU)
deaths_untx_per_case_samples <- deaths_untx_samples/(prev_samples*pop)
deaths_untx_per_case_samples <- deaths_untx_per_case_samples[deaths_untx_per_case_samples>=0] #throw out negatives (very few if any)
deaths_untx_per_case <- mean(deaths_untx_per_case_samples)
deaths_untx_per_case_lb <- quantile(deaths_untx_per_case_samples, 0.025)[[1]]
deaths_untx_per_case_ub <- quantile(deaths_untx_per_case_samples, 0.975)[[1]]
#convert to prob of each number
mort_samples_tmp <- data.frame(table(round(deaths_untx_per_case_samples*1000))/length(deaths_untx_per_case_samples))
mort_samples <- mort_samples_tmp$Freq
names(mort_samples) <- mort_samples_tmp$Var1
#update targets
targets_all[["deaths_tb"]] <- deaths_untx_per_case
targets_all_lb[["deaths_tb"]] <- deaths_untx_per_case_lb
targets_all_ub[["deaths_tb"]] <- deaths_untx_per_case_ub

#percent of notifications that are smear-positive
#adjust for smear-neg Xpert positives among lab-confirmed notifications
notif_xpert <- (26255+19916)/2 #newinc_rdx from 2017 and 2018 in WHO notifications spreadsheet
pxpertpos_testtreat_samples <- rbeta(n=100000, shape1=5.5, shape2=0.69) #to match mean 90%, LB 57%, UB 100%
psmearpos_xpertpos_samples <- rbeta(n=100000, shape1=24, shape2=16)
#adjust for clinical diagnoses that aren't truly TB
pnosmear_clindx_samples <- rbeta(n=100000, shape1=46/5, shape2=(257-46)/5) #based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7157739/
psmearpos_TBnosmearclindx_samples <- rbeta(n=100000, shape1=35/4.5, shape2=44/4.5) #lb=smear+ in prev survey, ub=mean on (notif_labconf - notif_xpert_smearneg / notif_TB)
#induce negative correlation between psmearpos_TBnosmearclindx_samples and pTB_nosmearclindx_samples based on correlation below
notif_xpert_smearneg <- notif_xpert*pxpertpos_testtreat_samples*(1-psmearpos_xpertpos_samples)
notif_smearpos_tested_samples <- notif_labconf - 
  notif_xpert*pxpertpos_testtreat_samples*(1-psmearpos_xpertpos_samples)
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
#calculate resulting estimates of notifications that were clinically diagnosed but are actually smear+
notif_clindx_smearpos <- notif_clindx*pnosmear_clindx_samples*pTB_clindx_sorted*psmearpos_TBnosmearclindx_sorted
notif_smearpos <- notif_labconf - notif_xpert_smearneg + notif_clindx_smearpos
notif_TB <- notif_labconf + notif_clindx*prop_clindx_tb_samples
prop_m_notif <- notif_smearpos/notif_TB
#add random noise to achieve 5% widening of the 2.5 and 97.5th CIs given uncertainty in this target
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
m_ac <- 1-exp(-1*-log(1-m_ac_annual)/12)
#update params
params$m_ac <- m_ac
params_fixed_prev$m_ac <- m_ac

#combine w/ historical cohort targets
targets_all <- c(targets_all, 
                 "tb_ms_dead_5yr"=0.5763, #percent symptomatic smear+ dead after 5 years (historical cohort) - from stata output
                 "tb_s_dead_5yr"=0.1272, #percent symptomatic smear- dead after 5 years (historical cohort) - from stata output
                 "tb_ms_dead_10yr"=0.7096, #percent symptomatic smear+ dead after 10 years (historical cohort) - from stata output
                 "tb_s_dead_10yr"=0.2105 #percent symptomatic smear- dead after 10 years (historical cohort) - from stata output
)
targets_all_lb <- c(targets_all_lb,
                    "tb_ms_dead_5yr"=0.51, #percent symptomatic smear+ dead after 5 years (historical cohort) - from stata output
                    "tb_s_dead_5yr"=0.09, #percent symptomatic smear- dead after 5 years (historical cohort) - from stata output
                    "tb_ms_dead_10yr"=0.65, #percent symptomatic smear+ dead after 10 years (historical cohort) - from stata output
                    "tb_s_dead_10yr"=0.15 #percent symptomatic smear- dead after 10 years (historical cohort) - from stata output
)
targets_all_ub <- c(targets_all_ub,
                    "tb_ms_dead_5yr"=0.64, #percent symptomatic smear+ dead after 5 years (historical cohort) - from stata output
                    "tb_s_dead_5yr"=0.16, #percent symptomatic smear- dead after 5 years (historical cohort) - from stata output
                    "tb_ms_dead_10yr"=0.77, #percent symptomatic smear+ dead after 10 years (historical cohort) - from stata output
                    "tb_s_dead_10yr"=0.26 #percent symptomatic smear- dead after 10 years (historical cohort) - from stata output
)

#save params and targets to RDA file for faster loading in the future/on MARCC
save(mort_samples, prop_m_notif_smooth,
     targets_all, targets_all_lb, targets_all_ub, 
     pnr_params, prev_cases,
     file="data/targets_vietnam.Rda")
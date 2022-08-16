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
prev_cases <- 314 #282 w/ no history of tx, 308 not currently on tx
m_cases <- 103 #89 w/ no history of tx, 99 not currently on tx
s_cases <- 93 #info on current/tx history not given
ms_cases <- 45 #info on current/tx history not given
m_only_cases <- m_cases - ms_cases
s_only_cases <- s_cases - ms_cases
targets_all[["prop_m"]] <- m_only_cases/prev_cases
targets_all[["prop_s"]] <- s_only_cases/prev_cases
targets_all[["prop_ms"]] <- ms_cases/prev_cases
targets_all_lb[["prop_m"]] <- qbeta(p=0.025, shape1=m_only_cases, shape2=prev_cases-m_only_cases)
targets_all_lb[["prop_s"]] <- qbeta(p=0.025, shape1=s_only_cases, shape2=prev_cases-s_only_cases)
targets_all_lb[["prop_ms"]] <- qbeta(p=0.025, shape1=ms_cases, shape2=prev_cases-ms_cases)
targets_all_ub[["prop_m"]] <- qbeta(p=0.975, shape1=m_only_cases, shape2=prev_cases-m_only_cases)
targets_all_ub[["prop_s"]] <- qbeta(p=0.975, shape1=s_only_cases, shape2=prev_cases-s_only_cases)
targets_all_ub[["prop_ms"]] <- qbeta(p=0.975, shape1=ms_cases, shape2=prev_cases-ms_cases)

#smear-positive prevalence to notification ratio
notif_m_naive <- 15812 #new smear-positive pulmonary notifications, 2011
prev_m_naive <- 271*(1-0.04)/100000 #smear-positive prev from survey, adjusting for 4% currently on tx
prev_m_naive_lb <- 212*(1-0.04)/100000 #lb smear-positive prev from survey, adjusting for 4% currently on tx
prev_m_naive_ub <- 348*(1-0.04)/100000 #ub smear-positive prev from survey, adjusting for 4% currently on tx
prev_m_naive_samples <- rbinom(n=100000, size=18000, prob=prev_m_naive)/18000 #sample size chosen to match prev_lb and prev_ub
notif_adj <- rbeta(n=100000, shape1=57, shape2=16) #adjustment factor for those on tx that don't get reported - based on % tx in public sector from prev survey
pop_adult <- 9753000 #adult pop in 2011 from WPP
pnr_m_samples <- prev_m_naive_samples/((notif_m_naive/notif_adj)/pop_adult)
#add to targets
targets_all[["pnr_m_all"]] <- mean(pnr_m_samples)
targets_all_lb[["pnr_m_all"]] <- quantile(pnr_m_samples, 0.025)
targets_all_ub[["pnr_m_all"]] <- quantile(pnr_m_samples, 0.975)
#calculate quantities used in parameterizing the likelihood for the PNR target - using dampack gamma_params
sd_gamma <- (targets_all_ub[["pnr_m_all"]]-targets_all_lb[["pnr_m_all"]])/3.9 
pnr_gamma_shape <- gamma_params(targets_all[["pnr_m_all"]], sd_gamma)$shape 
pnr_gamma_scale <- gamma_params(targets_all[["pnr_m_all"]], sd_gamma)$scale
pnr_params <- list("pnr_gamma_shape"=pnr_gamma_shape, "pnr_gamma_scale"=pnr_gamma_scale)

#percent of notifications that are smear-positive 
#split already given in 2011, so fewer adjustments required
notif_sp <- 15812
notif_sn <- 7686 #new smear-negative pulmonary notifications, 2011 (new_sn) - these are clinical diagnoses
all_cxr_s <- 710 #total ppl screened w/ +CXR and +symptom screen
case_cxr_s <- 88 #bac-confirmed TB cases w/ +CXR and +symptom screen
#adjust for % clinical diagnoses that aren't truly TB
prop_clindx_tb_lb <- case_cxr_s/all_cxr_s
prop_clindx_tb_mean <- 0.254 #Cambodia-specific, PPV from cough>2wks in table 3: from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018502#s3
prop_clindx_tb_ub <- 0.37 #ub from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018502#s3
prop_clindx_tb_samples <- rbeta(n=100000, shape1=25/3, shape2=(100-25)/3) #to fit mean and CI above
notif_sn_adj <- notif_sn*prop_clindx_tb_samples
prop_m_notif <- notif_sp/(notif_sn_adj+notif_sp)
#add random noise to achieve 5% widening of the 2.5 and 97.5th CIs given uncertainty in this target
prop_m_notif <- prop_m_notif + rnorm(1000000, mean=0, sd=0.05) 
#truncate at 0% and 100% 
prop_m_notif[prop_m_notif>1] <- 1
prop_m_notif[prop_m_notif<0] <- 0
#convert to prob of each probability and smooth to use as empirical distribution
prop_m_notif_tmp <- data.frame(table(round(prop_m_notif*100))/length(prop_m_notif))
prop_m_notif_samples <- prop_m_notif_tmp$Freq
names(prop_m_notif_samples) <- prop_m_notif_tmp$Var1
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

#mortality target: parameterize distributions around each component
#1. treated CFR - use lognormal distribution given skew
cfr_samples <- rlnorm(n=100000, meanlog=log(0.0232), sdlog=0.5) #from treatment outcomes reported to WHO
#2. total deaths estimate = overall CFR numerator: deaths - use normal distribution (matches WHO CI)
deaths_samples <- rnorm(n=100000, mean=4600, sd=(5900-2600)/(2*1.96)) #WHO-estimated total TB deaths
#3. overall CFR denominator: incident cases - use gamma distribution
inc_samples <- rgamma(n=100000, shape=gamma_params(62000,(88000-40000)/(2*1.96))$shape, 
                      scale=gamma_params(62000,(88000-40000)/(2*1.96))$scale)
#4. all ages untreated prevalence - estimated, parameterize normal distribution
pop <- 14540000
pop_prop_15plus <- 9753000/pop #% of pop that is 15+
notif_ratio_15plus <- (35169-28294)/28294 #ratio of TB notifications that are < 15 to 15+
ontx_prop <- 6/314 #% of prevalent TB cases that were on tx
prev <- 831/100000 #prev. mean and CI reported in survey among 15+
prev_lb <- 707/100000
prev_ub <- 977/100000
prev_adj <- prev*(1-ontx_prop)*(pop_prop_15plus+(1-pop_prop_15plus)*notif_ratio_15plus) #adjust for treatment history & estimate rel.# cases among < 15
prev_adj_lb <- prev_lb*(1-ontx_prop)*(pop_prop_15plus+(1-pop_prop_15plus)*notif_ratio_15plus)
prev_adj_ub <- prev_ub*(1-ontx_prop)*(pop_prop_15plus+(1-pop_prop_15plus)*notif_ratio_15plus)
prev_samples <- rnorm(n=100000, mean=prev_adj, sd=(prev_adj_ub-prev_adj_lb)/(2*1.96))
cases_failLTFU <- 82 + 1439 #numbers treated that didn't fail/LTFU and numbers treated that did fail/LTRU
cases_tx_nofailLTFU <- 37809 - cases_failLTFU
#combine all components to calculate targets
deaths_untx_samples <- deaths_samples - (1/notif_adj)*(cfr_samples*cases_tx_nofailLTFU + 
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
     file="data/targets_cambodia.Rda")
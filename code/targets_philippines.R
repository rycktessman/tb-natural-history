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
prev_cases_m <- 359 #symptom status available; only those w/out history of/currently on treatment
prev_cases_s <- 356 #smear status available; only those w/out history of/currently on treatment
prev_cases <- 356 #use slightly smaller sample for likelihood functions (more uncertainty)
m_cases <- 131
s_cases <- 104 
ms_cases <- 61
m_only_cases <- m_cases - ms_cases
s_only_cases <- s_cases - ms_cases
targets_all[["prop_m"]] <- m_only_cases/prev_cases_m
targets_all[["prop_s"]] <- s_only_cases/prev_cases_s
targets_all[["prop_ms"]] <- ms_cases/prev_cases_s
targets_all_lb[["prop_m"]] <- qbeta(p=0.025, shape1=m_only_cases, shape2=prev_cases_m-m_only_cases)
targets_all_lb[["prop_s"]] <- qbeta(p=0.025, shape1=s_only_cases, shape2=prev_cases_s-s_only_cases)
targets_all_lb[["prop_ms"]] <- qbeta(p=0.025, shape1=ms_cases, shape2=prev_cases_s-ms_cases)
targets_all_ub[["prop_m"]] <- qbeta(p=0.975, shape1=m_only_cases, shape2=prev_cases_m-m_only_cases)
targets_all_ub[["prop_s"]] <- qbeta(p=0.975, shape1=s_only_cases, shape2=prev_cases_s-s_only_cases)
targets_all_ub[["prop_ms"]] <- qbeta(p=0.975, shape1=ms_cases, shape2=prev_cases_s-ms_cases)

#prevalence to notification ratio - adjusted to remove those w/ history of/currently on treatment
prop_tx <- 0.24 #% of smear+ in prev survey that were already/had already been on treatment
prev_m_naive <- 4.34*(1-prop_tx) #prevalence of smear+ TB reported in prev survey (among 15+ only), removing those on tx
prev_m_naive_lb <- 3.5*(1-prop_tx) #lower bound reported, removing those on tx
prev_m_naive_ub <- 5.18*(1-prop_tx) #upper bound reported, removing those on tx
notif_m <- 4.34/3.1 #number of smear+ notifications implied by prev survey PNR and prevalence
notif_m_naive <- notif_m*(106516-12850)/106516 #remove % lab conf cases that were return/relapse (from notification data)
#add uncertainty related to % notifications reported
notif_adj <- 1/0.75 #75% is based on % in survey that claim they would seek treatment in the public sector (and is halfway between 50% below and 100% upper bound)
notif_adj_lb <- 1/0.5 #50% is based on % currently treated from survey that were in notifications data
notif_adj_ub <- 1 #upper bound - notifications represent all true cases treated
targets_all[["pnr_m_all"]] <- (prev_m_naive/notif_m_naive)/notif_adj
targets_all_lb[["pnr_m_all"]] <- (prev_m_naive_lb/notif_m_naive)/notif_adj_lb
targets_all_ub[["pnr_m_all"]] <- (prev_m_naive_ub/notif_m_naive)/notif_adj_ub

#calculate quantities used in parameterizing the likelihood for the PNR target - using dampack gamma_params
sd_gamma <- (targets_all_ub[["pnr_m_all"]]-targets_all_lb[["pnr_m_all"]])/3.9 
pnr_gamma_shape <- gamma_params(targets_all[["pnr_m_all"]], sd_gamma)$shape 
pnr_gamma_scale <- gamma_params(targets_all[["pnr_m_all"]], sd_gamma)$scale
pnr_params <- list("pnr_gamma_shape"=pnr_gamma_shape, "pnr_gamma_scale"=pnr_gamma_scale)

#mortality target: parameterize distributions around each component
#1. treated CFR - use lognormal distribution given skew
cfr_samples <- rlnorm(n=100000, meanlog=log(0.0243), sdlog=0.5) #from treatment outcomes reported to WHO
#2. total deaths estimate = overall CFR numerator: deaths - use normal distribution (matches WHO CI)
deaths_samples <- rnorm(n=100000, mean=29000, sd=(33000-25000)/(2*1.96)) #WHO-estimated total TB deaths
#3. overall CFR denominator: incident cases - use gamma distribution
inc_samples <- rgamma(n=100000, shape=gamma_params(574000,(898000-322000)/(2*1.96))$shape, 
                      scale=gamma_params(574000,(898000-322000)/(2*1.96))$scale)
#4. all ages untreated prevalence - estimated, parameterize normal distribution
pop <- 103664000
pop_prop_15plus <- 70573000/pop #% of pop that is 15+
notif_ratio_15plus <- (332941-284242)/284242 #ratio of TB notifications that are < 15 to 15+
ontx_prop <- 0.0558 #% of prevalent TB cases that were on tx
prev <- 1159/100000 #prev. mean and CI reported in survey among 15+
prev_lb <- 1016/100000
prev_ub <- 1301/100000
prev_adj <- prev*(1-ontx_prop)*(pop_prop_15plus+(1-pop_prop_15plus)*notif_ratio_15plus) #adjust for treatment history & estimate rel.# cases among < 15
prev_adj_lb <- prev_lb*(1-ontx_prop)*(pop_prop_15plus+(1-pop_prop_15plus)*notif_ratio_15plus)
prev_adj_ub <- prev_ub*(1-ontx_prop)*(pop_prop_15plus+(1-pop_prop_15plus)*notif_ratio_15plus)
prev_samples <- rnorm(n=100000, mean=prev_adj, sd=(prev_adj_ub-prev_adj_lb)/(2*1.96))
cases_failLTFU <- 1480 + 13105 #numbers treated that didn't fail/LTFU and numbers treated that did fail/LTRU
cases_tx_nofailLTFU <- 332308 - cases_failLTFU
notif_adj <- 1/0.75 #CI: 1/0.5 - 1
#combine all components to calculate targets
deaths_untx_samples <- deaths_samples - notif_adj*(cfr_samples*cases_tx_nofailLTFU + 
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

#percent of notifications that are smear+
notif_labconf <- 115200 #no uncertainty
notif_xpert <- 41759 #no uncertainty
notif_clindx <- 200572
#adjust for bacteriologically-confirmed that are smear negative (from Xpert)
pxpertpos_testtreat_samples <- rbeta(n=1000000, shape1=5.5, shape2=0.69) #to match mean 90%, LB 57%, UB 100%
psmearpos_xpertpos_samples <- rbeta(n=1000000, shape1=10, shape2=12.27) #to match mean 44%, LB 26%, UB 68%
notif_xpert_smearneg <- notif_xpert*pxpertpos_testtreat_samples*(1-psmearpos_xpertpos_samples)
#adjust for private sector
pprivate_samples <- rbeta(n=1000000, shape1=6.5, shape2=50) #to match mean 10-11%, LB 5-6%, UB 20-21% 
pnosmear_privateclindx_samples <- rbeta(n=1000000, shape1=22, shape2=(38-22)) #to match binomial CI from prev survey 58% [42-74%]
ppublic_samples <- 1-pprivate_samples
pnosmear_publicclindx_samples <- rnorm(n=1000000, mean=0.19, sd=0.19/2)
pnosmear_publicclindx_samples[pnosmear_publicclindx_samples<0] <- 0 #truncate at 0%
pTB_clindx_samples <- rbeta(n=1000000, shape1=3.4, shape2=7) #30ish% [10-61%] 61% comes from 61% fpr_clindx_est above and adding noise
psmearpos_TBnosmearclindx_samples <- rbeta(n=1000000, shape1=54/3, shape2=(100-54)/3)
#induce negative correlation between psmearpos_TBnosmearclindx_samples and pTB_nosmearclindx_samples based on correlation below
notif_smearpos_tested_samples <- notif_labconf - 
  notif_xpert*pxpertpos_testtreat_samples*(1-psmearpos_xpertpos_samples)
notif_TB_tested_samples <- notif_labconf + 
  notif_clindx*pTB_clindx_samples*((1-pnosmear_privateclindx_samples)*pprivate_samples +
                                     (1-pnosmear_publicclindx_samples)*ppublic_samples)
psmearpos_TBnosmearclindx_samples_alt1 <- notif_smearpos_tested_samples/notif_TB_tested_samples #alternative upper bound
rho <- cor(psmearpos_TBnosmearclindx_samples_alt1, pTB_clindx_samples)
means <- rep(0, 2)
cov_pd <- diag(nrow=2, ncol=2)
cov_pd[cov_pd==0] <- rho
norms <- mvrnorm(1000000, means, cov_pd)
samples_sorted <- sorted_rank(cbind(pTB_clindx_samples, psmearpos_TBnosmearclindx_samples),  norms)
pTB_clindx_sorted <- samples_sorted[,1]
psmearpos_TBnosmearclindx_sorted <- samples_sorted[,2]
#calculate resulting estimates of notifications that were clinically diagnosed but are actually smear+
notif_private_clindx_smearpos <- notif_clindx*pprivate_samples*pnosmear_privateclindx_samples*
  pTB_clindx_sorted*psmearpos_TBnosmearclindx_sorted
notif_public_clindx_smearpos <- notif_clindx*ppublic_samples*pnosmear_publicclindx_samples*
  pTB_clindx_sorted*psmearpos_TBnosmearclindx_sorted
#calculate estimate of total smear-positive and total true TB notifications
notif_smearpos <- notif_labconf - notif_xpert_smearneg + notif_private_clindx_smearpos +
  notif_public_clindx_smearpos
notif_TB <- notif_labconf + notif_clindx*pTB_clindx_samples
prop_m_notif <- notif_smearpos/notif_TB
#add random noise to achieve 5% widening of the 2.5 and 97.5th CIs given uncertainty in this target
prop_m_notif <- prop_m_notif + rnorm(1000000, mean=0, sd=0.05) 
#truncate at 0% and 100% (this doesn't actually affect any samples)
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
#update targets
targets_all[["prop_m_notif"]] <- mean(prop_m_notif)
targets_all_lb[["prop_m_notif"]] <- quantile(prop_m_notif, 0.025)
targets_all_ub[["prop_m_notif"]] <- quantile(prop_m_notif, 0.975)

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

#all-cause mortality also varies by country
m_ac_annual <- 0.005 #all-cause/non-TB mortality probability (applied to symptom- TB too - annual avg. mortality probability among 15-64 years in Philippines)
m_ac <- 1-exp(-1*-log(1-m_ac_annual)/12) #convert annual background mortality probabilities to monthly

#save params and targets to RDA file for faster loading in the future/on MARCC
save(mort_samples, prop_m_notif_smooth,
     targets_all, targets_all_lb, targets_all_ub, 
     pnr_params, prev_cases,
     file="data/targets_philippines.Rda")
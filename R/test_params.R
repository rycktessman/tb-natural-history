#generate RDA file with parameter and target info

setwd("C:/Users/Tess/OneDrive - Johns Hopkins/TB/Natural History Modeling")
library(MASS) #used to induce correlations in targets
library(tidyverse) #load dplyr after MASS so that "select" works as expected
library(readxl)
library(dampack) #gamma_params used to define target distributions

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

#PARAMETERS: specify parameters - these correspond mostly to model diagram
#progression/regression probabilities: p=progress, r=regress, m=smear microscopy, s=symptoms, a's are multipliers
p_m <- 0.04 #probability of progression from smear- symptom- to smear+ symptom-
p_s <- 0.05 #probability of progression from smear- symptom- to smear- symptom+
a_p_m <- 2 #multiplier on probability of progression from smear- symptom+ to smear+ symptom+ (vs. smear- symptom- to smear+ symptom-)
a_p_s <- a_p_m #multiplier on probability of progression from smear+ symptom- to smear+ symptom+ (vs. smear- symptom- to smear- symptom+)
r_m <- 0.02 #probability of regression from smear+ symptom- to smear- symptom-
r_s <- 0.04 #probability of regression from smear- symptom+ to smear- symptom-
a_r_m <- 0.3 #multiplier on probability of regression from smear+ symptom+ to smear- symptom+ (vs. smear+ symptom- to smear- symptom-)
a_r_s <- a_r_m #multiplier on probability of regression from smear+ symptom+ to smear+ symptom- (vs. smear- symptom+ to smear- symptom-)
p_c <- 0 #probability of progression from spontaneous cure to smear- symptom-
#cure/diagnosis/treatment probabilities: c=cure, a's are multipliers
c_sp <- 0.005 #probability of spontaneous cure (from smear- symptom-)
c_tx <- 0.15 #probability of diagnosis and treatment (from smear- symptom+)
a_tx <- 2 #multiplier on probability of diagnosis and treatment if smear+ symptom+ (vs. smear- symptom+)
#mortality probabilities: m=mortality, a's are multipliers
m_ac_annual <- 0.005 #all-cause/non-TB mortality probability (applied to symptom- TB too - annual avg. mortality probability among 15-64 years in Philippines)
m_ac_hist_annual <- 0.01 #historical all-cause non-TB mortality probability in early 20th century Europe (from Ragonnet)
#convert annual background mortality probabilities to monthly
m_ac <- 1-exp(-1*-log(1-m_ac_annual)/12)
m_ac_hist <- 1-exp(-1*-log(1-m_ac_hist_annual)/12)
m_tb <- 0.1/12 #mortality probability from smear- symptom+ TB 
a_m <- a_tx #multiplier on mortality from symptom+ smear+TB (vs. symptom+ smear-)
#relative infectiousness
i <- 1 #infectiousness of smear- symptom-
i_m <- 2 #infectiousness of smear+ symptom-
i_s <- i #infectiousness of smear- symptom+
i_ms <- i_m #infectiousness of smear+ symptom+
#other parameters used in sensitivity analyses
#true number of TB cases treated relative to the number notified (> 1 = notifications are an overestimate, < 1 = notifications are an underestimate)
notif_adj <- 1/0.75 #75% is based on % in survey that claim they would seek treatment in the public sector (and is halfway between 50% below and 100% upper bound)
notif_adj_lb <- 1/0.5 #50% is based on % currently treated from survey that were in notifications data
notif_adj_ub <- 1 #upper bound - notifications represent all true cases treated
#other specifications
n <- 10000 #assume population of n people
verbose <- 1 #1/0 to display tracker for microsim model progress
cyc_len <- 1/12 #monthly time steps

#all params
params <- list(
  "p_m"=p_m,
  "p_s"=p_s,
  "a_p_m"=a_p_m,
  "a_p_s"=a_p_s,
  "r_m"=r_m,
  "r_s"=r_s,
  "a_r_m"=a_r_m,
  "a_r_s"=a_r_s,
  "p_c"=p_c,
  "c_sp"=c_sp,
  "c_tx"=c_tx,
  "a_tx"=a_tx,
  "m_ac"=m_ac,
  "m_tb"=m_tb,
  "a_m"=a_m,
  "i"=i,
  "i_m"=i_m,
  "i_s"=i_s,
  "i_ms"=i_ms,
  "inflows"=1
)

#names of parameters that are being calibrated
names_params_calib <- c("p_m"="Smear Progression",
                        "p_s"="Symptom Progression",
                        "a_p_m"="Progression Multiplier",
                        "r_m"="Smear Regression",
                        "r_s"="Symptom Regression",
                        "a_r_m"="Regression Multiplier",
                        "c_sp"="Spontaneous Cure",
                        "c_tx"="Treatment",
                        "a_m"="Treat. & Mort. Multiplier",
                        "m_tb"="TB Mortality")
names_params_calib_mult <- c(names_params_calib,
                             "a_tx"="Treatment Multiplier")
names_params_calib_mult["a_m"] <- "Mortality Multiplier"

#parameters that are being calibrated (dependent params removed) - calibration to prevalence data (starting values for optim)
params_calib_prev <- list(
  "p_m"=p_m,
  "p_s"=p_s,
  "a_p_m"=a_p_m,
  "r_m"=r_m,
  "r_s"=r_s,
  "a_r_m"=a_r_m,
  "c_sp"=c_sp,
  "c_tx"=c_tx,
  "a_m"=a_m,
  "m_tb"=m_tb
)

#expanded mult version has a_tx and a_m as 2 distinct parameters
params_calib_prev_mult <- params_calib_prev
params_calib_prev_mult[["a_tx"]] <- a_tx

#parameters that stay fixed (doesn't include dependent params) - calibration to prevalence data
params_fixed_prev <- list (
  "p_c"=p_c, #no spontaneous cure back to TB -  fixed for now
  "m_ac"=m_ac, #all-cause mortality is fixed
  "i"=i, #relative infectiousness is fixed
  "i_m"=i_m,
  "i_s"=i_s,
  "i_ms"=i_ms,
  "inflows"=1
)


#parameters that are being calibrated (dependent params removed) - calibration to historical cohort data (starting values for optim)
params_calib_hist <- params_calib_prev[names(params_calib_prev)!="c_tx"]

#parameters that stay fixed (doesn't include dependent params) - calibration to historical cohort data
params_fixed_hist <- list(
  "p_c"=p_c, #no spontaneous cure back to TB -  fixed for now
  "m_ac"=m_ac_hist, #all-cause mortality is fixed
  "i"=i, #relative infectiousness is fixed
  "i_m"=i_m,
  "i_s"=i_s,
  "i_ms"=i_ms,
  "inflows"=0, #no inflows in historical cohort version
  "c_tx"=0 #no treatment in historical cohort version
)
#need to add a_tx=0 to mult_expand=1 version
params_fixed_hist_mult <- params_fixed_hist
params_fixed_hist_mult[["a_tx"]] <- 0

#parameters that depend on other parameters (need to handle separately in calibration)
#same for both prevalence and cohort versions of calibration
params_depend <- list(
  "a_p_s"=a_p_s, 
  "a_r_s"=a_r_s,
  "a_tx"=a_tx
)

#PRIOR PARAMETER DISTRIBUTIONS - lower and upper bounds for uniform sampling
p_m_lb <- 0 #probability
p_m_ub <- 0.2
p_s_lb <- 0 #probability
p_s_ub <- 0.2
a_p_m_lb <- 1 #multiplier > 1
a_p_m_ub <- 10
r_m_lb <- 0 #probability
r_m_ub <- 0.2
r_s_lb <- 0 #probability
r_s_ub <- 0.5
a_r_m_lb <- 0 #multiplier < 1
a_r_m_ub <- 1
c_sp_lb <- 0 #probability
c_sp_ub <- 0.5
c_tx_lb <- 0 #probability
c_tx_ub <- 0.25
a_m_lb <- 1 #multiplier > 1
a_m_ub <- 20
a_tx_lb <- a_m_lb
a_tx_ub <- a_m_ub
m_tb_lb <- 0 #probability
m_tb_ub <- 0.1

priors_prev_lb <- list(
  "p_m"=p_m_lb,
  "p_s"=p_s_lb,
  "a_p_m"=a_p_m_lb,
  "r_m"=r_m_lb,
  "r_s"=r_s_lb,
  "a_r_m"=a_r_m_lb,
  "c_sp"=c_sp_lb,
  "c_tx"=c_tx_lb,
  "a_m"=a_m_lb,
  "m_tb"=m_tb_lb
)

priors_prev_ub <- list(
  "p_m"=p_m_ub,
  "p_s"=p_s_ub,
  "a_p_m"=a_p_m_ub,
  "r_m"=r_m_ub,
  "r_s"=r_s_ub,
  "a_r_m"=a_r_m_ub,
  "c_sp"=c_sp_ub,
  "c_tx"=c_tx_ub,
  "a_m"=a_m_ub,
  "m_tb"=m_tb_ub
)

priors_hist_lb <- priors_prev_lb[names(priors_prev_lb)!="c_tx"]
priors_hist_ub <- priors_prev_ub[names(priors_prev_ub)!="c_tx"]

priors_prev_mult_lb <- priors_prev_lb
priors_prev_mult_lb[["a_tx"]] <- a_tx_lb
priors_prev_mult_ub <- priors_prev_ub
priors_prev_mult_ub[["a_tx"]] <- a_tx_ub

#CALIBRATION TARGETS

#prevalence targets
targets_prev <- read_xlsx("data/country_targets.xlsx", range="A1:D1000")
#remove non-text params and change to numeric
#targets_prev <- targets_prev %>% filter(!str_detect(param, "details"))
targets_prev <- targets_prev %>% filter(target>0) #1=directly used as targets, 2=used to calculate target(s)
targets_prev <- targets_prev %>% mutate(philippines=as.numeric(philippines))

#calculate mortality target CI
#parameterize distributions around each component
set.seed(4)
#1. treated CFR - use lognormal distribution given skew
cfr_samples <- rlnorm(n=100000, meanlog=log(targets_prev %>% filter(param=="treat_tb_cfr") %>% pull(philippines)), 
                      sdlog=0.5)
#2. total deaths estimate = overall CFR numerator: deaths - use normal distribution
deaths_samples <- rnorm(n=100000, 
                        mean=targets_prev %>% filter(param=="deaths_all_count") %>% pull(philippines), 
                        sd=(targets_prev %>% filter(param=="deaths_all_count_ub") %>% pull(philippines) - 
                              targets_prev %>% filter(param=="deaths_all_count_lb") %>% pull(philippines))/(2*1.96))
#3. overall CFR denominator: notifications - use normal distribution
notif_samples <- rgamma(n=100000, 
                        shape=gamma_params(targets_prev %>% filter(param=="tb_inc_all") %>% pull(philippines), 
                                           (targets_prev %>% filter(param=="tb_inc_all_ub") %>% pull(philippines) -
                                              targets_prev %>% filter(param=="tb_inc_all_lb") %>% pull(philippines))/
                                             (2*1.96))$shape, 
                        scale=gamma_params(targets_prev %>% filter(param=="tb_inc_all") %>% pull(philippines), 
                                           (targets_prev %>% filter(param=="tb_inc_all_ub") %>% pull(philippines) -
                                              targets_prev %>% filter(param=="tb_inc_all_lb") %>% pull(philippines))/
                                             (2*1.96))$scale)
#4. all ages prevalence
prev_samples <- rnorm(n=100000, mean=targets_prev %>% filter(param=="prev_all_nt") %>% pull(philippines),
                      sd=(targets_prev %>% filter(param=="prev_all_nt_ub") %>% pull(philippines) -
                            targets_prev %>% filter(param=="prev_all_nt_lb") %>% pull(philippines))/(2*1.96))
#calculate numbers treated that didn't fail/LTFU and numbers treated that did fail/LTRU
cases_failLTFU <- targets_prev %>% filter(param=="failure_tx_count") %>% pull(philippines) +
  targets_prev %>% filter(param=="ltfu_tx_count") %>% pull(philippines)
cases_tx_nofailLTFU <- targets_prev %>% filter(param=="tx_count_all") %>% pull(philippines) - cases_failLTFU
#combine all components to calculate targets
deaths_untx_samples <- deaths_samples - notif_adj*(cfr_samples*cases_tx_nofailLTFU + 
                                             (deaths_samples/notif_samples)*cases_failLTFU)
deaths_untx_per_case_samples <- deaths_untx_samples/
  (prev_samples*targets_prev %>% filter(param=="pop") %>% pull(philippines)/1000)

deaths_untx_per_case <- mean(deaths_untx_per_case_samples)
deaths_untx_per_case_lb <- quantile(deaths_untx_per_case_samples, 0.025)[[1]]
deaths_untx_per_case_ub <- quantile(deaths_untx_per_case_samples, 0.975)[[1]]
#throw out negatives
deaths_untx_per_case_samples <- deaths_untx_per_case_samples[deaths_untx_per_case_samples>=0]
#convert to prob of each number
mort_samples_tmp <- data.frame(table(round(deaths_untx_per_case_samples*1000))/length(deaths_untx_per_case_samples))
mort_samples <- mort_samples_tmp$Freq
names(mort_samples) <- mort_samples_tmp$Var1

#calculate smear+ notifications target CI (new target)
notif_labconf <- 115200 #no uncertainty
notif_xpert <- 41759 #no uncertainty
notif_clindx <- 200572
#notif_trueTB_est <- 250000*(notif_labconf+notif_clindx)/330000 #from https://academic.oup.com/jid/article/216/suppl_7/S740/4595556
#fpr_clindx_est <- (notif_labconf + notif_clindx - notif_trueTB_est)/notif_clindx #assume all FPs are clinical diagnosis
pxpertpos_testtreat_samples <- rbeta(n=1000000, shape1=5.5, shape2=0.69) #to match mean 90%, LB 57%, UB 100%
psmearpos_xpertpos_samples <- rbeta(n=1000000, shape1=10, shape2=12.27) #to match mean 44%, LB 26%, UB 68%
notif_xpert_smearneg <- notif_xpert*pxpertpos_testtreat_samples*(1-psmearpos_xpertpos_samples)


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

notif_private_clindx_smearpos <- notif_clindx*pprivate_samples*pnosmear_privateclindx_samples*
  pTB_clindx_sorted*psmearpos_TBnosmearclindx_sorted
notif_public_clindx_smearpos <- notif_clindx*ppublic_samples*pnosmear_publicclindx_samples*
  pTB_clindx_sorted*psmearpos_TBnosmearclindx_sorted

#psmearpos_clindxprivate <- pnosmear_privateclindx_samples*pTB_clindx_sorted*psmearpos_TBnosmearclindx_sorted
#psmearpos_clindxpublic <- pnosmear_publicclindx_samples*pTB_clindx_sorted*psmearpos_TBnosmearclindx_sorted

notif_smearpos <- notif_labconf - notif_xpert_smearneg + notif_private_clindx_smearpos +
  notif_public_clindx_smearpos
notif_TB <- notif_labconf + notif_clindx*pTB_clindx_samples

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
set.seed(NULL)

#combine add historical cohort targets - note - mort_samples is used separately
targets_all <- c("prop_m_all"=targets_prev %>% filter(param=="prop_15_smear_adj") %>% pull(philippines), #proportion of active TB that are smear+
             "prop_s_all"=targets_prev %>% filter(param=="prop_15_symptom_adj") %>% pull(philippines), #proportion of active TB that are symptom+
             "prop_ms"=targets_prev %>% filter(param=="prop_15_smear_symptom_adj") %>% pull(philippines), #proportion of active TB that are smear+ & symptom+
             "pnr_m_all"=(targets_prev %>% filter(param=="pnr_smear_15_adj") %>% pull(philippines))/notif_adj, #smear+ prevalence to smear+ notifications
             "deaths_tb"=deaths_untx_per_case, #TB non-treated deaths per non-treated cases
             "tb_ms_dead_5yr"=0.5763, #percent symptomatic smear+ dead after 5 years (historical cohort) - from stata output
             "tb_s_dead_5yr"=0.1272, #percent symptomatic smear- dead after 5 years (historical cohort) - from stata output
             "tb_ms_dead_10yr"=0.7096, #percent symptomatic smear+ dead after 10 years (historical cohort) - from stata output
             "tb_s_dead_10yr"=0.2105, #percent symptomatic smear- dead after 10 years (historical cohort) - from stata output
             "prop_m_notif"=0.45 #percent of notifications (that are true TB) that are smear-positive (optional prevalence target)
             )
targets_all_lb <- c("prop_m_all"=targets_prev %>% filter(param=="prop_15_smear_adj_lb") %>% pull(philippines), #proportion of active TB that are smear+
                    "prop_s_all"=targets_prev %>% filter(param=="prop_15_symptom_adj_lb") %>% pull(philippines), #proportion of active TB that are symptom+
                    "prop_ms"=targets_prev %>% filter(param=="prop_15_smear_symptom_adj_lb") %>% pull(philippines), #proportion of active TB that are smear+ & symptom+
                    "pnr_m_all"=(targets_prev %>% filter(param=="pnr_smear_15_adj_lb") %>% pull(philippines))/notif_adj_lb, #smear+ prevalence to smear+ notifications
                    "deaths_tb"=deaths_untx_per_case_lb, #TB non-treated deaths per non-treated cases
                    "tb_ms_dead_5yr"=0.51, #percent symptomatic smear+ dead after 5 years (historical cohort) - from stata output
                    "tb_s_dead_5yr"=0.09, #percent symptomatic smear- dead after 5 years (historical cohort) - from stata output
                    "tb_ms_dead_10yr"=0.65, #percent symptomatic smear+ dead after 10 years (historical cohort) - from stata output
                    "tb_s_dead_10yr"=0.15, #percent symptomatic smear- dead after 10 years (historical cohort) - from stata output
                    "prop_m_notif"=0.32 #percent of notifications (that are true TB) that are smear-positive (optional prevalence target)
)
targets_all_ub <- c("prop_m_all"=targets_prev %>% filter(param=="prop_15_smear_adj_ub") %>% pull(philippines), #proportion of active TB that are smear+
                    "prop_s_all"=targets_prev %>% filter(param=="prop_15_symptom_adj_ub") %>% pull(philippines), #proportion of active TB that are symptom+
                    "prop_ms"=targets_prev %>% filter(param=="prop_15_smear_symptom_adj_ub") %>% pull(philippines), #proportion of active TB that are smear+ & symptom+
                    "pnr_m_all"=(targets_prev %>% filter(param=="pnr_smear_15_adj_ub") %>% pull(philippines))/notif_adj_ub, #smear+ prevalence to smear+ notifications
                    "deaths_tb"=deaths_untx_per_case_ub, #TB non-treated deaths per non-treated cases
                    "tb_ms_dead_5yr"=0.64, #percent symptomatic smear+ dead after 5 years (historical cohort) - from stata output
                    "tb_s_dead_5yr"=0.16, #percent symptomatic smear- dead after 5 years (historical cohort) - from stata output
                    "tb_ms_dead_10yr"=0.77, #percent symptomatic smear+ dead after 10 years (historical cohort) - from stata output
                    "tb_s_dead_10yr"=0.26, #percent symptomatic smear- dead after 10 years (historical cohort) - from stata output
                    "prop_m_notif"=0.62 #percent of notifications (that are true TB) that are smear-positive (optional prevalence target)
)

#calculate quantities used in parameterizing the likelihood for the PNR target - using dampack gamma_params
sd_gamma <- (targets_all_ub[["pnr_m_all"]]-targets_all_lb[["pnr_m_all"]])/3.9 
pnr_gamma_shape <- gamma_params(targets_all[["pnr_m_all"]], sd_gamma)$shape 
pnr_gamma_scale <- gamma_params(targets_all[["pnr_m_all"]], sd_gamma)$scale
pnr_params <- list("pnr_gamma_shape"=pnr_gamma_shape, "pnr_gamma_scale"=pnr_gamma_scale)

#save params and targets to RDA file for faster loading in the future/on MARCC
save(n, verbose, cyc_len, params, names_params_calib, names_params_calib_mult,
     params_calib_prev, params_fixed_prev, params_calib_prev_mult,
     params_calib_hist, params_fixed_hist, params_fixed_hist_mult, 
     params_depend, 
     priors_prev_lb, priors_prev_ub, priors_prev_mult_lb, priors_prev_mult_ub,
     priors_hist_lb, priors_hist_ub, mort_samples, prop_m_notif_smooth,
     targets_all, targets_all_lb, targets_all_ub, pnr_params,
     file="data/params_targets.Rda")
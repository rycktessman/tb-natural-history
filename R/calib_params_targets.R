setwd("~/GitHub/tb-natural-history")

library(MASS) #used to induce correlations in targets
library(tidyverse) #load dplyr after MASS so that "select" works as expected
library(dampack)

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

#define objects needed to run calibration

#define calibration lower and upper bounds
params_lb <- list(
  "p_m"=0, #smear progression (from subclinical)
  "p_s"=0, #symptom progression (from smear-negative)
  "r_m"=0, #smear regression (from subclinical)
  "r_s"=0, #symptom regression (from smear-negative)
  "c_sp"=0, #spontaneous resolution (from subclinical smear-negative)
  "c_tx"=0, #treatment (from symptomatic smear-negative)
  "m_tb"=0, #TB mortality (from symptomatic smear-positive)
  "a_p_m"=1, #RR of smear/symptom progression if symptomatic/smear-pos (vs. subclinical/smear-neg)
  "a_r_m"=0, #RR of smear/symptom regression if symptomatic/smear-pos (vs. subclinical/smear-neg)
  "a_tx"=1, #RR of treatment if symptomatic smear-positive (vs. symptomatic smear-negative)
  "a_m"=1 #RR of TB mortality if symptomatic smear-positive (vs. symptomatic smear-negative)
)

params_ub <- list(
  "p_m"=0.2, #smear progression (from subclinical)
  "p_s"=0.2, #symptom progression (from smear-negative)
  "r_m"=0.2, #smear regression (from subclinical)
  "r_s"=0.5, #symptom regression (from smear-negative)
  "c_sp"=0.5, #spontaneous resolution (from subclinical smear-negative)
  "c_tx"=0.25, #treatment (from symptomatic smear-negative)
  "m_tb"=0.1, #TB mortality (from symptomatic smear-positive)
  "a_p_m"=10, #RR of smear/symptom progression if symptomatic/smear-pos (vs. subclinical/smear-neg)
  "a_r_m"=1, #RR of smear/symptom regression if symptomatic/smear-pos (vs. subclinical/smear-neg)
  "a_tx"=20, #RR of treatment if symptomatic smear-positive (vs. symptomatic smear-negative)
  "a_m"=20 #RR of TB mortality if symptomatic smear-positive (vs. symptomatic smear-negative)
)

#define uncalibrated (fixed) parameters
params_fixed_prev <- list("p_c"=0, #no spontaneous cure back to TB 
                          "m_ac"=0.000418, #all-cause mortality is fixed at present-day rates
                          "inflows"=1
)

params_fixed_hist <- list("p_c"=0, #no spontaneous cure back to TB
                          "m_ac"=0.000837, #all-cause mortality is fixed
                          "inflows"=0, #no inflows in historical cohort version
                          "c_tx"=0 #no treatment in historical cohort version
)

#define calibration targets
#1. % prevalent TB smear-positive 
tb_m <- 0.368
tb_m_lb <- 0.319
tb_m_ub <- 0.419
#2. % prevalence TB symptomatic
tb_s <- 0.290
tb_s_lb <- 0.244
tb_s_ub <- 0.338
#3. % prevalence TB smear-positive and symptomatic
tb_ms <- 0.171
tb_ms_lb <- 0.134
tb_ms_ub <- 0.230
#4. smear-positive prevalence to smear-positive notifications ratio (PNR)
pnr_m_all <- 1.95
pnr_m_all_lb <- 1.05
pnr_m_all_ub <- 3.10
#5. untreated TB deaths per case
#5a. treated CFR - use lognormal distribution given skew
cfr_samples <- rlnorm(n=100000, meanlog=log(0.0243), sdlog=0.5)
#5b. total deaths estimate = overall CFR numerator: deaths - use normal distribution
deaths_samples <- rnorm(n=100000, mean=29000, sd=2041)
#5c. overall CFR denominator: notifications - use normal distribution
notif_samples <- rgamma(n=100000, shape=15.26, scale=37615)
#5d. all ages prevalence
prev_samples <- rnorm(n=100000, mean=8.05, sd=0.505)
#5e. calculate numbers treated that didn't fail/LTFU and numbers treated that did fail/LTRU
cases_failLTFU <- 14585
cases_tx_nofailLTFU <- 332308 - cases_failLTFU
#5f. combine all components to calculate target
deaths_untx_samples <- deaths_samples - 
  (1/0.75)*(cfr_samples*cases_tx_nofailLTFU + 
               (deaths_samples/notif_samples)*cases_failLTFU)
pop <- 103664
deaths_untx_per_case_samples <- deaths_untx_samples/(prev_samples*pop)
#5g. convert to vector of empirical probabilities
deaths_untx_per_case <- mean(deaths_untx_per_case_samples)
deaths_untx_per_case_lb <- quantile(deaths_untx_per_case_samples, 0.025)[[1]]
deaths_untx_per_case_ub <- quantile(deaths_untx_per_case_samples, 0.975)[[1]]
deaths_untx_per_case_samples <- deaths_untx_per_case_samples[deaths_untx_per_case_samples>=0]
mort_samples_tmp <- data.frame(table(round(deaths_untx_per_case_samples*1000))/length(deaths_untx_per_case_samples))
mort_samples <- mort_samples_tmp$Freq
names(mort_samples) <- mort_samples_tmp$Var1
#6. % notifications smear-positive
notif_labconf <- 115200 
notif_xpert <- 41759 
notif_clindx <- 200572
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
#smooth to generate empirical probability distribution
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

#combine add historical cohort targets - note - mort_samples is used separately
targets_all <- c("prop_m_all"=tb_m,
                 "prop_s_all"=tb_s,
                 "prop_ms"=tb_ms,
                 "pnr_m_all"=pnr_m_all,
                 "deaths_tb"=deaths_untx_per_case, #TB non-treated deaths per non-treated cases
                 "tb_ms_dead_5yr"=0.5763, #percent symptomatic smear+ dead after 5 years (historical cohort) - from stata output
                 "tb_s_dead_5yr"=0.1272, #percent symptomatic smear- dead after 5 years (historical cohort) - from stata output
                 "tb_ms_dead_10yr"=0.7096, #percent symptomatic smear+ dead after 10 years (historical cohort) - from stata output
                 "tb_s_dead_10yr"=0.2105, #percent symptomatic smear- dead after 10 years (historical cohort) - from stata output
                 "prop_m_notif"=0.45 #percent of notifications (that are true TB) that are smear-positive (optional prevalence target)
)
targets_all_lb <- c("prop_m_all"=tb_m_lb,
                    "prop_s_all"=tb_s_lb, 
                    "prop_ms"=tb_ms_lb,
                    "pnr_m_all"=pnr_m_all_lb,
                    "deaths_tb"=deaths_untx_per_case_lb, #TB non-treated deaths per non-treated cases
                    "tb_ms_dead_5yr"=0.51, #percent symptomatic smear+ dead after 5 years (historical cohort) - from stata output
                    "tb_s_dead_5yr"=0.09, #percent symptomatic smear- dead after 5 years (historical cohort) - from stata output
                    "tb_ms_dead_10yr"=0.65, #percent symptomatic smear+ dead after 10 years (historical cohort) - from stata output
                    "tb_s_dead_10yr"=0.15, #percent symptomatic smear- dead after 10 years (historical cohort) - from stata output
                    "prop_m_notif"=0.32 #percent of notifications (that are true TB) that are smear-positive (optional prevalence target)
)
targets_all_ub <- c("prop_m_all"=tb_m_ub,
                    "prop_s_all"=tb_s_ub,
                    "prop_ms"=tb_ms_ub,
                    "pnr_m_all"=pnr_m_all_ub,
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

save(params_fixed_prev, params_fixed_hist, params_lb, params_ub, 
     targets_all, targets_all_lb, targets_all_ub, 
     mort_samples, prop_m_notif_smooth, pnr_params,
     file="calibration_inputs.Rda")
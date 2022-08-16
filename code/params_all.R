#generate RDA file with parameter info

setwd("~/GitHub/tb-natural-history")
library(MASS)
library(tidyverse)
library(dampack)
library(readxl)

cyc_len <- 1/52 #weekly timestep

#PARAMETERS: #progression/regression probabilities: p=progress, r=regress, m=smear microscopy, s=symptoms, a's are multipliers

#specify values for non-calibrated parameters
p_c <- 0.0 #probability of progression from spontaneous cure to smear- symptom-
#mortality probabilities: present-day vary by country
m_ac_bangladesh_annual <- 0.00389 #all-cause/non-TB mortality probability (applied to symptom- TB too - annual avg. mortality probability among 15-64 years in Bangladesh)
m_ac_cambodia_annual <- 0.00655 #all-cause/non-TB mortality probability (applied to symptom- TB too - annual avg. mortality probability among 15-64 years in Cambodia)
m_ac_nepal_annual <- 0.00483 #all-cause/non-TB mortality probability (applied to symptom- TB too - annual avg. mortality probability among 15-64 years in Nepal)
m_ac_philippines_annual <- 0.005 #all-cause/non-TB mortality probability (applied to symptom- TB too - annual avg. mortality probability among 15-64 years in Philippines)
m_ac_vietnam_annual <- 0.004 #all-cause/non-TB mortality probability (applied to symptom- TB too - annual avg. mortality probability among 15-64 years in Vietnam)
m_ac_present_annual <- c("Bangladesh"=m_ac_bangladesh_annual,
                  "Cambodia"=m_ac_cambodia_annual,
                  "Nepal"=m_ac_nepal_annual,
                  "Philippines"=m_ac_philippines_annual,
                  "Vietnam"=m_ac_vietnam_annual)
m_ac_hist_annual <- 0.01 #historical all-cause non-TB mortality probability in early 20th century Europe (from Ragonnet)

#convert annual background mortality probabilities to weekly
m_ac_present <- 1-exp(-1*-log(1-m_ac_present_annual)*cyc_len)
m_ac_hist <- 1-exp(-1*-log(1-m_ac_hist_annual)*cyc_len)

#starting values for calibrated params (just for testing)
p_m <- 0.04 #probability of progression from smear- symptom- to smear+ symptom-
p_s <- 0.05 #probability of progression from smear- symptom- to smear- symptom+
a_p_m <- 2 #multiplier on probability of progression from smear- symptom+ to smear+ symptom+ (vs. smear- symptom- to smear+ symptom-)
a_p_s <- a_p_m #multiplier on probability of progression from smear+ symptom- to smear+ symptom+ (vs. smear- symptom- to smear- symptom+)
r_m <- 0.02 #probability of regression from smear+ symptom- to smear- symptom-
r_s <- 0.04 #probability of regression from smear- symptom+ to smear- symptom-
a_r_m <- 0.3 #multiplier on probability of regression from smear+ symptom+ to smear- symptom+ (vs. smear+ symptom- to smear- symptom-)
a_r_s <- a_r_m #multiplier on probability of regression from smear+ symptom+ to smear+ symptom- (vs. smear- symptom+ to smear- symptom-)
#cure/diagnosis/treatment probabilities: c=cure, a's are multipliers
c_sp <- 0.005 #probability of spontaneous cure (from smear- symptom-)
c_tx <- 0.15 #probability of diagnosis and treatment (from smear- symptom+)
a_tx <- 2 #multiplier on probability of diagnosis and treatment if smear+ symptom+ (vs. smear- symptom+)
#TB mortality
m_tb <- 0.1/12 #mortality probability from smear- symptom+ TB 
a_m <- a_tx #multiplier on mortality from symptom+ smear+TB (vs. symptom+ smear-)

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
  "m_ac"=m_ac_hist,
  "m_tb"=m_tb,
  "a_m"=a_m,
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
                        "a_m"="Mortality Multiplier",
                        "m_tb"="TB Mortality",
                        "a_tx"="Treatment Multiplier")

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
  "m_tb"=m_tb,
  "a_tx"=a_tx
)

#parameters that stay fixed (doesn't include dependent params or country-specific mortality) - calibration to prevalence data
params_fixed_prev <- list (
  "p_c"=p_c, #no spontaneous cure back to TB -  fixed except in sensitivity analysis
  "inflows"=1 #inflows are turned on in present-day calib, off in historical
)

#parameters that are being calibrated (dependent params removed) - calibration to historical cohort data (starting values for optim)
params_calib_hist <- params_calib_prev[names(params_calib_prev)!="c_tx"]

#parameters that stay fixed (doesn't include dependent params) - calibration to historical cohort data
params_fixed_hist <- list(
  "p_c"=p_c, #no spontaneous cure back to TB -  fixed except in sensitivity analysis
  "m_ac"=m_ac_hist, #all-cause mortality is fixed
  "inflows"=0, #inflows are turned on in present-day calib, off in historical
  "c_tx"=0, #no treatment in historical cohort version
  "a_tx"=0 #no multiplier on tx in historical version
)

#parameters that depend on other parameters (need to handle separately in calibration)
#same for both prevalence and cohort versions of calibration
params_depend <- list(
  "a_p_s"=a_p_s, #equals a_p_m except in sensitivity analysis
  "a_r_s"=a_r_s #equals a_r_m except in sensitivity analysis
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
  "m_tb"=m_tb_lb,
  "a_tx"=a_tx_lb
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
  "m_tb"=m_tb_ub,
  "a_tx"=a_tx_ub
)

priors_hist_lb <- priors_prev_lb[names(priors_prev_lb)!="c_tx" & names(priors_prev_lb)!="a_tx"]
priors_hist_ub <- priors_prev_ub[names(priors_prev_ub)!="c_tx" & names(priors_prev_ub)!="a_tx"]

#save to RDA file
save(m_ac_present, params, names_params_calib, 
     params_calib_prev, params_fixed_prev,
     params_calib_hist, params_fixed_hist,
     params_depend, 
     priors_prev_lb, priors_prev_ub, 
     priors_hist_lb, priors_hist_ub, 
     file="data/params_all.Rda")

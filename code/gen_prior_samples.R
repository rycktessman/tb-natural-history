#script to generate prior distributions to use in calibration/upload to MARCC

library(lhs)
library(dplyr)
library(mvtnorm)
library(matrixStats)
source("code/calib_functions.R")
source("code/model_functions.R")
load("data/params_all.Rda")
path_out <- "output/main/"
cyc_len <- 1/12 #weekly or monthly timestep

country <- "Philippines" #doesn't matter except for likelihoods, which we don't use
load(paste0("data/targets_", tolower(country), ".Rda"))

priors_prev_lb[["a_p_s"]] <- priors_prev_lb[["a_p_m"]]
priors_prev_ub[["a_p_s"]] <- priors_prev_ub[["a_p_m"]]
priors_prev_lb[["a_r_s"]] <- priors_prev_lb[["a_r_m"]]
priors_prev_ub[["a_r_s"]] <- priors_prev_ub[["a_r_m"]]
params_fixed_prev[["m_ac"]] <- m_ac_present[[country]]
params_calib_prev <- c(params_calib_prev, 
                       "a_p_s"=params_calib_prev[["a_p_m"]],
                       "a_r_s"=params_calib_prev[["a_r_m"]])
params_calib_hist <- c(params_calib_hist, 
                       "a_p_s"=params_calib_hist[["a_p_m"]],
                       "a_r_s"=params_calib_hist[["a_r_m"]])
RR_free <- 1
smear_hist_calib <- 0
flag_symptom_dur <- 0
smear_notif_override <- NA
RR_regress_recip <- 1

if(RR_regress_recip==1) {
  priors_prev_lb[["a_r_m_recip"]] <- 1
  priors_prev_ub[["a_r_m_recip"]] <- 20
  priors_hist_lb[["a_r_m_recip"]] <- 1
  priors_hist_ub[["a_r_m_recip"]] <- 20
  params_calib_prev[["a_r_m_recip"]] <- 2
  params_calib_hist[["a_r_m_recip"]] <- 2
  priors_prev_lb[["a_r_m"]] <- NULL
  priors_prev_ub[["a_r_m"]] <- NULL
  priors_hist_lb[["a_r_m"]] <- NULL
  priors_hist_ub[["a_r_m"]] <- NULL
  params_calib_prev[["a_r_m"]] <- NULL
  params_calib_hist[["a_r_m"]] <- NULL
  if(RR_free==1) {
    priors_prev_lb[["a_r_s_recip"]] <- 1
    priors_prev_ub[["a_r_s_recip"]] <- 20
    priors_hist_lb[["a_r_s_recip"]] <- 1
    priors_hist_ub[["a_r_s_recip"]] <- 20
    params_calib_prev[["a_r_s_recip"]] <- 2
    params_calib_hist[["a_r_s_recip"]] <- 2
    priors_prev_lb[["a_r_s"]] <- NULL
    priors_prev_ub[["a_r_s"]] <- NULL
    priors_hist_lb[["a_r_s"]] <- NULL
    priors_hist_ub[["a_r_s"]] <- NULL
    params_calib_prev[["a_r_s"]] <- NULL
    params_calib_hist[["a_r_s"]] <- NULL
  }
}

#sample parameter sets
priors <- sample.prior(100000)
out_like <- output_like(priors)
out <- bind_cols(out_like$out_prev, out_like$out_hist_pos,
                 out_like$out_hist_neg)
out_priors <- bind_cols(as.data.frame(priors), out)
write.csv(out_priors, file=paste0(path_out, "out_priors_combined_recipprior.csv"))

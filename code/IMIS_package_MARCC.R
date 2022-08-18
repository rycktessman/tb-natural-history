library(lhs)
library(dplyr)
library(mvtnorm)
library(matrixStats)

source("code/model_functions.R")
source("code/calib_functions.R")

chain <- as.character(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#calibration options
RR_free <- 0 #4 free RR parameters in this version
spont_progress <- 0 #whether those who have spontaneously resolved can progress back to smear- symptom- TB
spont_prog <- 0.15 #what probability to use if spont_progress is 1
smear_hist_calib <- 0 #whether to include historical targets on bacillary status over time
no_10yr_hist <- 0 #whether to include 10 year historical survival as calibration targets
flag_symptom_dur <- 0 #don't exclude parameter sets where less than 80-90% spend >2 weeks symptomatic
country <- "Philippines"
cyc_len <- 1/12 #weekly or monthly timestep

#load files and implement options
load("data/params_all.Rda")
load(paste0("data/targets_", tolower(country), ".Rda"))
path_out <- paste0("output/", tolower(country))
if(RR_free==1) {
  priors_prev_lb[["a_p_s"]] <- priors_prev_lb[["a_p_m"]]
  priors_prev_ub[["a_p_s"]] <- priors_prev_ub[["a_p_m"]]
  priors_prev_lb[["a_r_s"]] <- priors_prev_lb[["a_r_m"]]
  priors_prev_ub[["a_r_s"]] <- priors_prev_ub[["a_r_m"]]
  params_calib_prev <- c(params_calib_prev, 
                         "a_p_s"=params_calib_prev[["a_p_m"]],
                         "a_r_s"=params_calib_prev[["a_r_m"]])
  params_calib_hist <- c(params_calib_hist, 
                         "a_p_s"=params_calib_hist[["a_p_m"]],
                         "a_r_s"=params_calib_hist[["a_r_m"]])
  path_out <- paste0(path_out, "_rrfree")
}
if(spont_progress==1) {
  params_fixed_prev[["p_c"]] <- 1-exp(log(1-spont_prog)*cyc_len) #weekly probability corresponding to annual probability of 15%
  params_fixed_hist[["p_c"]] <- params_fixed_prev[["p_c"]]
  path_out <- paste0(path_out, "_spontprog")
}
if(smear_hist_calib==1) {
  #switch to using calc_like_smear_hist instead of calc_like
  calc_like <- function(out, tr, tr_lb, tr_ub, mort_samples, prev_cases, 
                        prop_m_notif_smooth, pnr_params, calib_type, country) {
    like <- calc_like_smear_hist(out, tr, tr_lb, tr_ub, mort_samples, prev_cases,
                                 prop_m_notif_smooth, pnr_params, calib_type, country)
    return(like)
  }
  #sinding-larsen
  targets_all[["tb_smear_4yr_1"]] <- 61
  targets_all[["alive_4yr_1"]] <- 366
  #braeuning neisen
  targets_all[["tb_smear_4yr_2"]] <- 84
  targets_all[["alive_4yr_2"]] <- 166
  #griep
  targets_all[["tb_smear_4yr_3"]] <- 14
  targets_all[["alive_4yr_3"]] <- 57
  path_out <- paste0(path_out, "_smearhist")
}
if(no_10yr_hist==1) {
  #version without the 10-year mortality targets
  calc_like <- function(out, tr, tr_lb, tr_ub, mort_samples, prev_cases,
                        prop_m_notif_smooth, pnr_params, calib_type, country) { #outputs, targets, and upper/lower confidence bounds on targets
    like <- calc_like_no10(out, tr, tr_lb, tr_ub, mort_samples, prev_cases,
                           prop_m_notif_smooth, pnr_params, calib_type, country)
    return(like)
  }
  path_out <- paste0(path_out, "_no10")
}
if(RR_free==0 & spont_progress==0 & smear_hist_calib==0 & no_10yr_hist==0) {
  path_out <- paste0(path_out, "_base")
}
path_out <- paste0(path_out, "/")
params_fixed_prev[["m_ac"]] <- m_ac_present[[country]]

#define functions and arguments needed for IMIS package
B <- 10000 #10,000 samples per IMIS round (100,000 initially)
B.re <- 1000 #1000 samples of the posterior
number_k <- 12 #run 11 rounds of IMIS
D <- 0 #don't optimize first

#B <- 1000 #versions for testing
#B.re <- 1000
#number_k <- 5
#copy of IMIS package function by Raftery & Bao with additional outputs saved
IMIS_copy <- function(B, B.re, number_k, D) {
  B0 = B * 10
  X_all = X_k = sample.prior(B0)
  if (is.vector(X_all)) 
    Sig2_global = var(X_all)
  if (is.matrix(X_all)) 
    Sig2_global = cov(X_all)
  stat_all = matrix(NA, 6, number_k)
  center_all = prior_all = out_all = like_all = NULL
  sigma_all = list()
  if (D >= 1) 
    option.opt = 1
  if (D == 0) {
    option.opt = 0
    D = 1
  }
  for (k in 1:number_k) {
    ptm.like = proc.time()
    prior_all = c(prior_all, prior(X_k))
    out_like_all = output_like(X_k)
    out_all = rbind(out_all, cbind(out_like_all$out_prev, 
                                           out_like_all$out_hist_pos, 
                                           out_like_all$out_hist_neg, 
                                           "like_prev"=out_like_all$like_prev,
                                           "like_hist_pos"=out_like_all$like_hist_pos, 
                                           "like_hist_neg"=out_like_all$like_hist_neg,
                                           "log_like"=out_like_all$log_like, 
                                           "like"=out_like_all$like))
    like_all = c(like_all, out_like_all$like)
    ptm.use = (proc.time() - ptm.like)[3]
    if (k == 1) 
      print(paste(B0, "likelihoods are evaluated in", 
                  round(ptm.use/60, 2), "minutes"))
    if (k == 1) 
      envelop_all = prior_all
    if (k > 1) 
      envelop_all = apply(rbind(prior_all * B0/B, gaussian_all), 
                          2, sum)/(B0/B + D + (k - 2))
    Weights = prior_all * like_all/envelop_all
    stat_all[1, k] = log(mean(Weights))
    Weights = Weights/sum(Weights)
    stat_all[2, k] = sum(1 - (1 - Weights)^B.re)
    stat_all[3, k] = max(Weights)
    stat_all[4, k] = 1/sum(Weights^2)
    stat_all[5, k] = -sum(Weights * log(Weights), na.rm = TRUE)/log(length(Weights))
    stat_all[6, k] = var(Weights/mean(Weights))
    if (k == 1) 
      print("Stage   MargLike   UniquePoint   MaxWeight   ESS")
    print(c(k, round(stat_all[1:4, k], 3)))
    if (k == 1 & option.opt == 1) {
      if (is.matrix(X_all)) 
        Sig2_global = cov(X_all[which(like_all > min(like_all)), 
                                ])
      X_k = which_exclude = NULL
      label_weight = sort(Weights, decreasing = TRUE, index = TRUE)
      which_remain = which(Weights > label_weight$x[B0])
      size_remain = length(which_remain)
      for (i in 1:D) {
        important = NULL
        if (length(which_remain) > 0) 
          important = which_remain[which(Weights[which_remain] == 
                                           max(Weights[which_remain]))]
        if (length(important) > 1) 
          important = sample(important, 1)
        if (is.vector(X_all)) 
          X_imp = X_all[important]
        if (is.matrix(X_all)) 
          X_imp = X_all[important, ]
        which_exclude = union(which_exclude, important)
        which_remain = setdiff(which_remain, which_exclude)
        posterior = function(theta) {
          -log(prior(theta)) - log(likelihood(theta))
        }
        if (is.vector(X_all)) {
          if (length(important) == 0) 
            X_imp = center_all[1]
          optimizer = optim(X_imp, posterior, method = "BFGS", 
                            hessian = TRUE, control = list(parscale = sqrt(Sig2_global)/10, 
                                                           maxit = 5000))
          print(paste("maximum posterior=", round(-optimizer$value, 
                                                  2), ", likelihood=", round(log(likelihood(optimizer$par)), 
                                                                             2), ", prior=", round(log(prior(optimizer$par)), 
                                                                                                   2), ", time used=", round(ptm.use/60, 
                                                                                                                             2), "minutes, convergence=", optimizer$convergence))
          center_all = c(center_all, optimizer$par)
          sigma_all[[i]] = solve(optimizer$hessian)
          X_k = c(X_k, rnorm(B, optimizer$par, sqrt(sigma_all[[i]])))
          distance_remain = abs(X_all[which_remain] - 
                                  optimizer$par)
        }
        if (is.matrix(X_all)) {
          if (length(important) == 0) 
            X_imp = center_all[1, ]
          ptm.opt = proc.time()
          optimizer = optim(X_imp, posterior, method = "Nelder-Mead", 
                            control = list(maxit = 1000, parscale = sqrt(diag(Sig2_global))))
          theta.NM = optimizer$par
          optimizer = optim(theta.NM, posterior, method = "BFGS", 
                            hessian = TRUE, control = list(parscale = sqrt(diag(Sig2_global)), 
                                                           maxit = 1000))
          ptm.use = (proc.time() - ptm.opt)[3]
          print(paste("maximum posterior=", round(-optimizer$value, 
                                                  2), ", likelihood=", round(log(likelihood(optimizer$par)), 
                                                                             2), ", prior=", round(log(prior(optimizer$par)), 
                                                                                                   2), ", time used=", round(ptm.use/60, 
                                                                                                                             2), "minutes, convergence=", optimizer$convergence))
          center_all = rbind(center_all, optimizer$par)
          if (min(eigen(optimizer$hessian)$values) > 
              0) 
            sigma_all[[i]] = solve(optimizer$hessian)
          if (min(eigen(optimizer$hessian)$values) <= 
              0) {
            eigen.values = eigen(optimizer$hessian)$values
            eigen.values[which(eigen.values < 0)] = 0
            hessian = eigen(optimizer$hessian)$vectors %*% 
              diag(eigen.values) %*% t(eigen(optimizer$hessian)$vectors)
            sigma_all[[i]] = solve(hessian + diag(1/diag(Sig2_global)))
          }
          X_k = rbind(X_k, rmvnorm(B, optimizer$par, 
                                   sigma_all[[i]]))
          distance_remain = mahalanobis(X_all[which_remain, 
                                              ], optimizer$par, diag(diag(Sig2_global)))
        }
        label_dist = sort(distance_remain, decreasing = FALSE, 
                          index = TRUE)
        which_exclude = union(which_exclude, which_remain[label_dist$ix[1:floor(size_remain/D)]])
        which_remain = setdiff(which_remain, which_exclude)
      }
      if (is.matrix(X_all)) 
        X_all = rbind(X_all, X_k)
      if (is.vector(X_all)) 
        X_all = c(X_all, X_k)
    }
    if (k > 1 | option.opt == 0) {
      important = which(Weights == max(Weights))
      if (length(important) > 1) 
        important = important[1]
      if (is.matrix(X_all)) 
        X_imp = X_all[important, ]
      if (is.vector(X_all)) 
        X_imp = X_all[important]
      if (is.matrix(X_all)) 
        center_all = rbind(center_all, X_imp)
      if (is.vector(X_all)) 
        center_all = c(center_all, X_imp)
      if (is.matrix(X_all)) 
        distance_all = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)))
      if (is.vector(X_all)) 
        distance_all = abs(X_all - X_imp)
      label_nr = sort(distance_all, decreasing = FALSE, 
                      index = TRUE)
      which_var = label_nr$ix[1:B]
      if (is.matrix(X_all)) 
        Sig2 = cov.wt(X_all[which_var, ], wt = Weights[which_var] + 
                        1/length(Weights), cor = FALSE, center = X_imp, 
                      method = "unbias")$cov
      if (is.vector(X_all)) {
        Weights_var = Weights[which_var] + 1/length(X_all)
        Weights_var = Weights_var/sum(Weights_var)
        Sig2 = (X_all[which_var] - X_imp)^2 %*% Weights_var
      }
      sigma_all[[D + k - 1]] = Sig2
      if (is.matrix(X_all)) 
        X_k = rmvnorm(B, X_imp, Sig2)
      if (is.vector(X_all)) 
        X_k = rnorm(B, X_imp, sqrt(Sig2))
      if (is.matrix(X_all)) 
        X_all = rbind(X_all, X_k)
      if (is.vector(X_all)) 
        X_all = c(X_all, X_k)
    }
    if (k == 1) {
      gaussian_all = matrix(NA, D, B0 + D * B)
      for (i in 1:D) {
        if (is.matrix(X_all)) 
          gaussian_all[i, ] = dmvnorm(X_all, center_all[i, 
                                                        ], sigma_all[[i]])
        if (is.vector(X_all)) 
          gaussian_all[i, ] = dnorm(X_all, center_all[i], 
                                    sqrt(sigma_all[[i]]))
      }
    }
    if (k > 1) {
      if (is.vector(X_all)) 
        gaussian_new = matrix(0, D + k - 1, length(X_all))
      if (is.matrix(X_all)) 
        gaussian_new = matrix(0, D + k - 1, dim(X_all)[1])
      if (is.matrix(X_all)) {
        gaussian_new[1:(D + k - 2), 1:(dim(X_all)[1] - 
                                         B)] = gaussian_all
        gaussian_new[D + k - 1, ] = dmvnorm(X_all, X_imp, 
                                            sigma_all[[D + k - 1]])
        for (j in 1:(D + k - 2)) gaussian_new[j, (dim(X_all)[1] - 
                                                    B + 1):dim(X_all)[1]] = dmvnorm(X_k, center_all[j, 
                                                                                                    ], sigma_all[[j]])
      }
      if (is.vector(X_all)) {
        gaussian_new[1:(D + k - 2), 1:(length(X_all) - 
                                         B)] = gaussian_all
        gaussian_new[D + k - 1, ] = dnorm(X_all, X_imp, 
                                          sqrt(sigma_all[[D + k - 1]]))
        for (j in 1:(D + k - 2)) gaussian_new[j, (length(X_all) - 
                                                    B + 1):length(X_all)] = dnorm(X_k, center_all[j], 
                                                                                  sqrt(sigma_all[[j]]))
      }
      gaussian_all = gaussian_new
    }
    if (stat_all[2, k] > (1 - exp(-1)) * B.re) 
      break
  }
  nonzero = which(Weights > 0)
  which_X = sample(nonzero, B.re, replace = TRUE, prob = Weights[nonzero])
  if (is.matrix(X_all)) {
    resample_X = X_all[which_X, ]
    resample_out = out_all[which_X,]
  }
    
  if (is.vector(X_all)) {
    resample_X = X_all[which_X]
    resample_out = out_all[which_X]
  }
  return(list(stat = t(stat_all), resample = resample_X, center = center_all,
              out = resample_out, params_all=X_all, out_all=out_all, weights=Weights))
}

post <- IMIS_copy(B, B.re, number_k, D)

post_params <- bind_cols(data.frame(post$resample), post$out)
write.csv(post_params, file=paste0(path_out, "out_IMIS_", chain, ".csv"), row.names=F)

params_opt <- post$center
write.csv(params_opt, file=paste0(path_out, "IMIS_opt_each_", chain, ".csv"), row.names=F)

stats_out <- post$stat
write.csv(stats_out, file=paste0(path_out, "IMIS_stats_", chain, ".csv"), row.names=F)

out_all <- bind_cols(post$params_all[1:(B*10+(number_k-1)*B),], post$out_all, "weight"=post$weights)
save(out_all, file=paste0(path_out, "IMIS_out_all_", chain, ".RDA"))

weights <- post$weights
save(weights, file=paste0(path_out, "IMIS_wts_all_", chain, ".RDA"))
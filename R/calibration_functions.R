#define functions needed to run calibration

#"IMIS_copy" is a version of IMIS (package by Raftery & Bao: https://cran.r-project.org/package=IMIS)
#it returns additional output compared to the original "IMIS" function 
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
              out = resample_out))
}
#"IMIS_copy" requires "sample.prior", "prior", and "output_like" functions

#"sample.prior" samples from prior distributions
sample.prior <- function(n_samples) {
  i <- 0
  priors_all <- data.frame()
  while(i<n_samples) { #repeat until we have enough parameter samples
    #sample uniforms using LHS
    uniforms <- data.frame(randomLHS(n=n_samples, k=length(params_lb))) #sample using LHS
    names(uniforms) <- names(params_lb)
    #convert uniform[0,1] LHS samples to prior distribution samples
    priors <- uniforms %>% mutate(p_m=qunif(p_m, min=params_lb[["p_m"]], max=params_ub[["p_m"]]),
                                  p_s=qunif(p_s, min=params_lb[["p_s"]], max=params_ub[["p_s"]]),
                                  a_p_m=qunif(a_p_m, min=params_lb[["a_p_m"]], max=params_ub[["a_p_m"]]),
                                  r_m=qunif(r_m, min=params_lb[["r_m"]], max=params_ub[["r_m"]]),
                                  r_s=qunif(r_s, min=params_lb[["r_s"]], max=params_ub[["r_s"]]),
                                  a_r_m=qunif(a_r_m, min=params_lb[["a_r_m"]], max=params_ub[["a_r_m"]]),
                                  c_sp=qunif(c_sp, min=params_lb[["c_sp"]], max=params_ub[["c_sp"]]),
                                  c_tx=qunif(c_tx, min=params_lb[["c_tx"]], max=params_ub[["c_tx"]]),
                                  a_m=qunif(a_m, min=params_lb[["a_m"]], max=params_ub[["a_m"]]),
                                  m_tb=qunif(m_tb, min=params_lb[["m_tb"]], max=params_ub[["m_tb"]]),
                                  a_tx=qunif(a_tx, min=params_lb[["a_tx"]], max=params_ub[["a_tx"]])
    )
    #remove samples when probabilities sum to > 1
    priors <- priors %>% mutate(a_p_s=a_p_m, a_r_s=a_r_m, #dependencies
                                m_ac=params_fixed_hist$m_ac, #use hist bc its larger
                                p_c=params_fixed_hist$p_c) #prev and hist are same 
    priors <- priors %>% 
      mutate(flag1=1*((p_m + p_s + m_ac + c_sp) > 1),
             flag2=1*((r_m + a_p_s*p_s + m_ac) > 1),
             flag3=1*((r_s + a_p_m*p_m + m_tb + m_ac + c_tx) > 1),
             flag4=1*((a_r_s*r_s + a_r_m*r_m + a_m*m_tb + m_ac + a_tx*c_tx) > 1),
             flag5=1*((p_c + m_ac) > 1)
      )
    priors <- priors %>% select(-c(a_p_s, a_r_s, m_ac, p_c)) #remove dependencies
    priors <- priors %>% mutate(flag_sum=flag1+flag2+flag3+flag4+flag5)
    priors <- priors %>% filter(flag_sum==0) %>% select(-starts_with("flag"))
    priors_all <- bind_rows(priors, priors_all)
    i <- nrow(priors_all)
  }
  priors_all <- priors_all[1:n_samples,]
  priors_all <- as.matrix(priors_all) #IMIS package requires matrix or vector
  return(priors_all)
}

#"prior" calculates prior probabilities
prior <- function(params) {
  params <- data.frame(params)
  like <- sapply(names(params), 
                 function(x) dunif(params[[x]], params_lb[[x]], 
                                   params_ub[[x]]), simplify=F, USE.NAMES=T)
  like <- bind_cols(like)
  params <- params %>% mutate(a_p_s=a_p_m, a_r_s=a_r_m, #dependencies
                              m_ac=params_fixed_hist$m_ac, #use hist bc its larger
                              p_c=params_fixed_hist$p_c) #prev and hist are same 
  params <- params %>% 
    mutate(flag1=1*((p_m + p_s + m_ac + c_sp) > 1),
           flag2=1*((r_m + a_p_s*p_s + m_ac) > 1),
           flag3=1*((r_s + a_p_m*p_m + m_tb + m_ac + c_tx) > 1),
           flag4=1*((a_r_s*r_s + a_r_m*r_m + a_m*m_tb + m_ac + a_tx*c_tx) > 1),
           flag5=1*((p_c + m_ac) > 1)
    )
  flags <- params %>%  select(contains("flag"))
  #getting flagged = likelihood of 0, so change 1s to 0s and 0s to 1s, then take product
  flags <- flags*-1 + 1
  like <- cbind(like, flags)
  like <- rowProds(as.matrix(like))
  return(like)
}

#calib_out, and calc_like are helper functions used in output_like (below)

#"calib_out" runs the model with a set of params and returns model outputs (vs. targets) and penalties
calib_out <- function(params_calib, params_fixed, calib_type, t_end, cyc_len, sim_pop) {  
  p_use <- c(params_calib, params_fixed)
  #parameter dependencies
  params_depend <- c("a_p_s"=p_use[["a_p_m"]], 
                     "a_r_s"=p_use[["a_r_m"]]
                     )
  p_use <- c(p_use, params_depend)
  for(t in 1:t_end) {
    curr_pop <- nat_hist_markov(p_use, sim_pop[t,], t)
    sim_pop <- bind_rows(sim_pop, curr_pop)
  }
  if(calib_type=="prev") {
    #prevalence survey targets are cross-sectional - no need for cycle length adjustment. 
    #notifications/deaths need cycle length adjustment (modeled deaths are cumulative so adjust for that too)
    outputs <- c("prop_m_all"=(sim_pop[t+1,"tb_m"]+sim_pop[t+1, "tb_ms"])/sum(sim_pop[t, tb_states]),
                 "prop_s_all"=(sim_pop[t+1,"tb_s"]+sim_pop[t+1, "tb_ms"])/sum(sim_pop[t, tb_states]),
                 "prop_ms"=sim_pop[t+1,"tb_ms"]/sum(sim_pop[t+1,tb_states]),
                 "pnr_m_all"=(sim_pop[t+1,"tb_m"]+sim_pop[t+1,"tb_ms"])/
                   sum((sim_pop[((t+1)-(1/cyc_len)):(t+1),"tb_ms"])*p_use[["c_tx"]]*p_use[["a_tx"]]), #smear-positive divided by smear-positive that get treated
                 "deaths_tb"=(sim_pop[t+1,"died_tb"]-sim_pop[(t+1)-(1/cyc_len),"died_tb"])/sum(sim_pop[t+1,tb_states]), #incremental tb deaths divided by tb cases
                 "prop_m_notif"=(sim_pop[t+1, "tb_ms"]*p_use[["c_tx"]]*p_use[["a_tx"]])/
                   (sim_pop[t+1, "tb_ms"]*p_use[["c_tx"]]*p_use[["a_tx"]] + 
                      sim_pop[t+1, "tb_s"]*p_use[["c_tx"]]) #proportion of notifications that are smear-positive
    )
  } else if(calib_type=="hist_pos") {
    #adjust timing (5 yrs and 10 yrs) for cycle length
    dead_5yr=sim_pop[5/cyc_len, "died_tb"] + sim_pop[5/cyc_len, "died_nontb"]
    dead_10yr=sim_pop[10/cyc_len, "died_tb"] + sim_pop[10/cyc_len, "died_nontb"]
    outputs <- c("tb_ms_dead_5yr"=dead_5yr, "tb_ms_dead_10yr"=dead_10yr)
    } else if (calib_type=="hist_neg") {
      dead_5yr=sim_pop[5/cyc_len, "died_tb"] + sim_pop[5/cyc_len, "died_nontb"]
      dead_10yr=sim_pop[10/cyc_len, "died_tb"] + sim_pop[10/cyc_len, "died_nontb"]
      outputs <- c("tb_s_dead_5yr"=dead_5yr, "tb_s_dead_10yr"=dead_10yr)
    } 
  #apply penalties if compartment sizes are > 1 or < 0 or transitions out exceed 1 - add as separate variable in outputs
  penalties <- 999999*(sum(curr_pop[2:length(curr_pop)])<0) + 
    999999*(sum(curr_pop[alive_states])>1.02) +
    999999*((p_use$p_m + p_use$p_s + p_use$m_ac + p_use$c_sp) > 1) +
    999999*((p_use$r_m + p_use$a_p_s*p_use$p_s + p_use$m_ac) > 1) +
    999999*((p_use$r_s + p_use$a_p_m*p_use$p_m + p_use$m_tb + p_use$m_ac + p_use$c_tx) > 1) +
    999999*((p_use$a_r_s*p_use$r_s + p_use$a_r_m*p_use$r_m + p_use$a_m*p_use$m_tb + p_use$m_ac + p_use$a_tx*p_use$c_tx) > 1) +
    999999*((p_use$p_c + p_use$m_ac) > 1) +
    0
  if(is.na(penalties)) {
    penalties <- 999999*(is.na(sum(curr_pop))) +
      999999*((p_use$p_m + p_use$p_s + p_use$m_ac + p_use$c_sp) > 1) +
      999999*((p_use$r_m + p_use$a_p_s*p_use$p_s + p_use$m_ac) > 1) +
      999999*((p_use$r_s + p_use$a_p_m*p_use$p_m + p_use$m_tb + p_use$m_ac + p_use$c_tx) > 1) +
      999999*((p_use$a_r_s*p_use$r_s + p_use$a_r_m*p_use$r_m + p_use$a_m*p_use$m_tb + p_use$m_ac + p_use$a_tx*p_use$c_tx) > 1) +
      999999*((p_use$p_c + p_use$m_ac) > 1) +
      0
  }
  out_all <- list("outputs"=outputs, "penalties"=penalties, "sim_pop"=sim_pop)
  return(out_all)
}

#"calc_like" calculates likelihood based on model outputs and targets
calc_like <- function(out, targets_all, mort_samples, 
                      prop_m_notif_smooth, pnr_params, calib_type) {
  if(calib_type=="prev") {
    #prevalence survey targets proportion of infections by smear/symptom status - we have actual sample size
    prop_m_all <- dbinom(round(targets_all[["prop_m_all"]]*356), size=356, prob=out[["prop_m_all"]], log=T)
    prop_s_all <- dbinom(round(targets_all[["prop_s_all"]]*359), size=359, prob=out[["prop_s_all"]], log=T)
    prop_ms <- dbinom(round(targets_all[["prop_ms"]]*356), size=356, prob=out[["prop_ms"]], log=T)
    #prevalence to notification ratio: 0 to inf - gamma fits well (parameters estimated using dampack gamma_params)
    pnr_m_all <- dgamma(out[["pnr_m_all"]], shape=pnr_params$pnr_gamma_shape, scale=pnr_params$pnr_gamma_scale, log=T)
    #use empirical distribution for the TB mortality target
    deaths_tb <- unname(log(mort_samples[as.character(round(out[["deaths_tb"]]*1000))]))
    #make likelihood very very small (and decreasing) if > max of all mort samples (min is 0 so no need to do this on low end)
    deaths_tb[(is.na(deaths_tb)|deaths_tb==-Inf) & !is.na(out[["deaths_tb"]])] <- 
      unlist(log(1/abs(round(out[(is.na(deaths_tb)|deaths_tb==-Inf) & !is.na(out[["deaths_tb"]]), 
                                 "deaths_tb"]*1000)-max(unname(mort_samples)))/1000000)) 
    #proportion of notifications that are smear-positive - use empirical distribution
    prop_m_notif <- unname(log(prop_m_notif_smooth[as.character(round(out[["prop_m_notif"]]*100))]))
    log_like_all <- data.frame("prop_m_all"=prop_m_all,
                               "prop_s_all"=prop_s_all,
                               "prop_ms"=prop_ms,
                               "pnr_m_all"=pnr_m_all,
                               "deaths_tb"=deaths_tb,
                               "prop_m_notif"=prop_m_notif)
  } else if(calib_type=="hist_pos") {
    #historical mortality targets - sizes (n) of binomial distributions established to match CIs from meta-regression
    tb_ms_dead_5yr <- dbinom(round(targets_all[["tb_ms_dead_5yr"]]*220), size=220, prob=out[["tb_ms_dead_5yr"]], log=T)
    tb_ms_dead_10yr <- dbinom(round(targets_all[["tb_ms_dead_10yr"]]*220), size=220, prob=out[["tb_ms_dead_10yr"]], log=T)
    log_like_all <- data.frame("tb_ms_dead_5yr"=tb_ms_dead_5yr,
                               "tb_ms_dead_10yr"=tb_ms_dead_10yr)
  } else if(calib_type=="hist_neg") {
    #historical mortality targets - sizes (n) of binomial distributions established to match CIs from meta-regression
    tb_s_dead_5yr <- dbinom(round(targets_all[["tb_s_dead_5yr"]]*325), size=325, prob=out[["tb_s_dead_5yr"]], log=T)
    tb_s_dead_10yr <- dbinom(round(targets_all[["tb_s_dead_10yr"]]*200), size=200, prob=out[["tb_s_dead_10yr"]], log=T)
    log_like_all <- data.frame("tb_s_dead_5yr"=tb_s_dead_5yr,
                               "tb_s_dead_10yr"=tb_s_dead_10yr)
  } else {
    print("error: incorrect calibration type")
  }
  log_like <- rowSums(log_like_all)
  like <- rowProds(as.matrix(exp(log_like_all)))
  like_out <- list("like"=like, "log_like"=log_like)
  return(like_out)
}

#"output_like" calculates model outputs (to match calibration targets) and likelihoods
output_like <- function(params) {
  if(is.vector(params)) {
    params <- data.frame(t(params))
  } else {
    params <- data.frame(params)
  }
  t_end <- 200/cyc_len
  start_pop <- c("t"=0, 
                 "tb"=1, #smear- symptom- TB
                 "tb_m"=0, #smear+ symptom- TB
                 "tb_s"=0, #smear- symptom+ TB
                 "tb_ms"=0, #smear+ symptom+ TB
                 "sp_cure"=0, #spontaneously cured
                 "tx_cure"=0, #cured via diagnosis and treatment
                 "died_tb"=0, #TB death
                 "died_nontb"=0 #non-TB death
  )
  sim_pop <- data.frame() 
  sim_pop <- bind_rows(sim_pop, start_pop)
  params_use <- params %>% select(names(params_lb))
  out_prev <- lapply(1:nrow(params_use), function(x)
    calib_out(params_use[x,], params_fixed_prev, "prev", 
              t_end, cyc_len, sim_pop)$outputs)
  out_prev <- bind_rows(out_prev)
  like_prev <- calc_like(out_prev, targets_all, mort_samples, 
                         prop_m_notif_smooth, pnr_params, "prev")
  
  #HISTORICAL COHORT SMEAR POSITIVE TARGETS
  t_end <- 10/cyc_len
  start_pop <- c("t"=0, 
                 "tb"=0, #smear- symptom- TB
                 "tb_m"=0, #smear+ symptom- TB
                 "tb_s"=0, #smear- symptom+ TB
                 "tb_ms"=1, #smear+ symptom+ TB
                 "sp_cure"=0, #spontaneously cured
                 "tx_cure"=0, #cured via diagnosis and treatment
                 "died_tb"=0, #TB death
                 "died_nontb"=0, #non-TB death
                 "rel_inf"=0 #relative number of secondary infections generated
  )
  sim_pop <- data.frame() 
  sim_pop <- bind_rows(sim_pop, start_pop)
  params_use <- params %>% select(names(params_lb)) %>% 
    select(-c_tx) #no tx in historical version
  out_hist_pos <- lapply(1:nrow(params_use), function(x)
    calib_out(params_use[x,], params_fixed_hist, "hist_pos", 
              t_end, cyc_len, sim_pop)$outputs)
  out_hist_pos <- bind_rows(out_hist_pos)
  like_hist_pos <- calc_like(out_hist_pos, targets_all, mort_samples, 
                             prop_m_notif_smooth, pnr_params, "hist_pos")
  
  #HISTORICAL COHORT SMEAR NEGATIVE TARGETS
  start_pop <- c("t"=0, 
                 "tb"=0, #smear- symptom- TB
                 "tb_m"=0, #smear+ symptom- TB
                 "tb_s"=1, #smear- symptom+ TB
                 "tb_ms"=0, #smear+ symptom+ TB
                 "sp_cure"=0, #spontaneously cured
                 "tx_cure"=0, #cured via diagnosis and treatment
                 "died_tb"=0, #TB death
                 "died_nontb"=0, #non-TB death
                 "rel_inf"=0 #relative number of secondary infections generated
  )
  sim_pop <- data.frame() 
  sim_pop <- bind_rows(sim_pop, start_pop)
  out_hist_neg <- lapply(1:nrow(params_use), function(x)
    calib_out(params_use[x,], params_fixed_hist, "hist_neg", 
              t_end, cyc_len, sim_pop)$outputs)
  out_hist_neg <- bind_rows(out_hist_neg)
  like_hist_neg <- calc_like(out_hist_neg, targets_all, mort_samples, 
                             prop_m_notif_smooth, pnr_params, "hist_neg")
  
  #combine results of all 3 calibration types
  log_like_all <-  like_prev$log_like + like_hist_pos$log_like + like_hist_neg$log_like
  like <- exp(log_like_all)
  like[is.na(like)] <- 0
  out <- list("out_prev"=out_prev, 
              "out_hist_pos"=out_hist_pos, 
              "out_hist_neg"=out_hist_neg, 
              "like_prev"=like_prev$like, 
              "like_hist_pos"=like_hist_pos$like, 
              "like_hist_neg"=like_hist_neg$like,
              "log_like_all"=log_like_all, 
              "like"=like)
  return(out)
}





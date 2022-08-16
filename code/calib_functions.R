#FUNCTIONS USED IN MODEL CALIBRATION

#convert probabilities between different time steps
convert_prob <- function(prob_old, t_new) { #t_new is length of new timestep relatively to old timestep
  rate_old <- -log(1-prob_old)
  rate_new <- rate_old*t_new
  prob_new <- 1-exp(-rate_new)
  return(prob_new)
}

#pull calibration targets based on type of calibration
pull_targets <- function(calib_type, targets_all, country) {
  if(calib_type=="hist_pos") {
    targets <- targets_all[c("tb_ms_dead_5yr","tb_ms_dead_10yr")]
    names <- c("tb_ms_dead_5yr"="5-Year Mortality (%), Smear+ TB",
               "tb_ms_dead_10yr"="10-Year Mortality (%), Smear+ TB")
  } else if(calib_type=="hist_neg") {
    targets <- targets_all[c("tb_s_dead_5yr","tb_s_dead_10yr")]
    names <- c("tb_s_dead_5yr"="5-Year Mortality (%), Smear- TB",
               "tb_s_dead_10yr"="10-Year Mortality (%), Smear- TB")
  } else if(calib_type=="prev" & country %in% c("Philippines", "Cambodia")) {
    targets <- c(targets_all[1:5], targets_all[10])
    names <- c("prop_m"="% Cases Smear+ & Subclinical",
               "prop_s"="% Cases smear- & Symptomatic",
               "prop_ms"="% Cases Smear+ & Symptomatic",
               "pnr_m_all"="Smear+ Prevalence:Notifications",
               "deaths_tb"="Unreated TB Deaths/1000 cases",
               "prop_m_notif"="% Notifications Smear+")
  } else if(calib_type=="prev" & country %in% c("Vietnam", "Nepal", "Bangladesh")) {
    targets <- c(targets_all[1:4], targets_all[9:10])
    names <- c("prop_m"="% Cases Smear+ & Subclinical",
               "prop_s"="% Cases smear- & Symptomatic",
               "prop_ms"="% Cases Smear+ & Symptomatic",
               "pnr_all"="Prevalence:Notifications",
               "deaths_tb"="Unreated TB Deaths/1000 cases",
               "prop_m_notif"="% Notifications Smear+")
  } else {
    targets <- NA
    print("Error: incorrect calibration type specified")
  }
  return(list(targets, names))
}

#sample from prior distributions using latin hypercube sampling
sample_priors <- function(n_samples, priors_prev_lb, priors_prev_ub, 
                          params_fixed, RR_free, flag_symptom_dur) {
  uniforms <- data.frame(randomLHS(n=n_samples, k=length(priors_prev_lb))) #sample using LHS
  names(uniforms) <- names(priors_prev_lb)
  #convert uniform[0,1] LHS samples to prior distribution samples
  priors <- uniforms %>% mutate(p_m=qunif(p_m, min=priors_prev_lb[["p_m"]], max=priors_prev_ub[["p_m"]]),
                                p_s=qunif(p_s, min=priors_prev_lb[["p_s"]], max=priors_prev_ub[["p_s"]]),
                                a_p_m=qunif(a_p_m, min=priors_prev_lb[["a_p_m"]], max=priors_prev_ub[["a_p_m"]]),
                                r_m=qunif(r_m, min=priors_prev_lb[["r_m"]], max=priors_prev_ub[["r_m"]]),
                                r_s=qunif(r_s, min=priors_prev_lb[["r_s"]], max=priors_prev_ub[["r_s"]]),
                                a_r_m=qunif(a_r_m, min=priors_prev_lb[["a_r_m"]], max=priors_prev_ub[["a_r_m"]]),
                                c_sp=qunif(c_sp, min=priors_prev_lb[["c_sp"]], max=priors_prev_ub[["c_sp"]]),
                                c_tx=qunif(c_tx, min=priors_prev_lb[["c_tx"]], max=priors_prev_ub[["c_tx"]]),
                                a_m=qunif(a_m, min=priors_prev_lb[["a_m"]], max=priors_prev_ub[["a_m"]]),
                                m_tb=qunif(m_tb, min=priors_prev_lb[["m_tb"]], max=priors_prev_ub[["m_tb"]]),
                                a_tx=qunif(a_tx, min=priors_prev_lb[["a_tx"]], max=priors_prev_ub[["a_tx"]])
  )
  if(RR_free==1) {
    priors <- priors %>% mutate(a_p_s=qunif(a_p_s, min=priors_prev_lb[["a_p_s"]], max=priors_prev_ub[["a_p_s"]]),
                                  a_r_s=qunif(a_r_s, min=priors_prev_lb[["a_r_s"]], max=priors_prev_ub[["a_r_s"]])
                                  )
  }
  #remove samples when probabilities sum to > 1
  priors <- apply_flags(priors, params_fixed, RR_free, flag_symptom_dur)
  #note: with current priors, only flag2, flag3, and flag4 are binding constraints
  priors <- priors %>% mutate(flag_sum=flag1+flag2+flag3+flag4+flag5+flag6)
  priors <- priors %>% filter(flag_sum==0) %>% select(-starts_with("flag"))
  #this alters the marginal distributions - reflects prior belief that some of these probabilities will have to be smaller so their sum doesn't exceed 1
  return(priors)
}

#calculate percentage still symptom+ after 2 weeks (as param flag in sampling)
calc_prop_s <- function(samples_wk) {
  #2 timesteps - fixed for now, could add as variable later if needed
  t_end <- 2
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
  #run model
  for(t in 1:t_end) {
    curr_pop <- nat_hist_markov(samples_wk, sim_pop[t,], t)
    sim_pop <- bind_rows(sim_pop, curr_pop)
  }
  #calculate proportion still symptomatic in last week
  prop_s <- (sim_pop[(t_end+1),"tb_s"]+sim_pop[(t_end+1),"tb_ms"])/
    sum(sim_pop[(t_end+1), 2:6]) #treat those who die as censored - may have reached 2 weeks but didn't have opportunity
  return(prop_s)
}

#apply flags if transitions sum to > 1
apply_flags <- function(samples, params_fixed, RR_free, flag_symptom_dur) {
  #add dependent parameters so that constraints can be applied
  if(RR_free==1) {
    #sensitivity analysis allows a_p_s and a_r_s to be free parameters
    samples <- samples %>% mutate(m_ac=params_fixed$m_ac, #use hist bc its larger
                                  p_c=params_fixed$p_c) #prev and hist are same 
  } else {
    #in main analysis, a_r_s=a_r_m and a_p_s=a_p_m
    samples <- samples %>% mutate(a_p_s=a_p_m, a_r_s=a_r_m, 
                                  m_ac=params_fixed$m_ac, #use hist bc its larger
                                  p_c=params_fixed$p_c) #prev and hist are same 
  }
  if(flag_symptom_dur==1) {
    #convert all probabilities to weekly
    samples_wk <- samples
    samples_wk[!(names(samples_wk) %in% 
                   c("a_p_s", "a_p_m", "a_r_s", "a_r_m", "a_tx", "a_m"))] <- 
      convert_prob(samples_wk[!(names(samples_wk) %in% 
                                  c("a_p_s", "a_p_m", "a_r_s", "a_r_m", "a_tx", "a_m"))], 
                   12/52)
    samples_wk <- samples_wk %>% mutate(i=1, i_m=1, i_s=1, i_ms=1, inflows=0)
    #run model for 2 timesteps with everyone in smear- symptom+ to calculate % still symptom+ after 2 weeks
    tic()
    prop_s <- lapply(1:nrow(samples_wk), function(x) calc_prop_s(samples_wk[x,]))
    prop_s <- unlist(prop_s)
    toc()
  } else {
    prop_s <- rep(1, nrow(samples)) #value needed not to trigger the flag
  }
  #calculate flags
  samples <- samples %>% 
    mutate(flag1=1*((p_m + p_s + m_ac + c_sp) > 1),
           flag2=1*((r_m + a_p_s*p_s + m_ac) > 1),
           flag3=1*((r_s + a_p_m*p_m + m_tb + m_ac + c_tx) > 1),
           flag4=1*((a_r_s*r_s + a_r_m*r_m + a_m*m_tb + m_ac + a_tx*c_tx) > 1),
           flag5=1*((p_c + m_ac) > 1),
           flag6=1*(prop_s < 0.9)
           #flag6=1*(((1-convert_prob(r_s, 12/52))^2) < 0.9) #90% of ppl should spend 2+ weeks symptomatic
           )
  #remove dependent parameters
  if(RR_free==1) {
    samples <- samples %>% 
      select(-c(m_ac, p_c))
  } else {
    samples <- samples %>% 
      select(-c(a_p_s, a_r_s, m_ac, p_c))
  }
  return(samples)
}

#apply feasibility bounds to parameter set samples (multipliers >/< 1, probs in 0-1, transitions out < 1)
apply_bounds <- function(samples, params_fixed) {
  samples <- samples %>% 
    mutate(flag0=1*((p_m>1) + (p_s>1) + (r_m>1) + (r_s>1) + (c_sp>1) + (c_tx>1) + (m_tb>1) +
                      (p_m<0) + (p_s<0) + (r_m<0) + (r_s<0) + (c_sp<0) + (c_tx<0) + (m_tb<0) +
                      (a_p_m<1) + (a_r_m>1) + (a_m<1) + (a_r_m<0) + (a_m<0)))
  #apply flag for a_tx too
  samples <- samples %>% 
    mutate(flag0=if_else("a_tx" %in% colnames(samples), flag0+1*(a_tx<1), flag0))
  #if a_p_s and a_r_s are also sampled (e.g. when RR_free is 1) then apply flags for those 2 parameters too
  samples <- samples %>% 
    mutate(flag0=if_else("a_p_s" %in% colnames(samples), flag0+1*(a_p_s<1), flag0))
  samples <- samples %>% 
    mutate(flag0=if_else("a_r_s" %in% colnames(samples), flag0+1*((a_r_s>1) + (a_r_s<0)), flag0))
  samples <- apply_flags(samples)
  samples <- samples %>% mutate(flag_sum=flag0+flag1+flag2+flag3+flag4+flag5)
  samples <- samples %>% filter(flag_sum==0) %>% 
    select(-starts_with("flag"))
  return(samples)
}

#sampling function for use in IMIS package
sample.prior <- function(n_samples) {
  i <- 0
  priors_all <- data.frame()
  while(i<n_samples) {
    priors <- sample_priors(n_samples, priors_prev_lb, priors_prev_ub, 
                            params_fixed_prev, RR_free, flag_symptom_dur)
    priors_all <- bind_rows(priors, priors_all)
    i <- nrow(priors_all)
  }
  priors_all <- priors_all[1:n_samples,]
  priors_all <- as.matrix(priors_all) #IMIS package requires matrix or vector
  return(priors_all)
}

#calculate prior likelihood for use in IMIS package
prior <- function(params) {
  params <- data.frame(params)
  like <- sapply(names(params_calib_prev), 
                 function(x) dunif(params[[x]], priors_prev_lb[[x]], 
                                   priors_prev_ub[[x]]), simplify=F, USE.NAMES=T)
  like <- bind_cols(like)
  flags <- apply_flags(params, params_fixed_prev, RR_free, flag_symptom_dur) %>% select(starts_with("flag"))
  #getting flagged = likelihood of 0, so change 1s to 0s and 0s to 1s, then take product
  flags <- flags*-1 + 1
  like <- cbind(like, flags)
  like <- rowProds(as.matrix(like))
  return(like)
}

#calculate model outputs that correspond to prevalence survey calibration targets
#version used for Philippines and Cambodia that calculates smear+ PNR instead of bac+ PNR
calc_outputs_prev <- function(sim_pop, p_opt, t, cyc_len, RR_free) { #curr_pop is sim_pop from a single timestep - slightly speeds indexing?
  params_depend <- c()
  if(RR_free==0) {
    params_depend <- c(params_depend,
                       "a_p_s"=p_opt[["a_p_m"]], 
                       "a_r_s"=p_opt[["a_r_m"]]
    )
  }
  p <- c(p_opt, params_depend)
  #prevalence survey targets are cross-sectional - no need for cycle length adjustment. 
  #notifications/deaths need cycle length adjustment (modeled deaths are cumulative so adjust for that too)
  outputs <- c("prop_m"=sim_pop[t+1,"tb_m"]/sum(sim_pop[t, tb_states]),
               "prop_s"=sim_pop[t+1,"tb_s"]/sum(sim_pop[t, tb_states]),
               "prop_ms"=sim_pop[t+1,"tb_ms"]/sum(sim_pop[t+1,tb_states]),
               "pnr_m_all"=(sim_pop[t+1,"tb_m"]+sim_pop[t+1,"tb_ms"])/
                 sum((sim_pop[((t+1)-(1/cyc_len)):(t+1),"tb_ms"])*p[["c_tx"]]*p[["a_tx"]]), #smear-positive divided by smear-positive that get treated
               "deaths_tb"=(sim_pop[t+1,"died_tb"]-sim_pop[(t+1)-(1/cyc_len),"died_tb"])/sum(sim_pop[t+1,tb_states]), #incremental tb deaths divided by tb cases
               "prop_m_notif"=(sim_pop[t+1, "tb_ms"]*p[["c_tx"]]*p[["a_tx"]])/
                 (sim_pop[t+1, "tb_ms"]*p[["c_tx"]]*p[["a_tx"]] + 
                    sim_pop[t+1, "tb_s"]*p[["c_tx"]]) #proportion of notifications that are smear-positive
  )
  return(outputs)
}

#version used for Vietnam, Nepal, and Bangladesh that calculates bacteria+ PNR instead of smear+
calc_outputs_prev_vnm <- function(sim_pop, p_opt, t, cyc_len, RR_free) { #curr_pop is sim_pop from a single timestep - slightly speeds indexing?
  params_depend <- c()
  if(RR_free==0) {
    params_depend <- c(params_depend,
                       "a_p_s"=p_opt[["a_p_m"]], 
                       "a_r_s"=p_opt[["a_r_m"]]
    )
  }
  p <- c(p_opt, params_depend)
  #prevalence survey targets are cross-sectional - no need for cycle length adjustment. 
  #notifications/deaths need cycle length adjustment (modeled deaths are cumulative so adjust for that too)
  outputs <- c("prop_m"=sim_pop[t+1,"tb_m"]/sum(sim_pop[t, tb_states]),
               "prop_s"=sim_pop[t+1,"tb_s"]/sum(sim_pop[t, tb_states]),
               "prop_ms"=sim_pop[t+1,"tb_ms"]/sum(sim_pop[t+1,tb_states]),
               "pnr_all"=sum(sim_pop[t+1,tb_states])/
                 (sum((sim_pop[((t+1)-(1/cyc_len)):(t+1),"tb_ms"])*p[["c_tx"]]*p[["a_tx"]]) +
                    sum((sim_pop[((t+1)-(1/cyc_len)):(t+1),"tb_s"])*p[["c_tx"]])), #total divided by total that get treated
               "deaths_tb"=(sim_pop[t+1,"died_tb"]-sim_pop[(t+1)-(1/cyc_len),"died_tb"])/sum(sim_pop[t+1,tb_states]), #incremental tb deaths divided by tb cases
               "prop_m_notif"=(sim_pop[t+1, "tb_ms"]*p[["c_tx"]]*p[["a_tx"]])/
                 (sim_pop[t+1, "tb_ms"]*p[["c_tx"]]*p[["a_tx"]] + 
                    sim_pop[t+1, "tb_s"]*p[["c_tx"]]) #proportion of notifications that are smear-positive
  )
  return(outputs)
}

#calculate model outputs that correspond to historical cohort calibration targets
calc_outputs_hist <- function(sim_pop, cyc_len, calib_type) { 
  #adjust timing (5 yrs and 10 yrs) for cycle length
  dead_5yr=sim_pop[5/cyc_len, "died_tb"] + sim_pop[5/cyc_len, "died_nontb"]
  dead_10yr=sim_pop[10/cyc_len, "died_tb"] + sim_pop[10/cyc_len, "died_nontb"]
  if (calib_type=="hist_neg") {
    outputs <- c("tb_s_dead_5yr"=dead_5yr, "tb_s_dead_10yr"=dead_10yr)
  } else if(calib_type=="hist_pos") {
    outputs <- c("tb_ms_dead_5yr"=dead_5yr, "tb_ms_dead_10yr"=dead_10yr)
  }
  return(outputs)
}

#version that also calculte % still smear+ at 4 years (if calib_type is hist_pos only)
calc_outputs_hist_smear <- function(sim_pop, cyc_len, calib_type) { 
  #adjust timing (5 yrs and 10 yrs) for cycle length
  dead_5yr=sim_pop[5/cyc_len, "died_tb"] + sim_pop[5/cyc_len, "died_nontb"]
  dead_10yr=sim_pop[10/cyc_len, "died_tb"] + sim_pop[10/cyc_len, "died_nontb"]
  if (calib_type=="hist_neg") {
    outputs <- c("tb_s_dead_5yr"=dead_5yr, "tb_s_dead_10yr"=dead_10yr)
  } else if(calib_type=="hist_pos") {
    smear_4yr=(sim_pop[4/cyc_len, "tb_ms"] + sim_pop[4/cyc_len, "tb_m"])/
      sum(sim_pop[4/cyc_len, alive_states])
    outputs <- c("tb_ms_dead_5yr"=dead_5yr, "tb_ms_dead_10yr"=dead_10yr,
                 "tb_smear_4yr"=smear_4yr)
  }
  return(outputs)
}

#function that runs the model with a set of parameters and returns model outputs (vs. targets) and penalties
calib_out <- function(params_calib, params_fixed, calib_type,  
                      RR_free, smear_hist_calib, country, t_end, cyc_len, sim_pop) {  
  p_use <- c(params_calib, params_fixed)
  params_depend <- c()
  #parameter dependencies
  if(RR_free==0) {
    params_depend <- c(params_depend,
                       "a_p_s"=p_use[["a_p_m"]], 
                       "a_r_s"=p_use[["a_r_m"]]
    )
  }
  p_use <- c(p_use, params_depend)
  for(t in 1:t_end) {
    curr_pop <- nat_hist_markov(p_use, sim_pop[t,], t)
    sim_pop <- bind_rows(sim_pop, curr_pop)
  }
  if(calib_type=="prev" & country %in% c("Philippines", "Cambodia")) {
    outputs <- calc_outputs_prev(sim_pop, p_use, t, cyc_len, RR_free)
  } else if(calib_type=="prev" & country %in% c("Vietnam", "Nepal", "Bangladesh")) {
    outputs <- calc_outputs_prev_vnm(sim_pop, p_use, t, cyc_len, RR_free)
  } else if((calib_type=="hist_pos"|calib_type=="hist_neg") & smear_hist_calib==0) {
    outputs <- calc_outputs_hist(sim_pop, cyc_len, calib_type)
  } else if((calib_type=="hist_pos"|calib_type=="hist_neg") & smear_hist_calib==1) {
    outputs <- calc_outputs_hist_smear(sim_pop, cyc_len, calib_type)
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

#function to calculate the likelihood of observing calibration targets given model output
calc_like <- function(out, tr, tr_lb, tr_ub, mort_samples, prev_cases,
                      prop_m_notif_smooth, pnr_params, calib_type, country) { #outputs, targets, and upper/lower confidence bounds on targets
  if(calib_type=="prev") {
    #prevalence survey targets proportion of infections by smear/symptom status - we have actual sample size
    prop_m <- dbinom(round(tr[["prop_m"]]*prev_cases), size=prev_cases, prob=out[["prop_m"]], log=T)
    prop_s <- dbinom(round(tr[["prop_s"]]*prev_cases), size=prev_cases, prob=out[["prop_s"]], log=T)
    prop_ms <- dbinom(round(tr[["prop_ms"]]*prev_cases), size=prev_cases, prob=out[["prop_ms"]], log=T)
    #prevalence to notification ratio: 0 to inf - gamma fits well (parameters estimated using dampack gamma_params)
    if(country %in% c("Philippines", "Cambodia")) {
      pnr_m_all <- dgamma(out[["pnr_m_all"]], shape=pnr_params$pnr_gamma_shape, scale=pnr_params$pnr_gamma_scale, log=T)
    } else if(country %in% c("Vietnam", "Nepal", "Bangladesh")) {
      pnr_all <- dgamma(out[["pnr_all"]], shape=pnr_params$pnr_gamma_shape, scale=pnr_params$pnr_gamma_scale, log=T)
    }
    #use empirical distribution for the TB mortality target
    deaths_tb <- unname(log(mort_samples[as.character(round(out[["deaths_tb"]]*1000))]))
    #make likelihood very very small (and decreasing) if > max of all mort samples (min is 0 so no need to do this on low end)
    deaths_tb[(is.na(deaths_tb)|deaths_tb==-Inf) & !is.na(out[["deaths_tb"]])] <- 
        unlist(log(1/abs(round(out[(is.na(deaths_tb)|deaths_tb==-Inf) & !is.na(out[["deaths_tb"]]), 
                            "deaths_tb"]*1000)-max(unname(mort_samples)))/1000000)) 
    #proportion of notifications that are smear-positive - use empirical distribution
    prop_m_notif <- unname(log(prop_m_notif_smooth[as.character(round(out[["prop_m_notif"]]*100))]))
    
    if(country %in% c("Philippines", "Cambodia")) {
      log_like_all <- data.frame("prop_m"=prop_m,
                                 "prop_s"=prop_s,
                                 "prop_ms"=prop_ms,
                                 "pnr_m_all"=pnr_m_all,
                                 "deaths_tb"=deaths_tb,
                                 "prop_m_notif"=prop_m_notif)
    } else if(country %in% c("Vietnam", "Nepal", "Bangladesh")) {
      log_like_all <- data.frame("prop_m"=prop_m,
                                 "prop_s"=prop_s,
                                 "prop_ms"=prop_ms,
                                 "pnr_all"=pnr_all,
                                 "deaths_tb"=deaths_tb,
                                 "prop_m_notif"=prop_m_notif)
    }
  } else if(calib_type=="hist_pos") {
    #historical mortality targets - sizes (n) of binomial distributions established to match CIs from meta-regression
    tb_ms_dead_5yr <- dbinom(round(tr[["tb_ms_dead_5yr"]]*220), size=220, prob=out[["tb_ms_dead_5yr"]], log=T)
    tb_ms_dead_10yr <- dbinom(round(tr[["tb_ms_dead_10yr"]]*220), size=220, prob=out[["tb_ms_dead_10yr"]], log=T)
    log_like_all <- data.frame("tb_ms_dead_5yr"=tb_ms_dead_5yr,
                  "tb_ms_dead_10yr"=tb_ms_dead_10yr)
  } else if(calib_type=="hist_neg") {
    #historical mortality targets - sizes (n) of binomial distributions established to match CIs from meta-regression
    tb_s_dead_5yr <- dbinom(round(tr[["tb_s_dead_5yr"]]*325), size=325, prob=out[["tb_s_dead_5yr"]], log=T)
    tb_s_dead_10yr <- dbinom(round(tr[["tb_s_dead_10yr"]]*200), size=200, prob=out[["tb_s_dead_10yr"]], log=T)
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

#version without the 10-year mortality targets
#currently only works for Philippines
calc_like_no10 <- function(out, tr, tr_lb, tr_ub, mort_samples, prev_cases,
                           prop_m_notif_smooth, pnr_params, calib_type) { #outputs, targets, and upper/lower confidence bounds on targets
  if(calib_type=="prev") {
    #prevalence survey targets proportion of infections by smear/symptom status - we have actual sample size
    prop_m <- dbinom(round(tr[["prop_m"]]*prev_cases), size=prev_cases, prob=out[["prop_m"]], log=T)
    prop_s <- dbinom(round(tr[["prop_s"]]*prev_cases), size=prev_cases, prob=out[["prop_s"]], log=T)
    prop_ms <- dbinom(round(tr[["prop_ms"]]*prev_cases), size=prev_cases, prob=out[["prop_ms"]], log=T)
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
    
    log_like_all <- data.frame("prop_m"=prop_m,
                               "prop_s"=prop_s,
                               "prop_ms"=prop_ms,
                               "pnr_m_all"=pnr_m_all,
                               "deaths_tb"=deaths_tb,
                               "prop_m_notif"=prop_m_notif)
  } else if(calib_type=="hist_pos") {
    #historical mortality targets - sizes (n) of binomial distributions established to match CIs from meta-regression
    tb_ms_dead_5yr <- dbinom(round(tr[["tb_ms_dead_5yr"]]*220), size=220, prob=out[["tb_ms_dead_5yr"]], log=T)
    log_like_all <- data.frame("tb_ms_dead_5yr"=tb_ms_dead_5yr)
  } else if(calib_type=="hist_neg") {
    #historical mortality targets - sizes (n) of binomial distributions established to match CIs from meta-regression
    tb_s_dead_5yr <- dbinom(round(tr[["tb_s_dead_5yr"]]*325), size=325, prob=out[["tb_s_dead_5yr"]], log=T)
    log_like_all <- data.frame("tb_s_dead_5yr"=tb_s_dead_5yr)
  } else {
    print("error: incorrect calibration type")
  }
  log_like <- rowSums(log_like_all)
  like <- rowProds(as.matrix(exp(log_like_all)))
  like_out <- list("like"=like, "log_like"=log_like)
  return(like_out)
}

#version with historical target on bacillary status over time
calc_like_smear_hist <- function(out, tr, tr_lb, tr_ub, mort_samples, prev_cases,
                                 prop_m_notif_smooth, pnr_params, calib_type,
                                 country) { #outputs, targets, and upper/lower confidence bounds on targets
  if(calib_type=="prev") {
    #prevalence survey targets proportion of infections by smear/symptom status - we have actual sample size
    prop_m <- dbinom(round(tr[["prop_m"]]*prev_cases), size=prev_cases, prob=out[["prop_m"]], log=T)
    prop_s <- dbinom(round(tr[["prop_s"]]*prev_cases), size=prev_cases, prob=out[["prop_s"]], log=T)
    prop_ms <- dbinom(round(tr[["prop_ms"]]*prev_cases), size=prev_cases, prob=out[["prop_ms"]], log=T)
    #prevalence to notification ratio: 0 to inf - gamma fits well (parameters estimated using dampack gamma_params)
    if(country %in% c("Philippines", "Cambodia")) {
      pnr_m_all <- dgamma(out[["pnr_m_all"]], shape=pnr_params$pnr_gamma_shape, scale=pnr_params$pnr_gamma_scale, log=T)
    } else if(country %in% c("Vietnam", "Nepal", "Bangladesh")) {
      pnr_all <- dgamma(out[["pnr_all"]], shape=pnr_params$pnr_gamma_shape, scale=pnr_params$pnr_gamma_scale, log=T)
    }
    #use empirical distribution for the TB mortality target
    deaths_tb <- unname(log(mort_samples[as.character(round(out[["deaths_tb"]]*1000))]))
    #make likelihood very very small (and decreasing) if > max of all mort samples (min is 0 so no need to do this on low end)
    deaths_tb[(is.na(deaths_tb)|deaths_tb==-Inf) & !is.na(out[["deaths_tb"]])] <- 
      unlist(log(1/abs(round(out[(is.na(deaths_tb)|deaths_tb==-Inf) & !is.na(out[["deaths_tb"]]), 
                                 "deaths_tb"]*1000)-max(unname(mort_samples)))/1000000)) 
    #proportion of notifications that are smear-positive - use empirical distribution
    prop_m_notif <- unname(log(prop_m_notif_smooth[as.character(round(out[["prop_m_notif"]]*100))]))
    
    if(country %in% c("Philippines", "Cambodia")) {
      log_like_all <- data.frame("prop_m"=prop_m,
                                 "prop_s"=prop_s,
                                 "prop_ms"=prop_ms,
                                 "pnr_m_all"=pnr_m_all,
                                 "deaths_tb"=deaths_tb,
                                 "prop_m_notif"=prop_m_notif)
    } else if(country %in% c("Vietnam", "Nepal", "Bangladesh")) {
      log_like_all <- data.frame("prop_m"=prop_m,
                                 "prop_s"=prop_s,
                                 "prop_ms"=prop_ms,
                                 "pnr_all"=pnr_all,
                                 "deaths_tb"=deaths_tb,
                                 "prop_m_notif"=prop_m_notif)
    }
  } else if(calib_type=="hist_pos") {
    #historical mortality targets - sizes (n) of binomial distributions established to match CIs from meta-regression
    tb_ms_dead_5yr <- dbinom(round(tr[["tb_ms_dead_5yr"]]*220), size=220, prob=out[["tb_ms_dead_5yr"]], log=T)
    tb_ms_dead_10yr <- dbinom(round(tr[["tb_ms_dead_10yr"]]*220), size=220, prob=out[["tb_ms_dead_10yr"]], log=T)
    #historical target on bacillary status over time of those alive - combining sample sizes across the 4 studies
    tb_ms_smearneg_4yr_1 <- dbinom(tr[["tb_smear_4yr_1"]], size=tr[["alive_4yr_1"]], prob=out[["tb_smear_4yr"]], log=T)
    tb_ms_smearneg_4yr_2 <- dbinom(tr[["tb_smear_4yr_2"]], size=tr[["alive_4yr_2"]], prob=out[["tb_smear_4yr"]], log=T)
    tb_ms_smearneg_4yr_3 <- dbinom(tr[["tb_smear_4yr_3"]], size=tr[["alive_4yr_3"]], prob=out[["tb_smear_4yr"]], log=T)
    log_like_all <- data.frame("tb_ms_dead_5yr"=tb_ms_dead_5yr,
                               "tb_ms_dead_10yr"=tb_ms_dead_10yr,
                               "tb_ms_smearneg_4yr_1"=tb_ms_smearneg_4yr_1,
                               "tb_ms_smearneg_4yr_2"=tb_ms_smearneg_4yr_2,
                               "tb_ms_smearneg_4yr_3"=tb_ms_smearneg_4yr_3)
  } else if(calib_type=="hist_neg") {
    #historical mortality targets - sizes (n) of binomial distributions established to match CIs from meta-regression
    tb_s_dead_5yr <- dbinom(round(tr[["tb_s_dead_5yr"]]*325), size=325, prob=out[["tb_s_dead_5yr"]], log=T)
    tb_s_dead_10yr <- dbinom(round(tr[["tb_s_dead_10yr"]]*200), size=200, prob=out[["tb_s_dead_10yr"]], log=T)
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

#function to calculate sum of squared errors between model output and calibration targets (including penalties)
calc_sse <- function(out, tr, scale_fac) { #function of model outputs and calibration targets
  sse <- sum((out$outputs/scale_fac-tr/scale_fac)^2) + out$penalties
  #print(sse)
  if(sse==Inf|is.na(sse)) {
    sse <- out$penalties + 9999999
  }
  return(sse)
}


#calculate all model outputs and likelihoods (across the 3 types of calibrations)
output_like <- function(params) {
  if(is.vector(params)) {
    params <- data.frame(t(params))
  } else {
    params <- data.frame(params)
  }
  n_time <- 200
  t_end <- n_time/cyc_len
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
  targets <- pull_targets("prev", targets_all, country)[[1]]
  names <- pull_targets("prev", targets_all, country)[[2]]
  sim_pop <- data.frame() 
  sim_pop <- bind_rows(sim_pop, start_pop)
  params_use <- params %>% select(names(params_calib_prev))
  out_prev <- lapply(1:nrow(params_use), function(x)
    calib_out(params_use[x,], params_fixed_prev, "prev", RR_free, 
              smear_hist_calib, country,
              t_end, cyc_len, sim_pop)$outputs)
  out_prev <- bind_rows(out_prev)
  like_prev <- calc_like(out_prev, targets, 
                         targets_all_lb, targets_all_ub, 
                         mort_samples, prev_cases,
                         prop_m_notif_smooth,
                         pnr_params, "prev", country)
  
  #HISTORICAL COHORT SMEAR POSITIVE TARGETS
  n_time <- 10
  t_end <- n_time/cyc_len
  start_pop <- c("t"=0, 
                 "tb"=0, #smear- symptom- TB
                 "tb_m"=0, #smear+ symptom- TB
                 "tb_s"=0, #smear- symptom+ TB
                 "tb_ms"=1, #smear+ symptom+ TB
                 "sp_cure"=0, #spontaneously cured
                 "tx_cure"=0, #cured via diagnosis and treatment
                 "died_tb"=0, #TB death
                 "died_nontb"=0 #non-TB death
  )
  targets <- pull_targets("hist_pos", targets_all, country)[[1]]
  names <- pull_targets("hist_pos", targets_all, country)[[2]]
  if(smear_hist_calib==1) {
    targets[["tb_smear_4yr_1"]] <- targets_all[["tb_smear_4yr_1"]]
    targets[["alive_4yr_1"]] <- targets_all[["alive_4yr_1"]]
    targets[["tb_smear_4yr_2"]] <- targets_all[["tb_smear_4yr_2"]]
    targets[["alive_4yr_2"]] <- targets_all[["alive_4yr_2"]]
    targets[["tb_smear_4yr_3"]] <- targets_all[["tb_smear_4yr_3"]]
    targets[["alive_4yr_3"]] <- targets_all[["alive_4yr_3"]] 
  }
  sim_pop <- data.frame() 
  sim_pop <- bind_rows(sim_pop, start_pop)
  params_use <- params %>% select(names(params_calib_hist))
  out_hist_pos <- lapply(1:nrow(params_use), function(x)
    calib_out(params_use[x,], params_fixed_hist, "hist_pos", RR_free, 
              smear_hist_calib, country,
              t_end, cyc_len, sim_pop)$outputs)
  out_hist_pos <- bind_rows(out_hist_pos)
  like_hist_pos <- calc_like(out_hist_pos, targets, 
                             targets_all_lb, targets_all_ub, 
                             mort_samples, prev_cases, prop_m_notif_smooth,
                             pnr_params, "hist_pos", country)
  
  #HISTORICAL COHORT SMEAR NEGATIVE TARGETS
  n_time <- 10
  t_end <- n_time/cyc_len
  start_pop <- c("t"=0, 
                 "tb"=0, #smear- symptom- TB
                 "tb_m"=0, #smear+ symptom- TB
                 "tb_s"=1, #smear- symptom+ TB
                 "tb_ms"=0, #smear+ symptom+ TB
                 "sp_cure"=0, #spontaneously cured
                 "tx_cure"=0, #cured via diagnosis and treatment
                 "died_tb"=0, #TB death
                 "died_nontb"=0 #non-TB death
  )
  targets <- pull_targets("hist_neg", targets_all, country)[[1]]
  names <- pull_targets("hist_neg", targets_all, country)[[2]]
  sim_pop <- data.frame() 
  sim_pop <- bind_rows(sim_pop, start_pop)
  params_use <- params %>% select(names(params_calib_hist))
  out_hist_neg <- lapply(1:nrow(params_use), function(x)
    calib_out(params_use[x,], params_fixed_hist, "hist_neg", RR_free, 
              smear_hist_calib, country,
              t_end, cyc_len, sim_pop)$outputs)
  out_hist_neg <- bind_rows(out_hist_neg)
  like_hist_neg <- calc_like(out_hist_neg, targets, 
                             targets_all_lb, targets_all_ub, 
                             mort_samples, prev_cases, prop_m_notif_smooth,
                             pnr_params, "hist_neg", country)
  
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

#likelihood wrapper function for IMIS package
likelihood <- function(params) {
  out_all <- output_like(params)
  return(out_all$like)
}
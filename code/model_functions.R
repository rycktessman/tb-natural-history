#DEFINE STATES AND OUTCOMES
states <- c("tb", "tb_m", "tb_s", "tb_ms", "sp_cure", "tx_cure", "died_tb", "died_nontb")
alive_states <- c("tb", "tb_m", "tb_s", "tb_ms", "sp_cure", "tx_cure")
tb_states <- c("tb", "tb_m", "tb_s", "tb_ms")
outcomes <- c("died_tb", "died_nontb")

#GENERATE TRANSITION MATRIX FOR MARKOV MODEL - if inflows off, rowSums=1's
gen_trans_mat <- function(params) {
  trans_mat <- matrix(data=0, nrow=8, ncol=8) #from=rows, to=columns
  #order: TB deaths, non-TB deaths, treatment, progression & regression, inflows
  #transtitions out of smear-/symptom-
  trans_mat[1, 2] <- params$p_m*(1-params$m_ac)
  trans_mat[1, 3] <- params$p_s*(1-params$m_ac)
  trans_mat[1, 5] <- params$c_sp*(1-params$m_ac)
  trans_mat[1, 8] <- params$m_ac
  trans_mat[1, 1] <- 1 - ((1-params$m_ac)*(params$p_m + params$p_s + params$c_sp) + params$m_ac) + params$m_ac*(params$inflows==1)
  #transitions out of smear+/symptom-
  trans_mat[2, 1] <- params$r_m*(1-params$m_ac) + params$m_ac*(params$inflows==1)
  trans_mat[2, 4] <- params$p_s*params$a_p_s*(1-params$m_ac)
  trans_mat[2, 8] <- params$m_ac
  trans_mat[2, 2] <- 1- ((1-params$m_ac)*(params$r_m + params$p_s*params$a_p_s) + params$m_ac)
  #transitions out of smear-/symptom+
  trans_mat[3, 1] <- params$r_s*(1-params$m_tb)*(1-params$m_ac)*(1-params$c_tx) + 
    (params$m_tb + params$m_ac*(1-params$m_tb))*(params$inflows==1)
  trans_mat[3, 4] <- params$p_m*params$a_p_m*(1-params$m_tb)*(1-params$m_ac)*(1-params$c_tx)
  trans_mat[3, 6] <- params$c_tx*(1-params$m_tb)*(1-params$m_ac)
  trans_mat[3, 7] <- params$m_tb
  trans_mat[3, 8] <- params$m_ac*(1-params$m_tb)
  trans_mat[3, 3] <- 1 - ((params$r_s + params$p_m*params$a_p_m)*(1-params$m_tb)*(1-params$m_ac)*(1-params$c_tx) +
                            params$c_tx*(1-params$m_tb)*(1-params$m_ac) + 
                            params$m_tb + params$m_ac*(1-params$m_tb))
  #transitions out of smear+/symptom+
  trans_mat[4, 1] <- (params$m_tb*params$a_m + params$m_ac*(1-params$m_tb*params$a_m))*(params$inflows==1)
  trans_mat[4, 2] <- params$r_s*params$a_r_s*(1-params$m_tb*params$a_m)*(1-params$m_ac)*(1-params$c_tx*params$a_tx)
  trans_mat[4, 3] <- params$r_m*params$a_r_m*(1-params$m_tb*params$a_m)*(1-params$m_ac)*(1-params$c_tx*params$a_tx)
  trans_mat[4, 6] <- params$c_tx*params$a_tx*(1-params$m_tb*params$a_m)*(1-params$m_ac)
  trans_mat[4, 7] <- params$m_tb*params$a_m
  trans_mat[4, 8] <- params$m_ac*(1-params$m_tb*params$a_m)
  trans_mat[4, 4] <- 1- ((1-params$m_tb*params$a_m)*(1-params$m_ac)*(1-params$c_tx*params$a_tx)*(params$r_s*params$a_r_s + params$r_m*params$a_r_m) + 
                           params$c_tx*params$a_tx*(1-params$m_tb*params$a_m)*(1-params$m_ac) + 
                           params$m_tb*params$a_m + params$m_ac*(1-params$m_tb*params$a_m))
  #transitions out of spontaneously resolved
  trans_mat[5, 8] <- params$m_ac
  trans_mat[5, 1] <- params$p_c*(1-params$m_ac) + params$m_ac*(params$inflows==1)
  trans_mat[5, 5] <- 1 - (params$m_ac + params$p_c*(1-params$m_ac))
  #transitions out of detected/treated
  trans_mat[6, 8] <- params$m_ac
  trans_mat[6, 1] <- params$m_ac*(params$inflows==1)
  trans_mat[6, 6] <- 1 - params$m_ac 
  #deaths stay dead
  trans_mat[7, 7] <- 1
  trans_mat[8, 8] <- 1
  return(trans_mat)
}

#MICROSIM MODEL FUNCTIONS
#1=tb (smear- symptom-)
#2=tb_m (smear+ symptom-)
#3=tb_s (smear- symptom+)
#4=tb_ms (smear+ symptom+)
#5=sp_cure (spontaneously cured)
#6=tx_cure (cured w/ diagnosis & treatment)
#7=died_tb (TB death)
#8=died_nontb (non-TB death)

#function to sample all random numbers up front (for variance reduction)
micro_sample <- function(p, n, t) {
  #draw from multinomial distributions to determine who transitions during each transition
  #transitions from 1 (tb, smear- symptom-)
  tb_trans <- rmultinom(n=n*t, size=1, 
                        prob=c(p$p_m, p$p_s, p$c_sp, p$m_ac, 1-(p$p_m+p$p_s+p$c_sp+p$m_ac)))
  row.names(tb_trans) <- c("tb_to_tb_m", "tb_to_tb_s", "tb_to_sp_cure", "tb_to_died_nontb", "tb_stay")
  
  #transitions from 2 (tb_m, smear+ symptom-)
  tb_m_trans <- rmultinom(n=n*t, size=1, 
                          prob=c(p$r_m, p$p_s*p$a_p_s, p$m_ac, 1-(p$r_m+p$p_s*p$a_p_s+p$m_ac)))
  row.names(tb_m_trans) <- c("tb_m_to_tb", "tb_m_to_tb_ms", "tb_m_to_died_nontb", "tb_m_stay")
  
  #transitions from 3 (tb_s, smear- symptom+)
  tb_s_trans <- rmultinom(n=n*t, size=1,
                          prob=c(p$r_s, p$p_m*p$a_p_m, p$c_tx, p$m_ac, p$m_tb, 
                                 1-(p$r_s+p$p_m*p$a_p_m+p$c_tx+p$m_ac+p$m_tb)))
  row.names(tb_s_trans) <- c("tb_s_to_tb", "tb_s_to_tb_ms", "tb_s_to_tx_cure", "tb_s_to_died_nontb",
                             "tb_s_to_died_tb", "tb_s_stay")
  
  #transitions from 4 (tb_ms, smear- symptom+)
  tb_ms_trans <- rmultinom(n=n*t, size=1, 
                           prob=c(p$r_s*p$a_r_s, p$r_m*p$a_r_m, p$c_tx*p$a_tx, p$m_ac, p$m_tb*p$a_m,
                                  1-(p$r_s*p$a_r_s+p$r_m*p$a_r_m+p$c_tx*p$a_tx+p$m_ac+p$m_tb*p$a_m)))
  row.names(tb_ms_trans) <- c("tb_ms_to_tb_m", "tb_ms_to_tb_s", "tb_ms_to_tx_cure", 
                              "tb_ms_to_died_nontb", "tb_ms_to_died_tb", "tb_ms_stay")
  
  #transitions from 5 (sp_cure, spontaneously cured)
  sp_cure_trans <- rmultinom(n=n*t, size=1, prob=c(p$p_c, p$m_ac, 1-(p$p_c+p$m_ac)))
  row.names(sp_cure_trans) <- c("sp_cure_to_tb", "sp_cure_to_died_nontb", "sp_cure_stay")
  
  #transitions from 6 (tx_cure, cured via diagnosis and treatment)
  tx_cure_to_died_nontb <- rbinom(n=n*t, size=1, prob=p$m_ac)
  tx_cure_stay <- rep(1, n*t) - tx_cure_to_died_nontb
  
  #no transitions from 7 and 8 (dead states)
  #combine into dataframe for use in microsim model
  rands <- data.frame("tb_to_tb_m"=tb_trans["tb_to_tb_m",],
                      "tb_to_tb_s"=tb_trans["tb_to_tb_s",],
                      "tb_to_sp_cure" = tb_trans["tb_to_sp_cure",],
                      "tb_to_died_nontb" = tb_trans["tb_to_died_nontb",],
                      "tb_stay" = tb_trans["tb_stay",],
                      "tb_m_to_tb" = tb_m_trans["tb_m_to_tb",],
                      "tb_m_to_tb_ms" = tb_m_trans["tb_m_to_tb_ms",],
                      "tb_m_to_died_nontb" = tb_m_trans["tb_m_to_died_nontb",],
                      "tb_m_stay" = tb_m_trans["tb_m_stay",],
                      "tb_s_to_tb"= tb_s_trans["tb_s_to_tb",],
                      "tb_s_to_tb_ms" = tb_s_trans["tb_s_to_tb_ms",],
                      "tb_s_to_tx_cure" = tb_s_trans["tb_s_to_tx_cure",], 
                      "tb_s_to_died_nontb" = tb_s_trans["tb_s_to_died_nontb",],
                      "tb_s_to_died_tb" = tb_s_trans["tb_s_to_died_tb",],
                      "tb_s_stay" = tb_s_trans["tb_s_stay",],
                      "tb_ms_to_tb_m" = tb_ms_trans["tb_ms_to_tb_m",],
                      "tb_ms_to_tb_s" = tb_ms_trans["tb_ms_to_tb_s",],
                      "tb_ms_to_tx_cure" = tb_ms_trans["tb_ms_to_tx_cure",], 
                      "tb_ms_to_died_nontb" = tb_ms_trans["tb_ms_to_died_nontb",],
                      "tb_ms_to_died_tb" = tb_ms_trans["tb_ms_to_died_tb",],
                      "tb_ms_stay" = tb_ms_trans["tb_ms_stay",],
                      "sp_cure_to_tb" = sp_cure_trans["sp_cure_to_tb",],
                      "sp_cure_to_died_nontb" = sp_cure_trans["sp_cure_to_died_nontb",],
                      "sp_cure_stay" = sp_cure_trans["sp_cure_stay",],
                      "tx_cure_to_died_nontb" = tx_cure_to_died_nontb,
                      "tx_cure_stay" = tx_cure_stay
  )
  return(rands)
}

#function to run model for a single timestep - individual-level version
nat_hist_micro <- function(pop, r, p)  { #sim_pop row, random numbers subset, params
  curr_pop <- (pop==1)*(
    2*(r$tb_to_tb_m==1) +
      3*(r$tb_to_tb_s==1) +
      5*(r$tb_to_sp_cure==1) +
      8*(r$tb_to_died_nontb==1) +
      1*(r$tb_stay==1)
  ) + (pop==2)*( 
    1*(r$tb_m_to_tb==1) +
      4*(r$tb_m_to_tb_ms==1) +
      8*(r$tb_m_to_died_nontb==1) +
      2*(r$tb_m_stay==1)
  ) + (pop==3)*(
    1*(r$tb_s_to_tb==1) +
      4*(r$tb_s_to_tb_ms==1) +
      6*(r$tb_s_to_tx_cure==1) +
      7*(r$tb_s_to_died_tb==1) +
      8*(r$tb_s_to_died_nontb==1) +
      3*(r$tb_s_stay==1)
  ) + (pop==4)*(
    2*(r$tb_ms_to_tb_m==1) +
      3*(r$tb_ms_to_tb_s==1) +
      6*(r$tb_ms_to_tx_cure==1) +
      7*(r$tb_ms_to_died_tb==1) +
      8*(r$tb_ms_to_died_nontb==1) +
      4*(r$tb_ms_stay==1)
  ) + (pop==5)*(
    1*(r$sp_cure_to_tb==1) +
      8*(r$sp_cure_to_died_nontb==1) +
      5*(r$sp_cure_stay==1)
  ) + (pop==6)*(
    8*(r$tx_cure_to_died_nontb==1) +
      6*(r$tx_cure_stay==1)
  ) + 
    if_else(p$inflows==1, 1, 7)*(pop==7) + if_else(p$inflows==1, 1, 8)*(pop==8) #all deaths re-enter the model as inflows on the next timestep
  return(curr_pop)
}

#CLEANING FUNCTIONS (apply to both markov and individual versions)
process_output <- function(sim_pop) {
  sim_pop_long <- pivot_longer(sim_pop, cols=2:ncol(sim_pop), names_to="outcome", values_to="value")
  sim_pop_long <- sim_pop_long %>% mutate(name=case_when(outcome=="tb"~"Smear- Symptom- TB",
                                                         outcome=="tb_m"~"Smear+ Symptom- TB",
                                                         outcome=="tb_s"~"Smear- Symptom+ TB",
                                                         outcome=="tb_ms"~"Smear+ Symptom+ TB",
                                                         outcome=="sp_cure"~"Spontaneous Cures",
                                                         outcome=="tx_cure"~"Diagnosed & Treated",
                                                         outcome=="died_tb"~"TB Deaths",
                                                         outcome=="died_nontb"~"Non-TB Deaths"))
  sim_pop_long <- sim_pop_long %>% group_by(t) %>% 
    mutate(value_scale=case_when(outcome %in% alive_states~value/sum(value[outcome %in% alive_states]), #calculate as % of those still alive in time step
                                 outcome %in% c("died_tb", "died_nontb")~value)) #already cumulative
  sim_pop_long <- sim_pop_long %>% group_by(t) %>% 
    mutate(prop_tb=if_else(outcome %in% tb_states, value/sum(value[outcome %in% tb_states]), 0))  
}
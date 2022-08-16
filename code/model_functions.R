#packages and working directory
library(dplyr)

#DEFINE STATES AND OUTCOMES
states <- c("tb", "tb_m", "tb_s", "tb_ms", "sp_cure", "tx_cure", "died_tb", "died_nontb")
alive_states <- c("tb", "tb_m", "tb_s", "tb_ms", "sp_cure", "tx_cure")
tb_states <- c("tb", "tb_m", "tb_s", "tb_ms")
outcomes <- c("died_tb", "died_nontb", "rel_inf")

#MARKOV MODEL FUNCTION
#function to run model for a single timestep - markov model version (not individual level)
#generate data frame where each row is a timestep, each column is a state (track proportion of pop by state), relative # secondary infections is also a column
nat_hist_markov <- function(p, pop, t) { #p=params list, pop=simulated population, t=timestep
  tb <- pop$tb*(1-(p$p_m + p$p_s + p$m_ac + p$c_sp)) + #smear- symptom- stay
    pop$tb_m*p$r_m + #smear status regresses
    pop$tb_s*p$r_s + #symptom status regresses
    pop$sp_cure*p$p_c + #progress from spontaneous cure
    (pop$tb + pop$tb_m + pop$tb_s + pop$tb_ms + pop$sp_cure + pop$tx_cure)*p$m_ac*(p$inflows==1) + #new inflows to the model = deaths out of the model
    pop$tb_s*p$m_tb*(p$inflows==1) +  #new inflows to the model = deaths out of the model
    pop$tb_ms*p$a_m*p$m_tb*(p$inflows==1) #new inflows to the model = deaths out of the model
  tb_m <- pop$tb_m*(1-(p$r_m + p$a_p_s*p$p_s + p$m_ac)) + #smear+ symptom- stay
    pop$tb*p$p_m + #smear status progresses
    pop$tb_ms*p$a_r_s*p$r_s #symptom status regresses
  tb_s <- pop$tb_s*(1-(p$r_s + p$a_p_m*p$p_m + p$m_tb + p$m_ac + p$c_tx)) + #smear- symptom+ stay  
    pop$tb*p$p_s + #symptom status progresses
    pop$tb_ms*p$a_r_m*p$r_m  #smear status regresses
  tb_ms <- pop$tb_ms*(1-(p$a_r_s*p$r_s + p$a_r_m*p$r_m + p$a_m*p$m_tb + p$m_ac + p$a_tx*p$c_tx)) + #smear+ symptom+ stay
    pop$tb_s*p$a_p_m*p$p_m + #smear status progresses
    pop$tb_m*p$a_p_s*p$p_s  #symptom status progresses
  sp_cure <- pop$sp_cure*(1-(p$p_c + p$m_ac)) + #spontaneous cure stay
    pop$tb*p$c_sp #spontaneous cures
  tx_cure <- pop$tx_cure*(1-p$m_ac) + #treated cured stay
    pop$tb_s*p$c_tx + #symptom+ smear- treated
    pop$tb_ms*p$a_tx*p$c_tx #symptom+ smear+ treated
  died_tb <- pop$died_tb + #dead stay
    pop$tb_s*p$m_tb + #symptom+ smear- die of TB
    pop$tb_ms*p$a_m*p$m_tb #symptom+ smear+ died of TB
  died_nontb <- pop$died_nontb + #dead stay 
    (pop$tb + pop$tb_m + pop$tb_s + pop$tb_ms + pop$sp_cure + pop$tx_cure)*p$m_ac #non-TB deaths from all states
  rel_inf <- p$i*pop$tb + p$i_s*pop$tb_s + p$i_m*pop$tb_m + p$i_ms*pop$tb_ms
  curr_pop <- c("t"=t,
                "tb"=tb,
                "tb_m"=tb_m,
                "tb_s"=tb_s,
                "tb_ms"=tb_ms,
                "sp_cure"=sp_cure,
                "tx_cure"=tx_cure,
                "died_tb"=died_tb,
                "died_nontb"=died_nontb,
                "rel_inf"=rel_inf)
  return(curr_pop)
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
#function to track TB deaths and relative infections
micro_outcomes <- function(p, pop) {
  rel_inf <- p$i*sum(pop==1) + p$i_m*sum(pop==2) + p$i_s*sum(pop==3) + p$i_ms*sum(pop==4)
  return(rel_inf)
}
#may be able to remove deaths from non-TB compartments to improve performance - unless we need to track # in spontaneously cured, etc.

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
                                                         outcome=="died_nontb"~"Non-TB Deaths",
                                                         outcome=="rel_inf"~"Rel. Secondary Infections"))
  sim_pop_long <- sim_pop_long %>% group_by(t) %>% 
    mutate(value_scale=case_when(outcome %in% alive_states~value/sum(value[outcome %in% alive_states]), #calculate as % of those still alive in time step
                                 outcome=="rel_inf"~value, #transform separately
                                 outcome %in% c("died_tb", "died_nontb")~value*n)) #scale by # ppl (already cumulative)
  sim_pop_long <- sim_pop_long %>% ungroup() %>% mutate(value_scale=if_else(outcome=="rel_inf", cumsum(value), value_scale))
  sim_pop_long <- sim_pop_long %>% group_by(t) %>% 
    mutate(prop_tb=if_else(outcome %in% tb_states, value/sum(value[outcome %in% tb_states]), 0))  
}
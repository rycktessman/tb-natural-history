#load packages, parameters, model code
library(tidyverse)
library(data.table)
source("code/model_functions.R")
load("data/params_all.Rda")
cyc_len <- 1/12
n_time <- 20 #to reach steady state, run for 20 years/240 timesteps
t_end <- 20/cyc_len
verbose <- 1

#convert params to transition matrix
trans_mat <- gen_trans_mat(params)

#generate vector with initial state distribution (t tracks timesteps)
start_pop <- c(1, #smear- symptom- TB
               0, #smear+ symptom- TB
               0, #smear- symptom+ TB
               0, #smear+ symptom+ TB
               0, #spontaneously cured
               0, #cured via diagnosis and treatment
               0, #TB death
               0 #non-TB death
)
sim_pop <- list()
sim_pop[[1]] <- t(matrix(start_pop))
for(t in 1:t_end) {
  print(t)
  curr_pop <- sim_pop[[t]]%*%trans_mat
  sim_pop[[t+1]] <- curr_pop
}
sim_pop <- do.call(rbind, sim_pop)
sim_pop <- cbind(0:t_end, sim_pop)
sim_pop <- as.data.frame(sim_pop)
names(sim_pop) <- c("t", "tb", "tb_m", "tb_s", "tb_ms", "sp_cure", "tx_cure", "died_tb", "died_nontb")

#clean and graph output
sim_pop_long <- process_output(sim_pop)
#plot markov trace for all states (this only makes sense if inflows are turned off)
ggplot(sim_pop_long %>% filter(outcome!="rel_inf"),
       aes(x=t, y=value*100, color=name)) + geom_line() +
  labs(x="Time", y="Percent in State", color="States") + theme_bw() + theme(panel.grid=element_blank())
#plot markov trace for living states
ggplot(sim_pop_long %>% filter(outcome %in% alive_states), 
       aes(x=t, y=value_scale*100, color=name)) + geom_line() +
  labs(x="Time", y="Percent in State", color="States") + theme_bw() + theme(panel.grid=element_blank())
#plot markov trace for TB states
ggplot(sim_pop_long %>% filter(outcome %in% tb_states), 
       aes(x=t, y=prop_tb*100, color=name)) + geom_line() +
  labs(x="Time", y="Percent in State", color="States") + theme_bw() + theme(panel.grid=element_blank())
#plot cumulative TB deaths 
ggplot(sim_pop_long %>% filter(outcome=="died_tb"), 
       aes(x=t, y=value_scale*100, color=name)) + geom_line() +
  labs(x="Time", y="Cumulative TB Deaths", color="") + theme_bw() + theme(panel.grid=element_blank())

#calculate deaths scaled by pop growth (so comparable with microsim)
inc_deaths <- sim_pop_long %>% filter(outcome=="died_tb" & t>0) %>% pull(value) -
  sim_pop_long %>% filter(outcome=="died_tb" & t<t_end) %>% pull(value)
pop <- sim_pop_long %>% filter(t<t_end & outcome %in% alive_states) %>% group_by(t) %>% summarise(pop=sum(value)) %>% pull(pop)
tb_deaths_scale <- inc_deaths/pop

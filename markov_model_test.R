#script to test the Markov model

setwd("~/GitHub/tb-natural-history")
library(tidyverse)

#load model functions and test model parameters (posterior means from main analysis)
source("model_functions.R")
load("params_test.Rda")

t_end <- 240 #run for 20 years (240 months)
inflows <- 1 #run with inflows on (deaths re-enter the model as smear-negative subclinical TB)
params <- c(params, "inflows"=inflows)

#declare starting compartment sizes
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

#run Markov model
sim_pop <- data.frame() 
sim_pop <- bind_rows(sim_pop, start_pop)
#run for t timesteps
verbose <- 1 #show model progress
for(t in 1:t_end) {
  if(verbose==1) {
    cat('\r', paste(round(t/t_end * 100), "% done", sep = " ")) # display the progress of the simulation
  }
  curr_pop <- nat_hist_markov(params, sim_pop[t,], t)
  sim_pop <- bind_rows(sim_pop, curr_pop)
}

#clean output
sim_pop_long <- process_output(sim_pop)

#graph Markov trace curves
ggplot(sim_pop_long %>% filter(outcome %in% alive_states), 
       aes(x=t, y=value_scale*100, color=name)) + geom_line() +
  labs(x="Time", y="Percent in State", color="States") + 
  theme_bw() + theme(panel.grid=element_blank())


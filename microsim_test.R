#script to test the microsimulation version of the model

setwd("~/GitHub/tb-natural-history")
library(tidyverse)

#load model functions and test model parameters (posterior means from main analysis)
source("model_functions.R")
load("params_test.Rda")

t_end <- 240 #run for 20 years (240 months)
inflows <- 1 #run with inflows on (deaths re-enter the model as smear-negative subclinical TB)
params <- c(params, "inflows"=inflows)
n <- 10000 #run with 1000 individuals

#columns in "sim_ind" are timesteps, rows are individuals
#cells in "sim_ind" indicate the state an individual is in at each timestep:
#1=tb (smear- symptom-)
#2=tb_m (smear+ symptom-)
#3=tb_s (smear- symptom+)
#4=tb_ms (smear+ symptom+)
#5=sp_cure (spontaneously cured)
#6=tx_cure (cured w/ diagnosis & treatment)
#7=died_tb (TB death)
#8=died_nontb (non-TB death)

#VERSION 1: everyone starts out smear-negative, run to steady state
rands <- micro_sample(params, n, t_end) #sample random numbers up front
sim_ind <- data.frame("0"=rep(1, n)) 
verbose <- 1 #show model progress
for(t in 1:t_end) {
  if(verbose==1) {
    cat('\r', paste(round(t/t_end * 100), "% done", sep = " ")) # display the progress of the simulation
  }
  curr_ind <- nat_hist_micro(sim_ind[,t], rands[(n*(t-1)+1):(n*t),], params) 
  sim_ind <- cbind(sim_ind, curr_ind) 
}
names(sim_ind) <- 0:t_end 

#process output
sim_ind_long <- data.frame(sapply(1:length(states), function(x) colSums(sim_ind==x))) #convert to dataframe of counts by state over time. rowSums if transposed
names(sim_ind_long) <- states
sim_ind_long <- cbind("t"=0:t_end, sim_ind_long/n) #convert to proportions
if(params$inflows==1) { #replace with cumulatives since deaths are reseeded if inflows turned on
  sim_ind_long$died_tb <- cumsum(sim_ind_long$died_tb) 
  sim_ind_long$died_nontb <- cumsum(sim_ind_long$died_nontb)
}
sim_ind_long <- process_output(sim_ind_long)

#graph Markov trace curves
ggplot(sim_ind_long %>% filter(outcome %in% alive_states), 
       aes(x=t, y=value_scale*100, color=name)) + geom_line() +
  labs(x="Time", y="Percent in State", color="States") + 
  theme_bw() + theme(panel.grid=element_blank())



#VERSION 2: simulate cohort w/ Philippines TB prevalence (turn off inflows)
props <- c(0.51, 0.20, 0.12, 0.17)
params[["inflows"]] <- 0
sim_ind <- data.frame("0"=c(rep(1, round(props[[1]]*n)),
                            rep(2, round(props[[2]]*n)),
                            rep(3, round(props[[3]]*n)),
                            rep(4, round(props[[4]]*n))))
verbose <- 1 #show model progress
for(t in 1:t_end) {
  if(verbose==1) {
    cat('\r', paste(round(t/t_end * 100), "% done", sep = " ")) # display the progress of the simulation
  }
  curr_ind <- nat_hist_micro(sim_ind[,t], rands[(n*(t-1)+1):(n*t),], params) 
  sim_ind <- cbind(sim_ind, curr_ind) 
}
names(sim_ind) <- 0:t_end 

#process output
sim_ind_long <- data.frame(sapply(1:length(states), function(x) colSums(sim_ind==x))) #convert to dataframe of counts by state over time. rowSums if transposed
names(sim_ind_long) <- states
sim_ind_long <- cbind("t"=0:t_end, sim_ind_long/n) #convert to proportions
if(params$inflows==1) { #replace with cumulatives since deaths are reseeded if inflows turned on
  sim_ind_long$died_tb <- cumsum(sim_ind_long$died_tb) 
  sim_ind_long$died_nontb <- cumsum(sim_ind_long$died_nontb)
}
sim_ind_long <- process_output(sim_ind_long)

#graph Markov trace curves
ggplot(sim_ind_long %>% filter(outcome %in% alive_states), 
       aes(x=t, y=value_scale*100, color=name)) + geom_line() +
  labs(x="Time", y="Percent in State", color="States") + 
  theme_bw() + theme(panel.grid=element_blank())

#calculate average duration with TB by starting state
durations <- sim_ind %>% mutate_all(~if_else(. %in% 1:4, 1 ,0))
durations <- rowSums(durations)
durations <- data.frame("start_pop"=sim_ind[,"0"], 
                        "duration"=durations)
durations %>% group_by(start_pop) %>%
  summarise(duration=mean(duration))

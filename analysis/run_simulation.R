setwd("~/GitHub/tb-natural-history")

#load packages
library(tidyverse)

#load model functions
source("R/model_functions.R")

#define simulation specifications
n <- 100000 #simulate 100,000 individuals
t_end <- 60 #run simulations for 5 years (60 months)
start_pop <- 1 #run for 1 cohort: smear- sublinical (1), smear+ subclinical (2), smear- symptomatic (3), smear+ symptomatic (4)
path_out <- "output/"

#load objects used in simulations (some overlap with calibration)
load("data/calibration_inputs.Rda")
#out_IMIS_combined is a file with all of the "out_IMIS_[chain].csv" files from calibration appended into 1 CSV file
params_post <- read.csv(paste0(path_out, "out_IMIS_combined.csv"), stringsAsFactors=F) 

#most efficient to split params_post into multiple blocks and run each separately/in parallel
index_start <- 1
index_end <- 10
split <- 1 #specify which split of 1-10 this is
params_use <- params_post[index_start:index_end,]
n_params <- nrow(params_use)

#parameter dependencies
params_use <- params_use %>% mutate(a_p_s=a_p_m, a_r_s=a_r_m, inflows=0, p_c=0, m_ac=params_fixed_prev$m_ac)

#create objects to store output
sim_ind_all <- data.frame()
sim_ind_each <- list()

#run microsim for each parameter set in params_use
for(i in 1:n_params) {
  print(i)
  params <- params_use[i,]
  param_id <- params$id
  rands <- micro_sample(params, n, t_end) #sample random numbers up front - same across param sets for each run
  sim_ind <- data.frame("0"=rep(start_pop, n)) #everyone starts out in the same state (start_pop)
  for(t in 1:t_end) {
    curr_ind <- nat_hist_micro(sim_ind[,t], rands[(n*(t-1)+1):(n*t),], params) 
    sim_ind <- cbind(sim_ind, curr_ind) 
  }
  names(sim_ind) <- 0:t_end 
  
  #process output
  sim_ind_long <- data.frame(sapply(1:length(states), function(x) colSums(sim_ind==x))) #convert to dataframe of counts by state over time
  names(sim_ind_long) <- states
  sim_ind_long <- cbind("t"=0:t_end, sim_ind_long/n) #convert to proportions
  if(params$inflows==1) { #replace with cumulatives since deaths are reseeded if inflows turned on
    sim_ind_long$died_tb <- cumsum(sim_ind_long$died_tb) 
    sim_ind_long$died_nontb <- cumsum(sim_ind_long$died_nontb)
  }
  sim_ind_long <- process_output(sim_ind_long)
  sim_ind_long <- sim_ind_long %>% mutate(id=param_id)
  sim_ind_all <- bind_rows(sim_ind_all, sim_ind_long)
  
  #save individual level output too
  list_index <- paste0(param_id)
  sim_ind_each[[list_index]] <- sim_ind
}

#save raw output to file
save(sim_ind_each, file=paste0(path_out, "microsim_output_pop", start_pop, "_", split, ".Rda"))
#remerge with posterior params 
sim_ind_all <- inner_join(params_post %>% select(id), sim_ind_all)
write.csv(sim_ind_all, file=paste0(path_out, "sim_ind_all_", start_pop, "_", split, ".csv"), row.names=F)

#describe time spent in each state and proportion who reach each state
sim_ind_each_comb <- bind_rows(sim_ind_each, .id="id")
times_all <- list()
prop <- data.frame()
for(j in unique(sim_ind_each_comb$id)) { #loop through each parameter set
  print(j)
  sim_ind <- data.matrix(sim_ind_each_comb %>% filter(id==j) %>% select(-id))
  
  #all states
  states <- 1:8
  out <- list()
  for(i in states) {
    out[[i]] <- rowSums(sim_ind==i)
  }
  #any smear+ state, any TB state
  out[["smear"]] <- rowSums(sim_ind==2|sim_ind==4)
  out[["tb_any"]] <- rowSums(sim_ind==1|sim_ind==2|
                               sim_ind==3|sim_ind==4)
  out[["symptom"]] <- rowSums(sim_ind==3|sim_ind==4)
  
  out <- bind_cols(out) 
  colnames(out)[states] <- states
  
  out_cond <- out
  out_cond[out_cond==0] <- NA #ignore 0 times spent in a state
  
  #bind lists into dataframes
  mean <- data.frame("id"=j, out %>% summarise_all(~mean(.)))
  med <- data.frame("id"=j, out %>% summarise_all(~median(.)))
  q25 <- data.frame("id"=j, unname(out %>% summarise_all(~quantile(., 0.25))))
  q75 <- data.frame("id"=j, unname(out %>% summarise_all(~quantile(., 0.75))))
  lb <- data.frame("id"=j, unname(out %>% summarise_all(~quantile(., 0.025))))
  ub <- data.frame("id"=j, unname(out %>% summarise_all(~quantile(., 0.975))))
  
  prop <- bind_rows(prop,data.frame("id"=j, out_cond %>% summarise_all(~sum(!is.na(.))/n())))
  
  if(j==unique(sim_ind_each_comb$id)[1]) {
    times_all[["mean"]] <- mean
    times_all[["med"]] <- med
    times_all[["q25"]] <- q25
    times_all[["q75"]] <- q75
    times_all[["lb"]] <- lb
    times_all[["ub"]] <- ub
    } else {
    times_all[["mean"]] <- bind_rows(times_all[["mean"]], mean)
    times_all[["med"]] <- bind_rows(times_all[["med"]], med)
    times_all[["q25"]] <- bind_rows(times_all[["q25"]], q25)
    times_all[["q75"]] <- bind_rows(times_all[["q75"]], q75)
    times_all[["lb"]] <- bind_rows(times_all[["lb"]], lb)
    times_all[["ub"]] <- bind_rows(times_all[["ub"]], ub)
  }
}
save(times_all, file=paste0(path_out, "times_all_", start_pop, "_", split, ".Rda"))
save(prop, file=paste0(path_out, "props_", start_pop, "_", split, ".Rda"))


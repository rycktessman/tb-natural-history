#run the microsim using posterior parameters to characterize stability of states

#load packages, parameters, model code
library(dplyr)
library(tidyr)
library(stringr)

source("code/model_functions.R")
source("code/calib_functions.R")

chain_split <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) #chain/file to open and rows to run
print(chain_split)
country <- as.character(Sys.getenv('country'))
print(country)

#define microsim specs
cyc_len <- 1/12 #weekly timestep
t_end <- 5/cyc_len #5 years
n <- 50000 

#microsim options to vary
RR_free <- 1 #4 free RR parameters in this version
spont_progress <- 0 #whether those who have spontaneously resolved can progress back to smear- symptom- TB
spont_prog <- 0.15 #what probability to use if spont_progress is 1
smear_hist_calib <- 0 #whether to include historical targets on bacillary status over time
deaths_targets <- "base" #"base", or "ihme" or use ihme targets
no_10yr_hist <- 0 #whether to include 10 year historical survival as calibration targets
smear_notif_override <- NA #NA, or an alt estimate +/- 10% (uniformly distributed)
start_pop <- as.numeric(Sys.getenv('start_pop')) #1=smear-/symptom-, 2=smear+/symptom-, 3=smear-/symptom+, 4=smear+/symptom+

#load files 
load("data/params_all.Rda")
load(paste0("data/targets_", tolower(country), ".Rda"))
path_out <- paste0("output/", tolower(country))

#file path
if(RR_free==0) {
  path_out <- paste0(path_out, "_rrconstrain")
} 
if(spont_progress==1) {
  path_out <- paste0(path_out, "_spontprog")
}
if(smear_hist_calib==1) {
  path_out <- paste0(path_out, "_smearhist")
}
if(deaths_targets=="ihme") {
  path_out <- paste0(path_out, "_ihmedeaths")
}
if(no_10yr_hist==1) {
  path_out <- paste0(path_out, "_no10")
}
if(!is.na(smear_notif_override)) {
  path_out <- paste0(path_out, "_smearnotif", as.character(round(smear_notif_override*100)))
}
if(RR_free==1 & spont_progress==0 & smear_hist_calib==0 & no_10yr_hist==0 & deaths_targets=="base" &
   is.na(smear_notif_override)) {
  path_out <- paste0(path_out, "_base")
}
path_out <- paste0(path_out, "/")
print(path_out)

#load params from file
params_post <- read.csv(paste0(path_out, "out_IMIS_combined", ".csv"), stringsAsFactors=F)
params_post_unique <- unique(params_post) #no need to evaluate same param set twice
size <- ceiling(nrow(params_post_unique)/1500) #split into 1500 parts
index_start <- (chain_split - 1)*size + 1
index_end <- pmin(index_start + size - 1, nrow(params_post_unique)) #last array has fewer rows
if(index_end < index_start) {
  print("No indices left to run - ending")
  print(chain_split)
} else {
  print(paste0("Running indices ", index_start, " through ", index_end))
  print(paste("Everyone starts out in state", start_pop))
  
  params_post_unique <- params_post_unique[index_start:index_end,]
  #parameter dependencies
  params_post_unique <- params_post_unique %>% 
    mutate(inflows=0, p_c=0, 
           m_ac=m_ac_present[[country]])
  
  #implement options
  if(RR_free==0) {
    params_post_unique <- params_post_unique %>% 
      mutate(a_p_s=a_p_m, a_r_s=a_r_m)
  }
  if(spont_progress==1) {
    params_post_unique <- params_post_unique %>% 
      mutate(p_c=1-exp(log(1-spont_prog)*cyc_len))
  }
  
  params_use <- params_post_unique
  
  #run the model runs times for t_end timesteps
  #rows are individuals, columns are timesteps, values are states - easy to transpose - see notes
  sim_ind_all <- data.frame()
  sim_ind_each <- list()
  n_params <- nrow(params_use)
  
  #RUN MICROSIM
  for(i in 1:n_params) {
    print(i)
    params <- params_use[i,]
    param_id <- params$id
    
    rands <- micro_sample(params, n, t_end) #sample random numbers up front - same across param sets for each run
    sim_ind <- data.frame("0"=rep(start_pop, n)) #everyone starts out w/ smear- symptom- TB - data.frame(rep(1,n)) if transposed
    rel_inf <- 0 #track relative infections separately
    for(t in 1:t_end) {
      curr_ind <- nat_hist_micro(sim_ind[,t], rands[(n*(t-1)+1):(n*t),], params) #use sim_pop[t,] if transposed
      sim_ind <- cbind(sim_ind, curr_ind) #rbind if transposed
      rel_inf <- c(rel_inf, micro_outcomes(params, curr_ind))
    }
    names(sim_ind) <- 0:t_end #row.names if transposed
    
    #process output
    sim_ind_long <- data.frame(sapply(1:length(states), function(x) colSums(sim_ind==x))) #convert to dataframe of counts by state over time. rowSums if transposed
    names(sim_ind_long) <- states
    sim_ind_long <- cbind("t"=0:t_end, sim_ind_long/n, "rel_inf"=rel_inf/n) #convert to proportions
    if(params$inflows==1) { #replace with cumulatives since deaths are reseeded if inflows turned on
      sim_ind_long$died_tb <- cumsum(sim_ind_long$died_tb) 
      sim_ind_long$died_nontb <- cumsum(sim_ind_long$died_nontb)
    }
    sim_ind_long <- process_output(sim_ind_long)
    sim_ind_long <- sim_ind_long %>% mutate(id=param_id)
    sim_ind_all <- bind_rows(sim_ind_all, sim_ind_long)
    
    #save individual level output too
    list_index <- param_id
    sim_ind_each[[list_index]] <- sim_ind
  }
  
  #SAVE RAW OUTPUT TO FILE
  save(sim_ind_each, file=paste0(path_out, "microsim_output_pop", start_pop, "_", chain_split, ".Rda"))
  #remerge with full posterior params (subset relevant indices) for appropriate duplications
  sim_ind_all <- inner_join(params_post %>% select(id), sim_ind_all)
  write.csv(sim_ind_all, file=paste0(path_out, "sim_ind_all_", start_pop, "_", chain_split, ".csv"), row.names=F)
  #use this file to generate Markov state transition figs
  
  #TIME SPENT in each state
  #also proportion ever reaching a state
  sim_ind_each_comb <- bind_rows(sim_ind_each, .id="id_run")
  print("rows bound")
  sim_ind_each_comb <- sim_ind_each_comb %>% 
    mutate(id=sapply(str_split(id_run, pattern="[.]"), '[', 1))
  print("id created")
  #sim_ind_each_post <- inner_join(out_post %>% select(id), sim_ind_each_comb, by="id")
  weight <- sum(params_post$id %in% unique(sim_ind_each_comb$id)) #combined weight for this array
  print(weight)
  weights <- table(params_post %>% filter(id %in% unique(sim_ind_each_comb$id)) %>% pull(id))
  print(weights)
  sim_ind_each_post <- left_join(sim_ind_each_comb, data.frame(weights), by=c("id"="Var1")) 
  print("weights merged")
  sim_ind_each_post <- sim_ind_each_post %>% mutate(weight=Freq/nrow(params_post)) %>%
    select(-Freq) 
  times_all <- list()
  prop <- data.frame()
  for(j in unique(sim_ind_each_post$id)) { #loop through each parameter set
    print(j)
    weight_id <- unique(sim_ind_each_post %>% filter(id==j) %>% pull(weight))
    sim_ind <- data.matrix(sim_ind_each_post %>% filter(id==j) %>%
                             select(-c(id, weight, id_run)))
    
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
    mean <- data.frame("id"=j, "weight"=weight_id, out %>% summarise_all(~mean(.)))
    med <- data.frame("id"=j, "weight"=weight_id, out %>% summarise_all(~median(.)))
    q25 <- data.frame("id"=j, "weight"=weight_id, 
                      unname(out %>% summarise_all(~quantile(., 0.25))))
    q75 <- data.frame("id"=j, "weight"=weight_id, 
                      unname(out %>% summarise_all(~quantile(., 0.75))))
    lb <- data.frame("id"=j, "weight"=weight_id, 
                     unname(out %>% summarise_all(~quantile(., 0.025))))
    ub <- data.frame("id"=j, "weight"=weight_id, 
                     unname(out %>% summarise_all(~quantile(., 0.975))))
    
    prop <- bind_rows(prop, 
                      data.frame("id"=j, "weight"=weight_id, 
                                 out_cond %>% summarise_all(~sum(!is.na(.))/n())))
    
    if(j==unique(sim_ind_each_post$id)[1]) {
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
  save(times_all, file=paste0(path_out, "times_all_", start_pop, "_", chain_split, ".Rda"))
  save(prop, file=paste0(path_out, "props_", start_pop, "_", chain_split, ".Rda"))
}


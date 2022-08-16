#run the microsim using posterior parameters to characterize stability of states

#load packages, parameters, model code
library(dplyr)
library(tidyr)
library(stringr)
source("code/model_v3.R")
source("code/calib_functions2.R")

#define microsim specs
runs <- 1 #number of runs per parameter set - individuals don't interact at all, so don't need multiple runs, just enough sample size
t_end <- 60 #5 years
n <- 50000 
verbose <- 0
chain_split <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) #chain/file to open and rows to run
print(chain_split)
mult_expand <- 1
RR_free <- 1
spont_progress <- 0
spont_prog <- 0.15 #annual probability of progression from resolved
country <- "Cambodia"
start_pop <- as.numeric(Sys.getenv('start_pop')) #1=smear-/symptom-, 2=smear+/symptom-, 3=smear-/symptom+, 4=smear+/symptom+

if(country=="Philippines") {
  load("data/params_targets.Rda")
} else if (country=="Vietnam") {
  load("data/params_targets_vietnam.Rda")
} else if (country=="Cambodia") {
  load("data/params_targets_cambodia.Rda")
} else if (country=="Nepal") {
  load("data/params_targets_nepal.Rda")
} else if (country=="Bangladesh") {
  load("data/params_targets_bangladesh.Rda")
}

#load params from file
path_out <- "output/IMIS Nov2021 v4 cambodia RRfree/"
params_post <- read.csv(paste0(path_out, "out_IMIS_combined", ".csv"), stringsAsFactors=F)
params_post_unique <- unique(params_post) #no need to evaluate same param set twice
size <- ceiling(nrow(params_post_unique)/1000) #split into 1000 parts
index_start <- (chain_split - 1)*size + 1
index_end <- pmin(index_start + size - 1, nrow(params_post_unique)) #last array has fewer rows
if(index_end < index_start) {
  print("no indices left to run - ending")
  print(chain_split)
} else {
  print(paste0("running indices ", index_start, " through ", index_end))
  print(paste("Everyone starts out in state", start_pop))
  params_post_unique <- params_post_unique[index_start:index_end,]
  #parameter dependencies
  params_post_unique <- params_post_unique %>% 
    mutate(inflows=0, p_c=0, 
           m_ac=params_fixed_prev$m_ac, i=1, i_m=1, i_s=2, i_ms=1)
  if(mult_expand==0) {
    params_post_unique <- params_post_unique %>% mutate(a_tx=a_m)
  }
  if(RR_free==0) {
    params_post_unique <- params_post_unique %>% 
      mutate(a_p_s=a_p_m, a_r_s=a_r_m)
  }
  if(spont_progress==1) {
    params_post_unique <- params_post_unique %>% 
      mutate(p_c=1-exp(log(1-spont_prog)/12))
  }
  params_use <- params_post_unique
  #run this to use MAP instead
  if(FALSE) {
    params_use <- unique(params_post_unique %>% filter(like==max(like)))
  }
  
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
    for(j in 1:runs) {
      print(j)
      rands <- micro_sample(params, n, t_end) #sample random numbers up front - same across param sets for each run
      sim_ind <- data.frame("0"=rep(start_pop, n)) #everyone starts out w/ smear- symptom- TB - data.frame(rep(1,n)) if transposed
      rel_inf <- 0 #track relative infections separately
      for(t in 1:t_end) {
        if(verbose==1) {
          cat('\r', paste(round(t/t_end * 100), "% done", sep = " ")) # display the progress of the simulation
        }
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
      sim_ind_long <- sim_ind_long %>% mutate(run=j, id=param_id)
      sim_ind_all <- bind_rows(sim_ind_all, sim_ind_long)
      
      #save individual level output too
      list_index <- paste0(param_id, ".", j)
      sim_ind_each[[list_index]] <- sim_ind
    }
  }
  
  #SAVE RAW OUTPUT TO FILE
  save(sim_ind_each, file=paste0(path_out, "microsim_output_pop", start_pop, "_", chain_split, ".Rda"))
  #remerge with full posterior params (subset relevant indices) for appropriate duplications
  sim_ind_all <- inner_join(params_post %>% select(id), sim_ind_all)
  write.csv(sim_ind_all, file=paste0(path_out, "sim_ind_all_", start_pop, "_", chain_split, ".csv"), row.names=F)
  #use this file to generate Markov state transition figs
  
  #DESCRIBE PROPORTION OF TIME IN EACH STATE
  #merge with out_post for sampling weights
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
    select(-Freq) #sum(Freq)
  #write.csv(sim_ind_each_post, file=paste0(path_out, "sim_ind_each_post", start_pop, "_",
   #                                        chain_split, ".csv"), row.names=F)
  #too many rows to hold in memory - sample 100,000 times from each array
  #indices <- sample(1:nrow(sim_ind_each_post), size=100000, replace=T, prob=sim_ind_each_post$weight)
  #sim_ind_each_post <- sim_ind_each_post[indices,] %>% select(-weight)
  
  #calculations
  
  #TIME SPENT (for everyone and conditionally) in each state
  #also proportion ever reaching a state
  times_all <- list()
  times_cond <- list()
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
    
    mean_cond <- data.frame("id"=j, "weight"=weight_id, out_cond %>% summarise_all(~mean(., na.rm=T)))
    med_cond <- data.frame("id"=j, "weight"=weight_id, out_cond %>% summarise_all(~median(., na.rm=T)))
    q25_cond <- data.frame("id"=j, "weight"=weight_id, 
                           unname(out_cond %>% summarise_all(~quantile(., 0.25, na.rm=T))))
    q75_cond <- data.frame("id"=j, "weight"=weight_id, 
                           unname(out_cond %>% summarise_all(~quantile(., 0.75, na.rm=T))))
    lb_cond <- data.frame("id"=j, "weight"=weight_id, 
                          unname(out_cond %>% summarise_all(~quantile(., 0.025, na.rm=T))))
    ub_cond <- data.frame("id"=j, "weight"=weight_id, 
                          unname(out_cond %>% summarise_all(~quantile(., 0.975, na.rm=T))))
    
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
      
      times_cond[["mean"]] <- mean_cond
      times_cond[["med"]] <- med_cond
      times_cond[["q25"]] <- q25_cond
      times_cond[["q75"]] <- q75_cond
      times_cond[["lb"]] <- lb_cond
      times_cond[["ub"]] <- ub_cond
      
    } else {
      times_all[["mean"]] <- bind_rows(times_all[["mean"]], mean)
      times_all[["med"]] <- bind_rows(times_all[["med"]], med)
      times_all[["q25"]] <- bind_rows(times_all[["q25"]], q25)
      times_all[["q75"]] <- bind_rows(times_all[["q75"]], q75)
      times_all[["lb"]] <- bind_rows(times_all[["lb"]], lb)
      times_all[["ub"]] <- bind_rows(times_all[["ub"]], ub)
      
      times_cond[["mean"]] <- bind_rows(times_cond[["mean"]], mean_cond)
      times_cond[["med"]] <- bind_rows(times_cond[["med"]], med_cond)
      times_cond[["q25"]] <- bind_rows(times_cond[["q25"]], q25_cond)
      times_cond[["q75"]] <- bind_rows(times_cond[["q75"]], q75_cond)
      times_cond[["lb"]] <- bind_rows(times_cond[["lb"]], lb_cond)
      times_cond[["ub"]] <- bind_rows(times_cond[["ub"]], ub_cond)
    }
  }
  save(times_all, file=paste0(path_out, "times_all_", start_pop, "_", chain_split, ".Rda"))
  save(times_cond, file=paste0(path_out, "times_cond_", start_pop, "_", chain_split, ".Rda"))
  save(prop, file=paste0(path_out, "props_", start_pop, "_", chain_split, ".Rda"))
  
  #DESCRIBE CONTINUOUS DURATION SPENT IN EACH STATE
  durations <- list()
  for(j in unique(sim_ind_each_post$id)) { #loop through each parameter set
    print(j)
    weight_id <- unique(sim_ind_each_post %>% filter(id==j) %>% pull(weight))
    sim_ind <- sim_ind_each_post %>% filter(id==j)
    
    mean <- list()
    med <- list()
    q25 <- list()
    q75 <- list()
    lb <- list()
    ub <- list()
    
    counters_all <- t(sapply(1:nrow(sim_ind), 
                             function(x) sequence(rle(as.character(sim_ind[x,]))$lengths)))
    for(i in 1:4)  { #each tb state
      print(i)
      #filter out rows where they never enter the state (improves performance)
      sim_sub <- sim_ind %>% filter(rowSums(sim_ind==i)!=0)
      #count how long someone has been in a given state
      if(i==start_pop) {
        counters <- counters_all
      } else {
        counters <- t(sapply(1:nrow(sim_sub), 
                             function(x) sequence(rle(as.character(sim_sub[x,]))$lengths)))
      }
      
      #isolate the state of interest
      counters_state <- counters*(sim_sub==i)
      #average over all episodes in a given state (not over rows)
      #find highest numbers that occur before a 0
      #append 0s to the end of each timeframe so that closing state counts toward this
      new_col <- rep(0, nrow(counters_state))
      counters_state <- cbind(counters_state, new_col)
      #convert to one long vector and just pull out end values in a sequence
      counters_state <- c(t(counters_state))
      counters_state <- counters_state[counters_state!=0 & lead(counters_state, 1)==0]
      
      mean[[as.character(i)]] <- mean(counters_state)
      med[[as.character(i)]] <- median(counters_state)
      q25[[as.character(i)]] <- quantile(counters_state, 0.25)
      q75[[as.character(i)]] <- quantile(counters_state, 0.75)
      lb[[as.character(i)]] <- quantile(counters_state, 0.025)
      ub[[as.character(i)]] <- quantile(counters_state, 0.975)
    }
    #smear+
    sim_sub <- sim_ind %>% filter(rowSums(sim_ind==2|sim_ind==4)!=0)
    sim_sub[sim_sub==2|sim_sub==4] <- 999 #new label for all smear-positive states
    #count how long someone has been in a given state
    counters <- t(sapply(1:nrow(sim_sub), function(x) sequence(rle(as.character(sim_sub[x,]))$lengths)))
    #isolate the state of interest
    counters_state <- counters*(sim_sub==999)
    #average over all episodes in a given state (not over rows)
    #find highest numbers that occur before a 0
    #append 0s to the end of each timeframe so that closing state counts toward this
    new_col <- rep(0, nrow(counters_state))
    counters_state <- cbind(counters_state, new_col)
    #convert to one long vector and just pull out end values in a sequence
    counters_state <- c(t(counters_state))
    counters_state <- counters_state[counters_state!=0 & lead(counters_state, 1)==0]
    mean[["smear"]] <- mean(counters_state)
    med[["smear"]] <- median(counters_state)
    q25[["smear"]] <- quantile(counters_state, 0.25)
    q75[["smear"]] <- quantile(counters_state, 0.75)
    lb[["smear"]] <- quantile(counters_state, 0.025)
    ub[["smear"]] <- quantile(counters_state, 0.975)
    
    #all TB
    sim_sub <- sim_ind
    sim_sub[sim_sub==1|sim_sub==2|sim_sub==3|sim_sub==4] <- 999 #new label for all smear-positive states
    #count how long someone has been in a given state
    counters <- t(sapply(1:nrow(sim_sub), function(x) sequence(rle(as.character(sim_sub[x,]))$lengths)))
    #isolate the state of interest
    counters_state <- counters*(sim_sub==999)
    #average over all episodes in a given state (not over rows)
    #find highest numbers that occur before a 0
    #append 0s to the end of each timeframe so that closing state counts toward this
    new_col <- rep(0, nrow(counters_state))
    counters_state <- cbind(counters_state, new_col)
    #convert to one long vector and just pull out end values in a sequence
    counters_state <- c(t(counters_state))
    counters_state <- counters_state[counters_state!=0 & lead(counters_state, 1)==0]
    mean[["tb"]] <- mean(counters_state)
    med[["tb"]] <- median(counters_state)
    q25[["tb"]] <- quantile(counters_state, 0.25)
    q75[["tb"]] <- quantile(counters_state, 0.75)
    lb[["tb"]] <- quantile(counters_state, 0.025)
    ub[["tb"]] <- quantile(counters_state, 0.975)
    
    #symptomatic
    sim_sub <- sim_ind %>% filter(rowSums(sim_ind==3|sim_ind==4)!=0)
    sim_sub[sim_sub==3|sim_sub==4] <- 999 #new label for all symptomatic states
    #count how long someone has been in a given state
    counters <- t(sapply(1:nrow(sim_sub), function(x) sequence(rle(as.character(sim_sub[x,]))$lengths)))
    #isolate the state of interest
    counters_state <- counters*(sim_sub==999)
    #average over all episodes in a given state (not over rows)
    #find highest numbers that occur before a 0
    #append 0s to the end of each timeframe so that closing state counts toward this
    new_col <- rep(0, nrow(counters_state))
    counters_state <- cbind(counters_state, new_col)
    #convert to one long vector and just pull out end values in a sequence
    counters_state <- c(t(counters_state))
    counters_state <- counters_state[counters_state!=0 & lead(counters_state, 1)==0]
    mean[["symptom"]] <- mean(counters_state)
    med[["symptom"]] <- median(counters_state)
    q25[["symptom"]] <- quantile(counters_state, 0.25)
    q75[["symptom"]] <- quantile(counters_state, 0.75)
    lb[["symptom"]] <- quantile(counters_state, 0.025)
    ub[["symptom"]] <- quantile(counters_state, 0.975)
    
    #bind lists into dataframes
    mean <- data.frame("id"=j, "weight"=weight_id, bind_cols(mean))
    med <- data.frame("id"=j, "weight"=weight_id, bind_cols(med))
    q25 <- data.frame("id"=j, "weight"=weight_id, bind_cols(q25))
    q75 <- data.frame("id"=j, "weight"=weight_id, bind_cols(q75))
    lb <- data.frame("id"=j, "weight"=weight_id, bind_cols(lb))
    ub <- data.frame("id"=j, "weight"=weight_id, bind_cols(ub))
    
    #add to durations list
    if(j==unique(sim_ind_each_post$id)[1]) {
      durations[["mean"]] <- mean
      durations[["med"]] <- med
      durations[["q25"]] <- q25
      durations[["q75"]] <- q75
      durations[["lb"]] <- lb
      durations[["ub"]] <- ub
      
    } else {
      durations[["mean"]] <- bind_rows(durations[["mean"]], mean)
      durations[["med"]] <- bind_rows(durations[["med"]], med)
      durations[["q25"]] <- bind_rows(durations[["q25"]], q25)
      durations[["q75"]] <- bind_rows(durations[["q75"]], q75)
      durations[["lb"]] <- bind_rows(durations[["lb"]], lb)
      durations[["ub"]] <- bind_rows(durations[["ub"]], ub)
      
    }
  }
  #save to file
  save(durations, file=paste0(path_out, "durations_", start_pop, "_", chain_split, ".Rda"))
}





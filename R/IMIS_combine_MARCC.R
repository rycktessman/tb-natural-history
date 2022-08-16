#combine multiple calibration runs - including output (update since v2) - from MARCC
library(dplyr)
library(matrixStats)
library(lhs)
#setwd("C:/Users/Tess/OneDrive - Johns Hopkins/TB/Natural History Modeling")
source("code/model_v3.R")
source("code/calib_functions2.R")
load("data/params_targets.Rda")

#info to be updated each set of calibration runs
path_out <- "output/IMIS Nov2021 v4 nepal smearhist/"

#find completed chains
completed <- list.files(path=path_out, pattern="IMIS_opt_each")
chains <- as.numeric(gsub(".*?([0-9]+).*", "\\1", completed))
chains <- c(1:14, 16:51)
chains <- sort(chains)
rounds <- 12 #number of IMIS rounds

#load output from IMIS
best_params_rounds <- data.frame() #optimal parameter set from each IMIS round within each chain
stats_rounds <- data.frame() #IMIS stats from each round within each chain
out_post_chains <- data.frame() #posterior parameter sets and corresponding output from each chain
for(i in chains) {
  print(i)
  if(i %in% c(as.character(9999))) {
    best_params <- read.csv(paste0(path_out, "IMIS_opt_each_", i, ".csv")) %>% 
      mutate(chain=i, round=1:11)
  } else {
    best_params <- read.csv(paste0(path_out, "IMIS_opt_each_", i, ".csv")) %>% 
      mutate(chain=i, round=1:rounds)
  }
  stats <- read.csv(paste0(path_out, "IMIS_stats_", i, ".csv")) %>% 
    mutate(chain=i, round=1:rounds)
  out_post <- read.csv(paste0(path_out, "out_IMIS_", i, ".csv")) %>%
    mutate(chain=i) %>% group_by_all() %>% 
    mutate(chain_id=cur_group_id(), id=paste0(chain, "_", chain_id)) #assign unique id to each param set
  best_params_rounds <- bind_rows(best_params_rounds, best_params)
  stats_rounds <- bind_rows(stats_rounds, stats)
  out_post_chains <- bind_rows(out_post_chains, out_post)
} 

#write to file
write.csv(out_post_chains, file=paste0(path_out, "out_IMIS_combined.csv"), row.names=F)
write.csv(stats_rounds, file=paste0(path_out, "stats_rounds_combined.csv"), row.names=F)
write.csv(best_params_rounds, file=paste0(path_out, "best_params_rounds_combined.csv"), row.names=F)

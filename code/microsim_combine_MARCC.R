#combine processed output from microsim_post_process_MARCC

#load packages, parameters, model code
library(dplyr)
library(tidyr)
library(data.table)

country <- as.character(Sys.getenv('country'))
analysis <- as.character(Sys.getenv('analysis'))
start_pop <- as.numeric(Sys.getenv('start_pop')) #1=smear-/symptom-, 2=smear+/symptom-, 3=smear-/symptom+, 4=smear+/symptom+
print(country)
print(analysis)
print(start_pop)

#update path
path_out <- paste0("output/", tolower(country), "_", analysis, "/")

rda2list <- function(file) {
  e <- new.env()
  load(file, envir=e)
  (as.list(e)[[1]])
}

stats <- c("mean", "med", "q25", "q75", "lb", "ub")

#find, load, combine files
setwd(path_out)

#1. time_all
completed <- list.files(pattern=paste0("times_all_", start_pop, "_"))
print(completed)
times_all_all <- Map(rda2list, completed)
times_all_comb <- list()
for(i in stats) {
  print(i)
  times_all_comb[[i]] <- bind_rows(lapply(1:length(completed), function(x) times_all_all[[x]][[i]]))
}
save(times_all_comb, file=paste0("times_all_comb", start_pop, ".Rda"))

#2. proportions
completed <- list.files(pattern=paste0("props_", start_pop, "_"))
print(completed)
props_all <- Map(rda2list, completed)
props_comb <- bind_rows(props_all)
write.csv(props_comb, file=paste0("props_comb", start_pop, ".csv"), row.names=F)

#3. state distribution for figures
completed <- list.files(pattern=paste0("sim_ind_all_", start_pop, "_"))
print(completed)
sim_ind_all_comb <- lapply(completed, fread)
save(sim_ind_all_comb, file=paste0("sim_ind_all_comb", start_pop, ".Rda"))

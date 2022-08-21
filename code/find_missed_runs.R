library(stringr) 

country <- as.character(Sys.getenv('country'))
analysis <- as.character(Sys.getenv('analysis'))
start_pop <- as.numeric(Sys.getenv('start_pop')) #1=smear-/symptom-, 2=smear+/symptom-, 3=smear-/symptom+, 4=smear+/symptom+
print(country)
print(analysis)
print(start_pop)

#update path
path_out <- paste0("output/", tolower(country), "_", analysis, "/")
print(path_out)

setwd(path_out)
completed_files <- list.files(pattern=paste0("props_", start_pop, "_"))
print(completed_files)
str(completed_files)
completed <- lapply(1:length(completed_files), function(x)
  str_split(completed_files[[x]], "_")[[1]][[3]])
completed <- lapply(1:length(completed), function(x)
  str_split(completed[[x]], "\\.")[[1]][[1]])
completed <- as.numeric(unlist(completed))
full <- 1:999
not_completed <- full[!(full %in% completed)]
out <- paste(not_completed, collapse=",")

fileConn <- file(paste0("missed_runs_", start_pop, ".txt"))
writeLines(out, fileConn)
close(fileConn)

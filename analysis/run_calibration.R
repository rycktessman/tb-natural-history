setwd("~/GitHub/tb-natural-history")

#load packages
library(lhs)
library(dplyr)
library(mvtnorm)
library(matrixStats)

#load objects used in calibration (parameter bounds, fixed parameters, and targets/related)
load("data/calibration_inputs.Rda")

#load model functions and functions used in calibration
source("R/model_functions.R")
source("R/calibration_functions.R")

cyc_len <- 1/12 #monthly time step
chain <- 1 #running 50 separate IMIS chains - which one is this?
path_out <- "output/"

#arguments needed for IMIS function - this will take about 3 days to run
B <- 10000 #10,000 samples per IMIS round (100,000 initially)
B.re <- 1000 #1000 samples of the posterior
number_k <- 12 #run 10 rounds of IMIS
D <- 0 #don't optimize first

#test using smaller number of samples and rounds (this takes about 30 minutes)
B <- 100
B.re <- 1000
number_k <- 5
D <- 0 #don't optimize first

post <- IMIS_copy(B, B.re, number_k, D)

post_params <- bind_cols(data.frame(post$resample), post$out)
write.csv(post_params, file=paste0(path_out, "out_IMIS_", chain, ".csv"), row.names=F)

params_opt <- post$center
write.csv(params_opt, file=paste0(path_out, "IMIS_opt_each_", chain, ".csv"), row.names=F)

stats_out <- post$stat
write.csv(stats_out, file=paste0(path_out, "IMIS_stats_", chain, ".csv"), row.names=F)

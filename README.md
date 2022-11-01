# tb-natural-history
TB natural history model

The files in this GitHub repository include data, R scripts, and bash scripts that are used to calibrate a model of TB natural, simulate individuals' TB trajectories, and estimate their contribution to TB transmission. 
A preprint of the paper is available at: https://doi.org/10.1101/2022.06.27.22276965
The repository is organized as follows:

**bash folder:**
The model code is written in R, but we used a high-performance computing cluster housed at a shared computing facility called MARCC (Maryland Advanced Research Computing Center) to run the analyses. 

This folder contains the bash scripts we used to submit jobs/run the model code on the server. They include (in the following order):

-run_IMIS_MARCC.sbatch: runs IMIS_run_MARCC.R; used to run calibration. Specifically, it runs 50 IMIS "chains" for a given country (specified in the bash script) and analysis (e.g., base case analysis or sensitivity analysis, specified in run_IMIS_MARCC.R)

-run_IMIS_combine_MARCC.sbatch: runs IMIS_combine_MARCC.R; combines calibration output from the 50 IMIS chains (from run_IMIS_MARCC.sbatch) for a given country (specified in the bash script) and analysis (specified in the bash script)

-run_microsim_MARCC.sbatch: runs microsim_run_MARCC.R; used to run individual simulations. Specifically, once the model is calibrated, for a given country (specified in the bash script), baseline TB state ("start_pop", specified in the bash script) and analysis (specified in microsim_run_MARCC.R), it simulates 50,000 individuals per posterior parameter set, all of whom start in that baseline TB state.

-run_microsim_combine_MARCC.sbatch : runs microsim_combine_MARCC.R; combines simulation output from run_microsim_MARCC.sbatch for a given country and analysis (both specified in the bash script), with output stratified by baseline TB state 

**code folder:**
Contains all the R scripts used to run the model. The approximate order to run them is:

-IMIS_run_MARCC.R: calibration script. Runs a single chain of IMIS for a specified country and analysis. This script has been written to run on the server, but manually entering a "chain" and "country" in lines 9-10 will allow it to be run on a desktop/similar. Lines 15-21 determine which analysis (e.g. main analysis or sensitivity analyses) will be run. 

-IMIS_combine_MARCC.R: calibration script. This script combines output from all 50 IMIS chains for a specified country and analysis. It has been written to run on the server.

-microsim_run_MARCC.R: simulation script. Splits the posterior parameters from calibration into manageable sizes (chunks) and runs 5-year individual simulations for a cohort of 50,000 people with TB all starting in the same baseline state for a given chunk, for a given country and analysis. This script has been written to run on the server but could be run on a desktop/similar by manually providing chain_split (specifies which chunk), country, and start_pop (baseline state) in lines 11, 13, and 15.

-microsim_combine_MARCC.R: simulation script. This script combines output from all 50,000 posterior parameters and 4 baseline TB states for a specified country and analysis. It has been written to run on the server.

-microsim_summary_relinf.R: simulation script. Calculates summary statistics from the simulations (output from microsim_combine_MARCC.R) and also conducts relative contribution to transmission calculations. It has been written to run on a desktop/equivalent (not the server).

Two scripts with functions that are used in the other scripts are also included:

-model_functions.R: function to create the transition matrix used in the Markov/state-transition form of the model and function to run the individual-level/microsimulation version of the model for 1 timestep. Also includes functions for processing output. 

-calib_functions.R: additional functions used in model calibration (IMIS_run_MARCC.R)

**data folder:**
Contains Rda files that are used to run the model. These include:

-targets_bangladesh.Rda, targets_cambodia.Rda, targets_nepal.Rda, targets_philippines.Rda, targets_vietnam.Rda: calibration targets for each country, including means and 95% CIs for each target, plus empirical distributions for the mortality-to-prevalence ratio target and the % of true TB notifications that are smear-positive (or would be if tested) target. 

-params_all.Rda: parameter info for all countries, including non-calibrated model parameters and lower/upper bounds for the uniform prior distributions for each calibrated parameter

-mort_to_prev_ihme.Rda: estimates of the mortality-to-prevalence target for all countries used in the sensitivity analysis in which we used IHME estimates of TB deaths (https://vizhub.healthdata.org/gbd-results/) instead of WHO estimates. 

**output folder:**
Contains the full sets of posterior parameters for each country from the main/base case calibration. 

Also includes a file with summary statistics from the main calibration (posterior_summary.csv)

The parameter labels correspond to the appendix model diagram in the paper (Figure S1)
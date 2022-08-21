#!/bin/bash

#SBATCH --output=R.%A.out
#SBATCH --error=R.%A.err 
#SBATCH --time=01:30:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tryckma1@jhmi.edu

######################
# Begin Work Section #
######################

#Run code
module purge
ml stack
echo "ml stack successful"
ml av
echo "ml av successful"
ml r/3.6.1
echo "ml r successful"

Rscript --vanilla code/gen_prior_samples_MARCC.R

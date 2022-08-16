#!/bin/bash

#SBATCH --output=R.%A_%a.out
#SBATCH --error=R.%A_%a.err 
#SBATCH --array=1-999
#SBATCH --time=4:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --mem=12000
##SBATCH --mem-per-cpu=4000
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

Rscript --vanilla code/microsim_posterior_MARCC4.R

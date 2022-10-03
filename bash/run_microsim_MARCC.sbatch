#!/bin/bash

#SBATCH --output=R.%A_%a.out
#SBATCH --error=R.%A_%a.err 
#SBATCH --array=1-1000
#SBATCH --time=2:00:00
#SBATCH -p defq
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tryckma1@jhmi.edu

######################
# Begin Work Section #
######################

#Run code
module purge
ml av
echo "ml av successful"
ml gcc/9.3.0
echo "ml gcc successful"
ml r/4.0.2
echo "ml r successful"

export country="Nepal"
export start_pop=1
Rscript --vanilla code/microsim_run_MARCC.R

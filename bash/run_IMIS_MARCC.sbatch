#!/bin/bash

#SBATCH --output=R.%A_%a.out
#SBATCH --error=R.%A_%a.err 
#SBATCH --array=1-50
#SBATCH --time=5:00:00
#SBATCH -p defq
#SBATCH --ntasks=1
##SBATCH --mem=10000
##SBATCH --mem-per-cpu=4000
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

export country="Nepal"
Rscript --vanilla code/IMIS_package_MARCC.R

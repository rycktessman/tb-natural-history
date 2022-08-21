#!/bin/bash

#SBATCH --output=R.%A_%a.out
#SBATCH --error=R.%A_%a.err 
#SBATCH --array=1,15,17,25,31,37,39,48,61,68,75,108,112,123,128,129,130,136,142,146,153,173,194,213,218,225,230,236,249,256,268,277,279,314,315,319,320,323,329,342,344,345,351,356,384,385,392,452,455,481,484,494,505,513,536,538,542,543,561,580,584,596,606,607,626,628,631,634,676,682,690,697,699,705,855,857,872,880,885,903,913,914,921,930,976,977,978,979,980,981,982,983,984,985,986,987,988,989,990,991,992,993,994,995,996,997,998,999
#SBATCH --time=4:00:00
#SBATCH -p defq
#SBATCH --ntasks=1
#SBATCH --mem=14000
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
echo "ml r successful"

export country="Nepal"
export start_pop=4
Rscript --vanilla code/microsim_posterior_MARCC.R

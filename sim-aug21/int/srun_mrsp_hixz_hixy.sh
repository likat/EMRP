#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=likat@umich.edu
#SBATCH --job-name=mmintcts
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
srun R-3.6 CMD BATCH --no-restore mmrp.R

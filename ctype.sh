#!/bin/bash
#SBATCH -J df_sml
#SBATCH --account=def-ubcxzh
#SBATCH --nodes=2
#SBATCH --cpus-per-task=20 
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-23:59

module load StdEnv/2020
module load r/4.0.2

export R_LIBS=~/R/library/4.0
Rscript main.R

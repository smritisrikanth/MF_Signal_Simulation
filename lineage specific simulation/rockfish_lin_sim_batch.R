#!/bin/bash

#SBATCH -J lin_spec_mut_sim
#SBATCH --partition=defq
#SBATCH -c 1                            # one CPU core per task
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript rockfish_lin_spec_simulation_script.R
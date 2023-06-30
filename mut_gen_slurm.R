#!/bin/bash

#SBATCH -J mut_gen_sim
#SBATCH --partition=defq
#SBATCH --array=1-100                    # how many tasks in the array
#SBATCH --cpus-per-task=4                          # one CPU core per task
#SBATCH -t 01:00:00
#SBATCH -o try-%j-%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript mut_gen_sim.R
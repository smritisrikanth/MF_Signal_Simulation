#!/bin/bash

#SBATCH -J three_cell_count_graph_sim
#SBATCH --partition=defq
#SBATCH --array=1-100                   # how many tasks in the array
#SBATCH --cpus-per-task=4                          # one CPU core per task
#SBATCH -t 02:00:00
#SBATCH -o a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript three_cell_count_graph.R

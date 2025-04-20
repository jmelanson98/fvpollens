#!/bin/bash
#SBATCH --job-name=dada2_workflow
#SBATCH --output=logs/dada2_%A_%a.out
#SBATCH --error=logs/dada2_%A_%a.err
#SBATCH --array=1-3                           
#SBATCH --ntasks=1                         
#SBATCH --cpus-per-task=4                  
#SBATCH --mem=640G                           
#SBATCH --time=06:00:00                  

# Set the R library path for the specific project
export R_LIBS_USER="/home/melanson/R/x86_64-pc-linux-gnu-library/4.3:/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/r/4.3.1/lib64/R/library"

# Load necessary modules (if needed, e.g., R version)
module load StdEnv/2023 r/4.3.1 python

# Run the R script
Rscript ~/projects/def-ckremen/melanson/fvpollens/4_bioinformatics/dada2_workflow.R


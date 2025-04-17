#!/bin/bash
#SBATCH --job-name=dada2_dereplication
#SBATCH --output=logs/dereplication_output.log
#SBATCH --array=1-3                           
#SBATCH --ntasks=1                         
#SBATCH --cpus-per-task=4                  
#SBATCH --mem=640G                           
#SBATCH --time=06:00:00                  

# Set the R library path for the specific project
export R_LIBS_USER="/home/melanson/R/x86_64-pc-linux-gnu-library/4.3:/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/r/4.3.1/lib64/R/library"

# Load necessary modules (if needed, e.g., R version)
module load StdEnv/2023 r/4.3.1 python/3.11.4

# Run the R script
Rscript ~/projects/def-ckremen/melanson/dada2_workflow.R


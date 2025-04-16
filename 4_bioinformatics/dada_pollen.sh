#!/bin/bash
#SBATCH --job-name=dada2_dereplication      # Job name
#SBATCH --output=dereplication_output.log    # Standard output and error log
#SBATCH --ntasks=1                           # Number of tasks (cores)
#SBATCH --cpus-per-task=4                   # Number of CPU cores per task (adjust based on your system)
#SBATCH --mem=640G                            # Memory per node (adjust as needed)
#SBATCH --time=06:00:00                      # Maximum run time (adjust as needed)

# Set the R library path for the specific project
export R_LIBS_USER="/home/melanson/R/x86_64-pc-linux-gnu-library/4.3:/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/r/4.3.1/lib64/R/library"

# Load necessary modules (if needed, e.g., R version)
module load StdEnv/2023 r/4.3.1 python/3.11.4

# Run the R script
Rscript ~/projects/def-ckremen/melanson/dada2_dereplicate_and_infer.R  # Replace with the path to your R script


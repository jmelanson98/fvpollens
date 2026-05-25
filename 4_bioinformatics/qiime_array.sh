#!/bin/bash
#SBATCH --job-name=qiime2_array
#SBATCH --output=logs/qiimearray_%A_%a.out
#SBATCH --error=logs/qiimearray_%A_%a.err
#SBATCH --array=1-3 
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00                

# Set the R library path for the specific project
export R_LIBS_USER="/home/melanson/R/x86_64-pc-linux-gnu-library/4.3:/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/r/4.3.1/lib64/R/library"

# Load necessary modules (if needed, e.g., R version)
module load StdEnv/2023 r/4.3.1 python

# Get task id
task_id=$(printf "%03d" $SLURM_ARRAY_TASK_ID)

# Command line code for Qiime2!
PROJECT_DIR=/project/6100170/melanson/def-ckremen/melanson/fvpollens

# make QIIME format
module load qiime2
qiime tools import \
--input-path ${PROJECT_DIR}/3_data/JeMe${task_id}_qiimeseqtab.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path ${PROJECT_DIR}/3_data/JeMe${task_id}_qiimeseqtab.qza

qiime tools import \
--input-path ${PROJECT_DIR}/3_data/JeMe${task_id}_ASVs.fasta \
--type 'FeatureData[Sequence]' \
--output-path ${PROJECT_DIR}/3_data/JeMe${task_id}_ASVs.qza
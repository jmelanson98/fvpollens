#!/bin/bash
#SBATCH --job-name=qiime_classify
#SBATCH --output=logs/qiimeclassify_%A_%a.out
#SBATCH --error=logs/qiimeclassify_%A_%a.err                         
#SBATCH --ntasks=1                         
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00               

# Set the R library path for the specific project
export R_LIBS_USER="/home/melanson/R/x86_64-pc-linux-gnu-library/4.3:/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/r/4.3.1/lib64/R/library"

# Load necessary modules (if needed, e.g., R version)
module load StdEnv/2023 r/4.3.1 python
PPROJECT_DIR=/project/6100170/melanson/fvpollens
SCRATCH_DIR=$SCRATCH/qiime_tmp
mkdir -p $SCRATCH_DIR

# Tell QIIME2 to use scratch for temp files
export TMPDIR=$SCRATCH_DIR
export MPLCONFIGDIR=$SCRATCH_DIR

module load qiime2


# Merge the 3 feature tables
qiime feature-table merge \
  --i-tables ${PROJECT_DIR}/3_data/JeMe001_qiimeseqtab.qza \
  --i-tables ${PROJECT_DIR}/3_data/JeMe002_qiimeseqtab.qza \
  --i-tables ${PROJECT_DIR}/3_data/JeMe003_qiimeseqtab.qza \
  --o-merged-table ${PROJECT_DIR}/3_data/merged_qiimeseqtab.qza

# Merge the ASV sequences
qiime feature-table merge-seqs \
  --i-data ${PROJECT_DIR}/3_data/JeMe001_ASVs.qza \
  --i-data ${PROJECT_DIR}/3_data/JeMe002_ASVs.qza \
  --i-data ${PROJECT_DIR}/3_data/JeMe003_ASVs.qza \
  --o-merged-data ${PROJECT_DIR}/3_data/merged_ASVs.qza

#get classifier from github
#wget https://github.com/apallavicini/PLANiTS/raw/master/qiime_classifiers/ITS2_classifier.qza

# also copy the classifier to scratch so it extracts there
cp ${PROJECT_DIR}/3_data/ITS2_classifier.qza $SCRATCH_DIR/

qiime feature-classifier classify-sklearn \
  --i-classifier ${SCRATCH_DIR}/ITS2_classifier.qza \
  --i-reads ${PROJECT_DIR}/3_data/merged_ASVs.qza \
  --o-classification ${PROJECT_DIR}/3_data/merged_taxonomy.qza \
  --p-n-jobs 8

# Export feature table to TSV
qiime tools export \
  --input-path ${PROJECT_DIR}/3_data/merged_qiimeseqtab.qza \
  --output-path ${PROJECT_DIR}/3_data/merged_exported

# Export taxonomy
qiime tools export \
  --input-path ${PROJECT_DIR}/3_data/merged_taxonomy.qza \
  --output-path ${PROJECT_DIR}/3_data/merged_exported_taxonomy

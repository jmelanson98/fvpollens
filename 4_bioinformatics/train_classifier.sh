#!/bin/bash
#SBATCH --job-name=train_classifier
#SBATCH --output=logs/trainclassifier_%A_%a.out
#SBATCH --error=logs/trainclassifier_%A_%a.err 
#SBATCH --job-name=train_classifier
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00


cd $SCRATCH
#wget https://github.com/apallavicini/PLANiTS/raw/master/PLANiTS_29-03-2020.zip
#unzip PLANiTS_29-03-2020.zip

export TMPDIR=$SCRATCH/tmp
mkdir -p $SCRATCH/tmp

module load qiime2

# Retrain classifier with updated QIIME2 version
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path $SCRATCH/ITS2.fasta \
  --output-path $SCRATCH/PLANiTS_ref_seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path $SCRATCH/ITS2_taxonomy \
  --input-format HeaderlessTSVTaxonomyFormat \
  --output-path $SCRATCH/PLANiTS_ref_tax.qza
  
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $SCRATCH/PLANiTS_ref_seqs.qza \
  --i-reference-taxonomy $SCRATCH/PLANiTS_ref_tax.qza \
  --o-classifier $SCRATCH/PLANiTS_classifier_2024.qza

# Copy classifier to working repo
cp $SCRATCH/PLANiTS_classifier_2024.qza /project/6100170/melanson/fvpollens/3_data/

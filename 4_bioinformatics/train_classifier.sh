#!/bin/bash
#SBATCH --job-name=train_keller_classifier
#SBATCH --output=logs/train_keller_%j.out
#SBATCH --error=logs/train_keller_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00

set -euo pipefail

mkdir -p logs
export TMPDIR=$SCRATCH/tmp
mkdir -p $SCRATCH/tmp

export INPUT_FA=$SCRATCH/its2.global.2023-01-17.curated.tax.mc.add.fa
export FASTA=$SCRATCH/keller_ITS2.fasta
export TAX=$SCRATCH/keller_ITS2_taxonomy.tsv

# --- Step 1: Parse combined FASTA into separate seqs + taxonomy files ---
module load StdEnv/2023 r/4.3.1 python

Rscript - <<'EOF'
input_fa  <- Sys.getenv("INPUT_FA")
out_fasta <- Sys.getenv("FASTA")
out_tax   <- Sys.getenv("TAX")

cat("Parsing", input_fa, "...\n")

con_in  <- file(input_fa, "r")
con_seq <- file(out_fasta, "w")
con_tax <- file(out_tax, "w")

n <- 0
while (TRUE) {
  line <- readLines(con_in, n = 1, warn = FALSE)
  if (length(line) == 0) break

  if (startsWith(line, ">")) {
    header      <- sub("^>", "", line)
    header      <- sub(";$", "", header)
    parts       <- strsplit(header, ";tax=")[[1]]
    seq_id      <- parts[1]
    tax_raw     <- parts[2]
    tax_levels  <- strsplit(tax_raw, ",")[[1]]
    tax_string  <- paste(sub(":", "__", tax_levels), collapse = ";")
    writeLines(paste0(">", seq_id), con_seq)
    writeLines(paste0(seq_id, "\t", tax_string), con_tax)
    n <- n + 1
  } else {
    writeLines(line, con_seq)
  }
}

close(con_in)
close(con_seq)
close(con_tax)
cat("Done. Parsed", n, "sequences.\n")
EOF

# --- Step 2: Train classifier ---
module load qiime2

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path $FASTA \
  --output-path $SCRATCH/keller_ref_seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path $TAX \
  --input-format HeaderlessTSVTaxonomyFormat \
  --output-path $SCRATCH/keller_ref_tax.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $SCRATCH/keller_ref_seqs.qza \
  --i-reference-taxonomy $SCRATCH/keller_ref_tax.qza \
  --o-classifier $SCRATCH/keller_classifier.qza \
  --p-n-jobs $SLURM_CPUS_PER_TASK

# --- Step 3: Copy to project space ---
cp $SCRATCH/keller_classifier.qza /project/6100170/melanson/fvpollens/3_data/

echo "All done!"
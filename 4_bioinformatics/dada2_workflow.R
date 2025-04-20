
setwd("/home/melanson/projects/def-ckremen/melanson/fvpollens")

library(dada2)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(ShortRead)
library(Biostrings)
library(seqinr)


####For Parallel Runs####
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))


####File Path Setup####
project_path <- sprintf("3_data/JeMe%03d", task_id)  # the directory containing the fastq files
scratch_path <- sprintf("/scratch/melanson/fvpollen_intermediates/JeMe%03d", task_id)

# list file names for forward and reverse reads
fnFs <- sort(list.files(project_path, pattern="_R1_001.fastq", full.names = TRUE)) 
fnRs <- sort(list.files(project_path, pattern="_R2_001.fastq", full.names = TRUE))

# split filenames to get unique sample name
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)


#### Filter and Trim Reads ####
# set output file paths
fnFs.filtN <- file.path(paste(scratch_path, "_filtN", sep = ""), basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(paste(scratch_path, "_filtN", sep = ""), basename(fnRs))

#filter and trim
print(getwd())
print(fnFs.filtN[1])
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, 
              trimLeft = c(0,0), 
              trimRight = c(0,0), 
              truncLen=c(250, 200), 
              maxN = 0, 
              multithread = 32, 
              compress = TRUE, 
              matchIDs=TRUE)



#### Prepare Functions for Primer Removal ####
FWD <- "ATGCGATACTTGGTGTGAAT"  ## ITS2 forward primer sequence
REV <- "TCCTCCGCTTATTGATATGC"  ## ITS2 reverse primer sequence

# create all orientations of the input sequence
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


# function for detecting primers in sequence data
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#index <- 6 #this is the index of the file we want to check for primers, within the lists "fn*s.filtN", it can be any number from 1 to N, where N is the number of samples you are processing
#rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[index]]), #the index of the sample you'd like to use for this test is used here (your first sample may be a blank/control and not have many sequences in it, be mindful of this)
#     FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[index]]), 
#     REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[index]]), 
#     REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[index]]))


#### Remove Primers ####
cutadapt <- "/home/melanson/projects/def-ckremen/melanson/fvpollens/4_bioinformatics/cutadapt-venv/bin/cutadapt" # cutadapt path on your machine
system2(cutadapt, args = "--version")

path.cut <- file.path(paste(scratch_path, "_cutadapt", sep = ""))
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

# run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-j", 32,# -n 2 required to remove FWD and REV from reads #-j sets no. threads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#sanity check, should report zero for all orientations and read sets
#index <- 6 #this is the index of the file we want to check for primers, within the lists "fn*s.cut", it can be any number from 1 to N, where N is the number of samples you are processing
#rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[index]]),
#      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[index]]),
#      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[index]]),
#      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[index]]))

# forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001", full.names = TRUE))



#### Filter and trim reads ... again! ####
#dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, 
#in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble 
#handling data with this many expected errors. it is best, after primer removal, to not truncate with 18s data, or with data from 
#any region in which the length is broadly variable. you may exclude organisms that have a shorter insert than the truncation length 
#(definitely possible, good example is giardia). defining a minimum sequence length is best

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(0,0), minLen = c(150,120),
                     maxN=c(0,0), maxEE=c(6,8), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=36)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
write.table(retained, 
            paste(sprintf("4_bioinformatics/dada2_output_files/JeMe%03d", task_id), "retained.reads.filterandtrim.txt", sep = ""), 
            sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)


#### Learn error rates ####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's
#sequence processing pipeline. read the dada2 paper and their online documentation

errF <- learnErrors(filtFs, multithread=32)
errR <- learnErrors(filtRs, multithread=32)

#assess this graph. it shows the error rates observed in your dataset. 
# strange or unexpected shapes in the plot should be considered before moving on.
pdf(paste(sprintf("4_bioinformatics/dada2_output_files/JeMe%03d", task_id), "error.rates.R1s.pdf"), width = 10, height = 10)
plotErrors(errF, nominalQ=TRUE) 
dev.off()

pdf(paste(sprintf("4_bioinformatics/dada2_output_files/JeMe%03d", task_id), "error.rates.R2s.pdf"), width = 10, height = 10)
plotErrors(errR, nominalQ=TRUE) 
dev.off()

#sometimes, samples with low sequence count don't produce any filtered output, 
#so we have to do something extra here to get rid of their entries in filtFs and filtRs
filtFs <- sort(list.files(paste(path.cut, "/filtered", sep = ""), pattern = "_R1_001.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(paste(path.cut, "/filtered", sep = ""), pattern = "_R2_001.fastq.gz", full.names = TRUE))

#now run dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# rewrite names of derepFs derepRs so that they match what's in sample.names 
# (just take the sampleID from the first element, split on _S)
names(derepFs) <- sapply(strsplit(names(derepFs), "_S"), `[`, 1)
names(derepRs) <- sapply(strsplit(names(derepRs), "_S"), `[`, 1)

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=32)
dadaRs <- dada(derepRs, err=errR, multithread=32)

#samples to keep
samples_to_keep <- as.numeric(out[,"reads.out"]) > 100
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples you have the option of removing
write.table(names(which(samples_to_keep == TRUE)), 
            paste(sprintf("4_bioinformatics/dada2_output_files/JeMe%03d", task_id), "samples.retained.txt", sep = ""), 
            row.names=FALSE, quote=F, sep="\n")
write.table(setdiff(sample.names, names(which(samples_to_keep == TRUE))),
            paste(sprintf("4_bioinformatics/dada2_output_files/JeMe%03d", task_id), "samples.removed.txt", sep = ""), 
            row.names=FALSE, quote=F, sep="\n")


#merge paired reads
mergers <- mergePairs(dadaFs[samples_to_keep], 
                      derepFs[samples_to_keep], 
                      dadaRs[samples_to_keep], 
                      derepRs[samples_to_keep], 
                      verbose=TRUE)

#make sequence table
seqtab <- makeSequenceTable(mergers)


#### View Sequence Length Distribution Post-Merging ####
# most useful with merged data. this plot will not show you much for forward reads only, 
# which should have a nearly uniform length distribution

#tabulate sequence length distribution
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab))))

pdf(paste(sprintf("4_bioinformatics/dada2_output_files/JeMe%03d", task_id), "merged.read.lengths.pdf", sep = ""), 
    width = 10, height = 8)
plot(x=length.histogram[,1], y=length.histogram[,2],
     xlab = "length (bp)",
     ylab = "frequency (# observed ASVs)",
     main = "ASV Length Histogram")
dev.off()


#### Remove low-count singleton ASVs ####
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

# generate counts of sample per ASV
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV

# create relative abundance table
otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x))

# convert to plain data frame
df <- as.data.frame(unclass(otus_rel_ab))

# set samples with no merged reads to 0!
df[is.na(df)] <- 0

# compute row sums (sum of relative abundance per ASV)
# for those present in only one sample, we can use this value to filter them for relative abundance on a per-sample basis
otus_rel_ab.rowsums <- rowSums(df)

# which ASVs are present in only one sample
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1)

# which OTUs are present below a relative abundance threshold?
b <- which(otus_rel_ab.rowsums <= 0.001)

# how many singleton ASVs fail on this relative abundance filter?
removed <- length(intersect(a,b))

# remove singletons with lower relative abundance than threshold
rows_to_remove <- intersect(a,b)

# filter these out of OTU table from earlier
if (removed > 0) {
otus_filt <- otus[-rows_to_remove,]
} else {
otus_filt <- otus #nothing to filter
}

#check how many ASVs you retained
dim(otus_filt)

# convert filtered OTU table back to sequence table
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt)))


#Start ASV report
ASV_report = paste(sprintf("4_bioinformatics/dada2_output_files/JeMe%03d", task_id), "ASVreport.txt", sep = "")

cat("Dimensions of unfiltered ASV table:", dim(otus), file = ASV_report, sep="\t", append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
cat("# singleton ASVs Removed:", length(intersect(a,b)), file= ASV_report, sep="\t", append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
cat("Dimensions of ASV table after singleton removal:", dim(otus_filt),file=ASV_report,sep="\t", append=TRUE)
cat("", file=ASV_report, sep="\n", append=TRUE)

#### Remove chimeras ####
#here we remove "bimeras" or chimeras with two sources. 
#look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=36,
 verbose=TRUE) 

# Append resuls to ASV table
cat("Dimensions of ASV table after chimera removal", dim(seqtab.nosingletons.nochim),file= ASV_report, sep="\t", append=TRUE)
cat("", file= ASV_report, sep="\n", append=TRUE)

# Calculate proportion of nonchimeras 
# should be relatively high after filtering out your singletons/low-count ASVs
# e.g., even if you lose a lot of ASVs, the number of reads lost should be quite low
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons)

# append to ASV table
cat("Proportion of non-chimeric to chimeric reads (0-1):", sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons), file= ASV_report, sep="\t", append=TRUE)

#### Track read retention through steps ####
#first remove samples in the "out" R object that aren't in samples_to_keep
row.names(out) <- as.character(t(as.data.frame(strsplit(row.names(out), "_S")))[,1])

getN <- function(x) sum(getUniques(x))

track <- cbind(out[which(row.names(out) %in% names(samples_to_keep)[which(samples_to_keep == TRUE)]),], 
               sapply(dadaFs[samples_to_keep], getN), 
               sapply(dadaRs[samples_to_keep], getN),
               sapply(mergers, getN),
               rowSums(seqtab.nosingletons),
               rowSums(seqtab.nosingletons.nochim))

track <- cbind(track, 
               100-track[,6]/track[,5]*100, 
               100-track[,7]/track[,6]*100, 
               track[,7]/track[,1]*100)
colnames(track) <- c("input", 
                     "filtered", 
                     "denoisedF", 
                     "denoisedR", 
                     "merged", 
                     "nosingletons", 
                     "nochimeras", 
                     "percent_singletons", 
                     "percent_chimeras", 
                     "percent_retained_of_total")


#### Save output from sequence table construction steps ####
write.table(data.frame("row_names"=rownames(track),track),
            paste(sprintf("4_bioinformatics/dada2_output_files/JeMe%03d", task_id), "/read.retention.merged.txt", sep = ""),
            row.names=FALSE, 
            quote=F, 
            sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim), seqtab.nosingletons.nochim),
            paste(sprintf("3_data/JeMe%03d", task_id), "_sequencetable.txt", sep = ""),
            row.names=FALSE, 
            quote=F, 
            sep="\t")


########################################
### ITS2 bioinformatics
### Jenna Melanson, Feb 2025
### based on script from Evan Morien
########################################

####Libraries####
library(dada2)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(ShortRead)
library(Biostrings)
library(seqinr)

####Environment Setup####
theme_set(theme_bw())
setwd("/Users/jenna1/Documents/UBC/Bombus Project/Pollen Metabarcoding")

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "JeMe001" # CHANGE ME to the directory containing the fastq files
list.files(path)

#change the pattern to match all your R1 files
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE)) 
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#change the delimiter in quotes and the number at the end of this command to 
#decide how to split up the file name, and which element to extract for a unique sample name
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)

#check your sample names
sample.names



####Quality Plots####
#if you have hundreds of samples in your dataset, it's easiest to just look at a few (10-50) fastq files
numplot <- 10
a <- sample(fnFs, numplot) #randomly select a set of N samples
b <- which(fnFs %in% a) #identify the indices of those samples
plotfnFs <- fnFs[b] #use the indices to create two lists of corresponding forward and reverse files
plotfnRs <- fnRs[b]

#plot quality profiles for each sample
pdf("quality_plots.dada2.R1s.pdf", width = 16, height = 9)
plotQualityProfile(plotfnFs)
dev.off()
pdf("quality_plots.dada2.R2s.pdf", width = 16, height = 9)
plotQualityProfile(plotfnRs)
dev.off()

####Primer Removal####
####here using bell primers####
FWD <- "ATGCGATACTTGGTGTGAAT"  ## CHANGE ME to your forward primer sequence
REV <- "TCCTCCGCTTATTGATATGC"  ## CHANGE ME to your reverse primer sequence

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

#filter and trim reads
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, trimLeft = c(0,0), trimRight = c(0,0), truncLen=c(250, 200), maxN = 0, multithread = 32, compress = TRUE, matchIDs=TRUE)

#create function for  
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
index <- 6 #this is the index of the file we want to check for primers, within the lists "fn*s.filtN", it can be any number from 1 to N, where N is the number of samples you are processing
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[index]]), #the index of the sample you'd like to use for this test is used here (your first sample may be a blank/control and not have many sequences in it, be mindful of this)
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[index]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[index]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[index]]))


#Primer removal
cutadapt <- "cutadapt-venv/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")

path.cut <- file.path("JeMe001_cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-j", 32,# -n 2 required to remove FWD and REV from reads #-j sets no. threads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#sanity check, should report zero for all orientations and read sets
index <- 6 #this is the index of the file we want to check for primers, within the lists "fn*s.cut", it can be any number from 1 to N, where N is the number of samples you are processing
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[index]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[index]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[index]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[index]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001", full.names = TRUE)) #remember to change this so it matches ALL your file names!
cutRs <- sort(list.files(path.cut, pattern = "_R2_001", full.names = TRUE)) #remember to change this so it matches ALL your file names!
  
  ####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, 
#in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble 
#handling data with this many expected errors. it is best, after primer removal, to not truncate with 18s data, or with data from 
#any region in which the length is broadly variable. you may exclude organisms that have a shorter insert than the truncation length 
#(definitely possible, good example is giardia). defining a minimum sequence length is best
.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(0,0), minLen = c(150,120),
                     maxN=c(0,0), maxEE=c(6,8), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=36)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
write.table(retained, "retained_reads.filterAndTrim_step.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's
#sequence processing pipeline. read the dada2 paper and their online documentation 
#(linked at top of this guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=32)
errR <- learnErrors(filtRs, multithread=32)

pdf("error_rates.dada2.R1s.pdf", width = 10, height = 10) # define plot width and height. completely up to user.
plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()
pdf("error_rates.dada2.R2s.pdf", width = 10, height = 10) # define plot width and height. completely up to user.
plotErrors(errR, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()

#sometimes, samples with low sequence count don't produce any filtered output, so we have to do something extra here to get rid of their entries in filtFs and filtRs
filtFs <- sort(list.files("JeMe001_cutadapt/filtered", pattern = "_R1_001.fastq.gz", full.names = TRUE)) #remember to change the pattern so it matches ALL your file names!
filtRs <- sort(list.files("JeMe001_cutadapt/filtered", pattern = "_R2_001.fastq.gz", full.names = TRUE)) #remember to change the pattern so it matches ALL your file names!

####dereplication####
#split filtFs in two because R runs out of memory otherwise!
filtFs_1 <- filtFs[seq(1, length(filtFs), by=3)]  
filtFs_2 <- filtFs[seq(2, length(filtFs), by=3)]
filtFs_2 <- filtFs[seq(3, length(filtFs), by=3)]

#dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# rewrite names of derepFs derepRs so that they match what's in sample.names (just take the sampleID from the first element, split on _)
names(derepFs) <- sapply(strsplit(names(derepFs), "_S"), `[`, 1)
names(derepRs) <- sapply(strsplit(names(derepRs), "_S"), `[`, 1)

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=32)
dadaRs <- dada(derepRs, err=errR, multithread=32)

dadaFs[[1]]
dadaRs[[1]]

####OPTIONAL: remove low-sequence samples before merging####
#a "subscript out of bounds" error at the next step (merging) may indicate that you aren't merging any reads in one or more samples.
#NB, NOT getting this error doesn't necessarily mean that all of your samples end up with more than 0 merged reads, as i found out 
#while processing a large 18s dataset. your guess is as good as mine as to why this error does or does not appear, but filtering 
#out the samples that cause it is necessary for completion of the pipeline.
samples_to_keep <- as.numeric(out[,"reads.out"]) > 100 #example of simple method used above after the
#filter and trim step. if you already did this but still got an error when merging, try the steps below


#getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
#track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
#samples_to_keep <- track[,4] > 1 #your threshold. try different ones to get the lowest one that will work. 
#this method accounts for dereplication/ASVs left after inference


samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples you have the option of removing
write.table(names(which(samples_to_keep == TRUE)), "samples_retained.txt", row.names=FALSE, quote=F, sep="\n")
write.table(setdiff(sample.names, names(which(samples_to_keep == TRUE))), "samples_removed.txt", row.names=FALSE, quote=F, sep="\n")

####merge paired reads####
#OPTION 2: modify command when removing low-sequence samples
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a nearly uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
pdf("length_histogram.merged_reads.pdf", width = 10, height = 8) # define plot width and height. completely up to user.
plot(x=length.histogram[,1], y=length.histogram[,2],
     xlab = "length (bp)",
     ylab = "frequency (# observed ASVs)",
     main = "ASV Length Histogram") #view length distribution plot
dev.off()

####remove low-count singleton ASVs####
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#generate counts of sample per ASV
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's 
#"transform_sample_counts" to create a relative abundance table, and then filter your ASVs by 
#choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)

otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step
#(a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need 
#to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those 
#only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
removed <- length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
if (removed > 0) {
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
} else {
otus_filt <- otus #nothing to filter
}
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

#Start ASV report
cat("dimensions of unfiltered ASV table:", dim(otus),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
cat("# singleton ASVs Removed:", length(intersect(a,b)),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
cat("dimensions of ASV table after singleton removal:", dim(otus_filt),file="ASV_report.txt",sep="\t",
    append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=36,
 verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
cat("dimensions of ASV table after chimera removal", dim(seqtab.nosingletons.nochim),file="ASV_report.
txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your 
#singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low
cat("proportion of non-chimeric to chimeric reads (0-1):", sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons),file="ASV_report.txt",sep="\t",append=TRUE)

####track read retention through steps####
#first remove samples in the "out" R object that aren't in samples_to_keep
row.names(out) <- as.character(t(as.data.frame(strsplit(row.names(out), "_S")))[,1]) 
#make sure that the output you're assigning to the row.names of "out" matches your sample names before proceeding. 
#select the correct delimiter (_ in exampe) and element (1 in example) to extract a unique name for each sample that corresponds to any metadata you will work with in downstream analyses
getN <- function(x) sum(getUniques(x))
track <- cbind(out[which(row.names(out) %in% names(samples_to_keep)[which(samples_to_keep == TRUE)]),]
, sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN),
 rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) withgetN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100, track[,7]/track[,1]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras", "percent_retained_of_total")

####save output from sequnce table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.16S_merged.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.ITS2.merged.txt", row.names=FALSE, quote=F, sep="\t")


#move into QIIME
#code modified from chatgpt
library(dada2)
library(biomformat)

# Load sequence table
seqtab <- read.delim("sequence_table.ITS2.merged.txt")
seqtab <- t(seqtab)#transpose the table
seqtab <- cbind('#OTUID' = rownames(seqtab), seqtab)#Add '#OTUID' to the header (required by biom)
write.table(seqtab, "seqtab.txt", sep='\t', row.names=FALSE, quote=FALSE)

#export representative sequences
asv_seqs <- DNAStringSet(colnames(seqtab)[-1])
names(asv_seqs) <- paste0("ASV", seq_along(asv_seqs))  # Assign ASV IDs

# Save as a FASTA file
writeXStringSet(asv_seqs, "rep_seqs.fasta")

####in CMD LINE NOT R
#make biom object
biom convert -i seqtab.txt -o seq_table.biom --table-type "OTU table" --to-hdf5

module load qiime2
qiime tools import \
--input-path /project/6100170/melanson/seq_table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path /project/6100170/melanson/seq_table.qza

qiime tools import \
--input-path rep_seqs.fasta \
--type 'FeatureData[Sequence]' \
--output-path rep_seqs.qza

#get classifier from github
#wget https://github.com/apallavicini/PLANiTS/blob/master/qiime_classifiers/ITS2_classifier.qza

qiime feature-classifier classify-sklearn \
--i-classifier plantITS_classifier.qza \
--i-reads rep_seqs.qza \
--o-classification taxonomy.qza


#to train classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref_seqs.qza \
--i-reference-taxonomy ref_tax.qza \
--o-classifier plantITS_classifier.qza

#convert asv to tsv
qiime tools export \
--input-path asv_table.qza \
--output-path exported-asv-table
biom convert -i exported-asv-table/feature-table.biom -o asv_table.tsv --to-tsv

#export tazonomy assignments
qiime tools export \
--input-path taxonomy.qza \
--output-path exported-taxonomy

#filter and visualize
qiime taxa barplot \
--i-table asv_table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file sample-metadata.tsv \
--o-visualization taxa-bar-plots.qzv




#view results
qiime tools view taxa-bar-plots.qzv




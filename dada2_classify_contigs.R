#!/usr/bin/env Rscript


# Name:		dada2_classify_contigs.R
# Author:	Tom van Wijk - RIVM Bilthoven
# Date:		12-12-2017
# Licence:	GNU General Public License v3.0 (copy provided in directory)

# This R script is a subscript of dada2_amplicon_pipeline.py
# For detailed information and instruction on how to install and use this software
# please view the included "README.md" file in this repository.

# The script can also be run separate from dada2_amplicon_pipeline.py
# with the following command, however, this is not recommended.

# Run with command: 'dada2_classify_contigs.R <input_dir> <reference_dir>
# input_dir:		FULL path! to input directory with .fastq files
# reference_dir:	FULL path! to durectory with reference files

# Requirements:
# - R v3.4.2 of later
# - Bioconductor v3.6 or later
# - dada 2 (R package)


# import R packages
library(dada2); packageVersion("dada2")

# Parse command line args
args = commandArgs(trailingOnly=TRUE)

# check if parameters are given
if (length(args) != 2) stop ("ERROR: Exactly 2 parameters must be given.")
setwd(args[1])
refpath = args[2]
list.files()

# get separate lists of forward and reverse reads
target_R1 <- c("_R1", ".fastq")
target_R2 <- c("_R2", ".fastq")
FR_files <- sort(list.files(pattern=target_R1, full.names = TRUE))
RR_files <- sort(list.files(pattern=target_R2, full.names = TRUE))

# check if no. of forward and reverse files are the same
if (length(FR_files) != length(RR_files)) stop ("ERROR: Number of forward and reverse files do not match.")

# check if forward and reverse files are properly paired by name
sample.FR_filenames <- sapply(strsplit(basename(FR_files), "_R1"), `[`, 1)
sample.RR_filenames <- sapply(strsplit(basename(RR_files), "_R2"), `[`, 1)
if (!identical(sample.FR_filenames, sample.RR_filenames)) stop ("ERROR: Forward and reverse files do not match.")

# adding names to forward and reverse files
names(FR_files) <- sample.FR_filenames
names(RR_files) <- sample.RR_filenames

# create quality rapport of raw input data
pdf(file = "fastq_quality_profiles_raw.pdf")
for (sample in sample.FR_filenames) {
       cat ("Generating quality profile:\t", sample, "\n")
       print(plotQualityProfile(list(FR_files[sample], RR_files[sample])))
}
dev.off()

# Create location for trimmed reads
FR_files_filt <- file.path(paste0("filtered/", sample.FR_filenames, "_R1_filt.fastq.gz"))
RR_files_filt <- file.path(paste0("filtered/", sample.RR_filenames, "_R2_filt.fastq.gz"))

# adding names to location of filtered forward and reverse files
names(FR_files_filt) <- sample.FR_filenames
names(RR_files_filt) <- sample.RR_filenames

# filtering and trimming the reads from FR_files to FR_files_filt
# WARNING: THESE PARAMETERS AREN'T OPTIMAL FOR ALL DATASETS, please inform "README.MD"
filterAndTrim(fwd=FR_files, filt=FR_files_filt,	rev=RR_files, filt.rev=RR_files_filt,
			  truncLen=c(235,185), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
			  compress=TRUE, verbose=TRUE, multithread=TRUE)

# create quality rapport of filtered and trimmed input data
pdf(file = "fastq_quality_profiles_filtered.pdf")
for (sample in sample.FR_filenames) {
       cat ("Generating quality profile:\t", sample, "\n")
       print(plotQualityProfile(list(FR_files_filt[sample], RR_files_filt[sample])))
}
dev.off()
	
set.seed(100)

# learn error rates of files
FR_err <- learnErrors(FR_files_filt, nread=1e6, multithread=TRUE)
RR_err <- learnErrors(RR_files_filt, nread=1e6, multithread=TRUE)

# create error plots of forward and reverse fastq files
pdf(file = "error_plot_FR.pdf")
plotErrors(FR_err, nominalQ=TRUE)
dev.off()
pdf(file = "error_plot_RR.pdf")
plotErrors(RR_err, nominalQ=TRUE)
dev.off()

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.FR_filenames))
names(mergers) <- sample.FR_filenames

# Merge samples
for (sample in sample.FR_filenames) {
	cat ("Processing:\t", sample, "\n")
	F_derep <- derepFastq(FR_files_filt[[sample]])
	F_dd <- dada(F_derep, err=FR_err, multithread=TRUE)
	R_derep <- derepFastq(RR_files_filt[[sample]])
	R_dd <- dada(R_derep, err=RR_err, multithread=TRUE)
	merger <- mergePairs(F_dd, F_derep, R_dd, R_derep)
	mergers[[sample]] <- merger
}

rm(F_derep)
rm(R_derep)

# Construct raw sequence table
seqtab <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove Chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
saveRDS(seqtab_nochim, "sequence_table.rds")
swapped_seqtab <- data.frame(t(seqtab_nochim))
write.table(swapped_seqtab, "sequence_table.txt", sep="\t", col.names=NA)

# Assign taxonomy
tax6 <- assignTaxonomy(seqtab_nochim, paste(refpath, "/silva_nr_v128_train_set.fa.gz", sep = ""), multithread=TRUE)
tax7 <- addSpecies(tax6, paste(refpath, "/silva_species_assignment_v128.fa.gz", sep = ""))
saveRDS(tax6, "taxonomy_level6_genus.rds")
saveRDS(tax7, "taxonomy_level7_species.rds")
write.table(tax6, "taxonomy_level6_genus.txt", sep="\t", col.names=NA)
write.table(tax7, "taxonomy_level7_species.txt", sep="\t", col.names=NA)

# check if number of entries in tax and seqtab files are the same
if (nrow(swapped_seqtab) != nrow(tax7)) stop ("ERROR: seqtab and taxonomy file to not match.")

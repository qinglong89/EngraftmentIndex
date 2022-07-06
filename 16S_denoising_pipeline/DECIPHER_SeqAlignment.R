#!/usr/bin/env Rscript

#usage: Rscript --vanilla DECIPHER_SeqAlignment.R /path-to-otu-seqs/
#alternative usage: R CMD BATCH IDTAXA.R (but need to specify the path directory in the script)

#define argumens from command line
args <- commandArgs(trailingOnly=TRUE)

#pass the work directory (where fastq files kept) from the command line to the script
#this only take the first argument from the command line
DirPath <- args[1]

#install DECIPHER which is designed to complement other tools that are part of BioConductor, a suite of packages for the R statistical programming language
#	install.packages("BiocManager")
#	BiocManager::install("DECIPHER")
	#get error: package ‘BiocManager’ is not available (for R version 3.4.4)
#this works for R version 3.4.4 in Ubuntu 18.0
#	source("https://bioconductor.org/biocLite.R")
#	biocLite("DECIPHER")

# load the DECIPHER library in R
library(DECIPHER)

# specify the path to the FASTA file (in quotes)
setwd(DirPath)
fas <- "Final_ASV-seq_DADA2.fasta"

# load the sequences from the file
# change "DNA" to "RNA" or "AA" for RNA and amino acid sequence alignment if necessary
seqs <- readDNAStringSet(fas)

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
#seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs, processors = 44)

# view the alignment in a browser (optional)
#BrowseSeqs(aligned, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned, file="Final_ASV-seq_DADA2_DECIPHIER-aligned.fasta")

#quit and save the R objects for future view
quit(save="no")

#!/usr/bin/env Rscript

#usage: Rscript --vanilla bash_DADA2+IDTAXA_workflow_Illumina-PairedEnd.R /path-to-fastq-files/
#alternative usage: R CMD BATCH bash_DADA2_script.R (but need to specify the path directory in the script)

#define argumens from command line
args <- commandArgs(trailingOnly=TRUE)

#pass the work directory (where fastq files kept) from the command line to the script
#this only take the first argument from the command line
DirPath <- args[1]

#alternative usage (need to deactivate above two commands)
#DirPath <- "/mnt/home1/qinglong/@16S-reanalysis/Illumina_data/Seattle_IBS/DADA2_process"


#this workflow is based on the tutorial of DADA2 version 1.8
library(dada2)
setwd(DirPath)
path <- DirPath
list.files(path)

###################################
###Inspect read quality profiles###
###################################

#TCMC forward and reverse fastq filenames have format: MS3894011-ID-v4-001-0714-TREFOIL-002-1_S1_L001_R1_001.fastq & MS3894011-ID-v4-001-0714-TREFOIL-002-1_S1_L001_R2_001.fastq 
fnFs <- sort(list.files(path, pattern=".fq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: MS3894011-ID-v4-001-0714-TREFOIL-002-1_S1_L001_R1_001.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

pdf("QualityProfile_forward.pdf")
plotQualityProfile(fnFs[1:2])
dev.off()

#pdf("QualityProfile_reverse.pdf")
#plotQualityProfile(fnRs[1:2])
#dev.off()

###################################
###########Filter and trim#########
###################################

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
#filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
#names(filtRs) <- sample.names

#quality filtering and trimming; this will also remove primers by position (trimLeft), filtered reads will have the length of (truncLen-trimLeft) 
#output is binary file, not traditional fastq file
out <- filterAndTrim(fnFs, filtFs, maxN=0, maxEE=2, rm.phix=TRUE, multithread=TRUE)

#this is to keep sequence sample with at least 20 or more reads for downstream, otherwise error "Error in outs[, 1] : subscript out of bounds" in mergePairs step
keep <- out[,"reads.out"] > 20
filtFs <- filtFs[keep]
#filtRs <- filtRs[keep]

###################################
#######Learn the Error Rates#######
###################################
errF <- learnErrors(filtFs, multithread=TRUE)
#errR <- learnErrors(filtRs, multithread=TRUE)

pdf("ErrorProfile_forward.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

#pdf("ErrorProfile_reverse.pdf")
#plotErrors(errR, nominalQ=TRUE)
#dev.off()

###################################
#########sample inference##########
###################################
#dereplicate sequence first
derepFs <- derepFastq(filtFs, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)

#name the derep-class objects by the sample names
names(derepFs) <- sample.names[keep]
#names(derepRs) <- sample.names[keep]

#apply the core sample inference algorithm to the dereplicated data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE) 
#dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

###################################
#########Merge paired reads########
###################################
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

###################################
####Construct sequence table#######
###################################
#seqtab <- makeSequenceTable(mergers)
seqtab <- makeSequenceTable(dadaFs)
#inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

###################################
#######Remove chimeras#############
###################################
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#check the rates of non-chimeras
sum(seqtab.nochim)/sum(seqtab)

###################################
##Track reads through the pipeline#
###################################
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
write.table(track, file = "Track_reads_through_DADA2_pipeline.txt", sep="\t")

###################################
##########Assign taxonomy##########
###################################
#this is default with RDP classifier
#taxa <- assignTaxonomy(seqtab.nochim, "/home/qinglong/16S_RefDB/DADA2_formatedDB/RefSeq-RDP16S_v2_May2018.fa.gz", multithread=TRUE)
#taxa <- addSpecies(taxa, "/home/qinglong/16S_RefDB/DADA2_formatedDB/RefSeq-RDP_dada2_assignment_species.fa.gz")

#suggest to use IDTAXA
library(DECIPHER)
#create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))
#load DECIPHIER pre-built training set
load("/home/qinglong/16S_RefDB/IDTAXA_TrainingSets/SILVA_SSU_r132_March2018.RData")
#load("/home/qinglong/16S_RefDB/IDTAXA_TrainingSets/RDP_v16-mod_March2018.RData")

ids <- IdTaxa(dna, trainingSet, strand="top", threshold=70, bootstraps=100, processors=20)

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

#convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

taxa <- taxid


###################################
########generate output############
###################################
ASVseq <- getSequences(seqtab.nochim)
write.table(ASVseq, file = "DADA2_ASV-seq.txt", sep="\t")
write.table(seqtab.nochim, file = "DADA2_ASV-feature-table.txt", sep="\t")
write.table(taxa, file = "DADA2_ASV-seq_taxa.txt", sep="\t")


#save all objects and class
quit(save="yes")


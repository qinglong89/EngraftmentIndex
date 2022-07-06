#!/bin/bash

WorkDir=/mnt/home1/qinglong/@Shotgun/@FMT_Kellermayer/Cdiff_isolates

THREADS=80

cd $WorkDir

for sample in $WorkDir/RawReadProcess/*
do

	cd $sample
	SampleName=`basename $sample`

	/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --fastq_filter \
                                                                  ''$SampleName'_R1.fastq.gz' \
                                                                  --fastqout ''$SampleName'_R1_maxEE1.fastq' \
                                                                  --fastq_truncee 1 --fastq_minlen 100 &

        /home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --fastq_filter \
                                                                  ''$SampleName'_R2.fastq.gz' \
                                                                  --fastqout ''$SampleName'_R2_maxEE1.fastq' \
                                                                  --fastq_truncee 1 --fastq_minlen 100

	wait

	#pairing the filtered forward and reverse reads
	/home/qinglong/softwares/bbmap-38.34/repair.sh in1=''$SampleName'_R1_maxEE1.fastq' \
        	                                       in2=''$SampleName'_R2_maxEE1.fastq' \
                                                       out1=''$SampleName'_R1_maxEE1_paired.fastq' \
                                                       out2=''$SampleName'_R2_maxEE1_paired.fastq' \
						       outs=''$SampleName'_maxEE1_unpaired.fastq'


	#de novo assembly
	#keep in mind that metaspades use k-mers of 21, 33 and 55, so "the length of reads should be not short than 55 nt" in theory.
	#only version 3.11.1 works for trimmed reads; version 3.12.0 and version 3.13.0 work for non-trimmed reads, which means the length of paired-end reads should be the same.
	/home/qinglong/softwares/SPAdes-3.11.1-Linux/bin/spades.py --only-assembler -t $THREADS -m 350 --sc \
 	                                                               -1 ''$SampleName'_R1_maxEE1_paired.fastq' \
 		                                                       -2 ''$SampleName'_R2_maxEE1_paired.fastq' \
	                                                               -s ''$SampleName'_maxEE1_unpaired.fastq' \
 	                                                               -o spades_assemblies


	#functional annotation
	source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/prokka
	cd spades_assemblies

	#rename contig headers (too long for prokka)
	/home/qinglong/softwares/bbmap-38.34/rename.sh in=scaffolds.fasta out=scaffolds_renamed.fasta prefix=contig

	prokka --cpus $THREADS -outdir prokka_annotation scaffolds_renamed.fasta
	source /home/DataAnalysis/miniconda2/bin/deactivate /home/qinglong/Programs/prokka

	rm $sample/*maxEE1*

done


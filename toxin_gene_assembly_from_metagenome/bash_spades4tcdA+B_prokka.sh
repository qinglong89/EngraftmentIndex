#!/bin/bash

WorkDir=/mnt/home1/qinglong/@Shotgun/@FMT_Kellermayer/

THREADS=66

cd $WorkDir

for sample in $WorkDir/RawReadProcess/*
do

	cd $sample
	SampleName=`basename $sample`

	#get reads mapped to custon database of tcdA and tcdB of many C. difficile isolates (NCBI)
	/home/qinglong/softwares/bowtie2-binary-2.3.4.3/bowtie2 -p $THREADS --very-sensitive \
		-x /mnt/home1/qinglong/Tutorials_Testings/Tutorial_Cdiff_toxin_genes_detection/Cdiff_toxin_genes/TcdAB_BT2-index \
		-1 ''$SampleName'_1.fq.gz' -2 ''$SampleName'_2.fq.gz' \
		--al-conc ''$SampleName'_TcdAB_paired'

	mv ''$SampleName'_TcdAB_paired.1' ''$SampleName'_TcdAB_paired.1.fq'
	mv ''$SampleName'_TcdAB_paired.2' ''$SampleName'_TcdAB_paired.2.fq'

	#run assemblies for tcdA and tcdB reads
	/home/qinglong/softwares/SPAdes-3.11.1-Linux/bin/spades.py \
		-1 ''$SampleName'_TcdAB_paired.1.fq' -2 ''$SampleName'_TcdAB_paired.2.fq' \
		--only-assembler -t $THREADS --sc -o TcdAB_assemblies

	#functional annotation
	source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/prokka

	cd TcdAB_assemblies

	/home/qinglong/softwares/bbmap-38.34/rename.sh in=scaffolds.fasta out=scaffolds_renamed.fasta prefix=contig

        prokka --cpus $THREADS -outdir prokka_annotation scaffolds_renamed.fasta

	source /home/DataAnalysis/miniconda2/bin/deactivate /home/qinglong/Programs/prokka

	rm $sample/*TcdAB*

done

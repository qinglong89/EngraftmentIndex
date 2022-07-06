#!/bin/bash

WorkDir=/mnt/home1/qinglong/@Shotgun/@FMT_Kellermayer/Cdiff_isolates/PediatricCdiffIsolates

THREADS=66

cd $WorkDir

mkdir -p Cdiff_ST11_snp


while read -r line
do
        cd $WorkDir/RawReadProcess/$line

	/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --fastq_filter ''$line'_1.fastq' \
								  --fastqout ''$line'_R1_maxEE1.fastq' \
								  --fastq_truncee 1 --fastq_minlen 100 &

	/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --fastq_filter ''$line'_2.fastq' \
                                                                  --fastqout ''$line'_R2_maxEE1.fastq' \
                                                                  --fastq_truncee 1 --fastq_minlen 100

	wait

	/home/qinglong/softwares/bbmap-38.34/repair.sh in1=''$line'_R1_maxEE1.fastq' \
						       in2=''$line'_R2_maxEE1.fastq' \
						       out1=''$line'_R1_maxEE1_paired.fastq' \
						       out2=''$line'_R2_maxEE1_paired.fastq' \
						       outs=''$line'_maxEE1_unpaired.fastq'

	source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/snippy

	snippy --cpus $THREADS --outdir snps_maxEE1_ref630 \
		--ref $WorkDir/CD630_NC009089.fa \
		--R1 ''$line'_R1_maxEE1_paired.fastq' \
		--R2 ''$line'_R2_maxEE1_paired.fastq'

	mv $WorkDir/RawReadProcess/$line/snps_maxEE1_ref630 $WorkDir/Cdiff_ST11_snp/''$line'_snps_maxEE1_ref630'

	source /home/DataAnalysis/miniconda2/bin/deactivate /home/qinglong/Programs/snippy


	rm *maxEE1*

done < $WorkDir/list_Cdiff_genomes_ST11.txt


#Core SNP phylogeny
cd $WorkDir/Cdiff_ST11_snp
snippy-core --prefix core ./*_snps_maxEE1_ref630 --ref $WorkDir/CD630_NC009089.fa
snippy-clean_full_aln core.full.aln > clean.full.aln
run_gubbins.py -p gubbins clean.full.aln
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
FastTree -gtr -nt clean.core.aln > clean.core.tree


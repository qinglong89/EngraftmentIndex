#!/bin/bash

WorkDir=/mnt/home1/qinglong/@16S-reanalysis/Illumina_data/V4_C26_IMTCDIFF
THREADS=44

mkdir $WorkDir/DADA2_process

for SampleDir in $WorkDir/RawReadProcess/*
do
	cd $SampleDir

	#"those variables might be changed if the input sequence is different"
	#for typical TCMC demultiplexed forward file: MS5413379-ID-v1v3-001-0817-TSCDIF2-UKCD1A-1_S1_L001_R1_001.fastq
	ForwardSeq=`ls *R1*`
	ReverseSeq=`ls *R2*`
	#SampleName=`ls *R1* | cut -d '_' -f 1`
	SampleName=`basename $SampleDir`

	#pairing the paired-end reads if they are disordered, sometimes possible in some centers
	/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --threads $THREADS \
								  --fastq_filter $ForwardSeq \
								  --fastq_trunclen 220 \
								  --fastqout ''$SampleName'_R1_trunc.fq'


        /home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --threads $THREADS \
                                                                  --fastq_filter $ReverseSeq \
                                                                  --fastq_trunclen 200 \
                                                                  --fastqout ''$SampleName'_R2_trunc.fq'

        /home/qinglong/softwares/bbmap-38.34/repair.sh in=''$SampleName'_R1_trunc.fq' \
                                                       in2=''$SampleName'_R2_trunc.fq' \
                                                       out=''$SampleName'_R1_matched.fq' \
                                                       out2=''$SampleName'_R2_matched.fq' \
                                                       outs=''$SampleName'_R1+R2_unmatched.fq'

	#merge the paired-end reads
	/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --threads $THREADS \
								  --fastq_mergepairs ''$SampleName'_R1_matched.fq' \
								  --reverse ''$SampleName'_R2_matched.fq' \
								  --fastqout ''$SampleName'_matched_merged.fq'
								  #--fastqout_notmerged_fwd ''$SampleName'_matched_unmerged_forward.fq'

	#quality filtering and trimming, and strip forward and reverse primers positionaly (pay attention to the order of parameters)
	/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --threads $THREADS \
								  --fastq_filter ''$SampleName'_matched_merged.fq' \
								  --fastqout ''$SampleName'_matched_merged_clean.fq' \
                                                                  --fastq_stripleft 31 \
                                                                  --fastq_stripright 32 \
								  --fastq_minlen 230 \
								  --fastq_maxlen 270

#        /home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --threads $THREADS \
#                                                                  --fastq_filter ''$SampleName'_matched_unmerged_forward.fq' \
#                                                                  --fastqout ''$SampleName'_matched_unmerged_forward_clean.fq' \
#                                                                  --fastq_stripleft 31 \
#                                                                  --fastq_minlen 200 \
#                                                                  --fastq_maxlen 270

#	cat ''$SampleName'_matched_merged_clean.fq' ''$SampleName'_matched_unmerged_forward_clean.fq' > ''$SampleName'_clean.fq'
#	mv ''$SampleName'_clean.fq' $WorkDir/DADA2_process

	mv ''$SampleName'_matched_merged_clean.fq' $WorkDir/DADA2_process

	rm *trunc* *matched*
done

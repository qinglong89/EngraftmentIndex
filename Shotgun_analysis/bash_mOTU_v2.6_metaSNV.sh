#!/bin/bash

WorkDir=/mnt/home1/@Untargeted_metagenomics/WGS_CDI/WGS_BGI_batch2/QW_SubAnalysis_FMT

MinLen=90

THREADS=80

cd $WorkDir

for folder in $WorkDir/RawReadProcess/*
do
	cd $folder
        SampleName=`basename $folder`

        R1=`ls *_1.fq.gz`    #forward read
        R2=`ls *_2.fq.gz`    #reverse read

        #perform expected accumulated error-based filtering
        /home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --threads $THREADS --fastq_qmax 100 \
                                                                  --fastq_filter $R1 \
                                                                  --fastqout ''$SampleName'_1_maxEE1.fq' \
                                                                  --fastq_truncee 1 \
                                                                  --fastq_stripleft 2 \
                                                                  --fastq_minlen $MinLen &

        /home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --threads $THREADS --fastq_qmax 100 \
                                                                  --fastq_filter $R2 \
                                                                  --fastqout ''$SampleName'_2_maxEE1.fq' \
                                                                  --fastq_truncee 1 \
                                                                  --fastq_stripleft 2 \
                                                                  --fastq_minlen $MinLen &

        wait

	#pairing the paired-end reads
	/home/qinglong/softwares/bbmap-38.34/repair.sh in=''$SampleName'_1_maxEE1.fq' \
						       in2=''$SampleName'_2_maxEE1.fq' \
						       out=''$SampleName'_R1_matched.fq' \
						       out2=''$SampleName'_R2_matched.fq' \
						       outs=''$SampleName'_R1+R2_unmatched.fq'


	python3 /home/qinglong/softwares/mOTUs_v2.6/motus map_snv -l $MinLen -t $THREADS \
								  -f ''$SampleName'_R1_matched.fq' \
								  -r ''$SampleName'_R2_matched.fq' \
								  -s ''$SampleName'_R1+R2_unmatched.fq' \
								  -o ''$SampleName'.bam'

        rm *maxEE* *matched*  #save storage
done

#merge individual profiles
mkdir $WorkDir/mOTU_v2.6_metaSNV_results && mkdir $WorkDir/mOTU_v2.6_metaSNV_results/sample-level_profile
mv $WorkDir/RawReadProcess/*/*.bam $WorkDir/mOTU_v2.6_metaSNV_results/sample-level_profile/

cd $WorkDir/mOTU_v2.6_metaSNV_results

#need to call metaSNV.py: only unique alignment was used for SNV profiling
source /home/DataAnalysis/miniconda2/bin/activate base
#python3 /home/qinglong/softwares/mOTUs_v2.6/motus snv_call -d $WorkDir/mOTU_v2.6_metaSNV_results/sample-level_profile \
#							   -t $THREADS -K \
#							   -o $WorkDir/mOTU_v2.6_metaSNV_results/SNV_output_default

#customized option based on per-sample per-species (mOTUs) vertical coverages (cov.tab; coverage depth: -fd) and horizontal coverage (perc.tab; coverage breadth: -fb) from above output
python3 /home/qinglong/softwares/mOTUs_v2.6/motus snv_call -d $WorkDir/mOTU_v2.6_metaSNV_results/sample-level_profile \
                                                           -fb 0.01 -fd 0.001 -fm 1 -fp 0.50 -fc 0.001 \
                                                           -t $THREADS -K \
                                                           -o $WorkDir/mOTU_v2.6_metaSNV_results/SNV_output_custom

source /home/DataAnalysis/miniconda2/bin/deactivate base

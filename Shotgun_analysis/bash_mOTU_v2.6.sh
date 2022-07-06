#!/bin/bash

WorkDir=/mnt/home1/@Untargeted_metagenomics/WGS_CDI/WGS_BGI_batch2/FMT-panphlan3

MinLen=90

THREADS=44

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
#	/home/qinglong/softwares/bbmap-38.34/repair.sh in=''$SampleName'_1_maxEE1.fq' \
#						       in2=''$SampleName'_2_maxEE1.fq' \
#						       out=''$SampleName'_R1_matched.fq' \
#						       out2=''$SampleName'_R2_matched.fq' \
#						       outs=''$SampleName'_R1+R2_unmatched.fq'

        #perform mOTU2 profiling: mOTU profile each read orientation separately, so not really necessary to have above pairing step
	python3 /home/qinglong/softwares/mOTUs_v2.6/motus profile -q -k mOTU -t $THREADS -n $SampleName \
								  -db /home/qinglong/softwares/mOTUs_v2.6/db_mOTU/ \
								  -f ''$SampleName'_1_maxEE1.fq' -r ''$SampleName'_2_maxEE1.fq' \
								  > ''$SampleName'_mOTU_v2.6_profile.txt'

        rm *maxEE*  #save storage
done

#merge individual profiles
mkdir $WorkDir/mOTU_v2.6_results && mkdir $WorkDir/mOTU_v2.6_results/sample-level_profile
mv $WorkDir/RawReadProcess/*/*_mOTU_v2.6_profile.txt $WorkDir/mOTU_v2.6_results/sample-level_profile/

cd $WorkDir/mOTU_v2.6_results/sample-level_profile
ProjectName=`basename $WorkDir`
BiomFileList=`ls -m | sed -e 's/, /,/g' | tr -d '\n'`
python3 /home/qinglong/softwares/mOTUs_v2.6/motus merge -i $BiomFileList -o $WorkDir/mOTU_v2.6_results/merged_mOTU_v2.6_profiles.txt


#parse
cd $WorkDir/mOTU_v2.6_results
cat merged_mOTU_v2.6_profiles.txt | sed '1,2d' | sed 's/#mOTU/#OTU ID/g' | sed 's/consensus_taxonomy/taxonomy/g' | cut -d$'\t' -f 1,3- > temp1
cat merged_mOTU_v2.6_profiles.txt | sed '1,2d' | sed 's/#mOTU/#OTU ID/g' | sed 's/consensus_taxonomy/taxonomy/g' | cut -d$'\t' -f 2 > temp2
paste temp1 temp2 > merged_mOTU_v2.6_profiles_parsed.txt && rm temp*

biom convert -i merged_mOTU_v2.6_profiles_parsed.txt -o merged_mOTU_v2.6_profiles_parsed.biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy
beta_diversity.py -i merged_mOTU_v2.6_profiles_parsed.biom -o beta-div_abund-jaccard -m abund_jaccard
principal_coordinates.py -i ./beta-div_abund-jaccard/* -o ./beta-div_abund-jaccard/pcoa_abund-jaccard.txt

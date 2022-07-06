#!/bin/bash

WorkDir=/mnt/home1/@Untargeted_metagenomics/WGS_CDI/WGS_BGI_batch2/QW_SubAnalysis_FMT

THREADS=44

MinLen=90

cd $WorkDir

mkdir $WorkDir/metaphlan3_results && mkdir $WorkDir/metaphlan3_results/sample-level_profile

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

        #perform metaphlan2
        cat ''$SampleName'_1_maxEE1.fq' ''$SampleName'_2_maxEE1.fq' > ''$SampleName'_combined_maxEE1.fastq'


	source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/biobakery3

	metaphlan ''$SampleName'_combined_maxEE1.fastq' \
		--input_type fastq --nproc $THREADS \
		--ignore_eukaryotes --ignore_archaea \
		> ''$SampleName'_metaphlan3_profile.txt'

	source /home/DataAnalysis/miniconda2/bin/deactivate /home/qinglong/Programs/biobakery3

	mv ''$SampleName'_metaphlan3_profile.txt' $WorkDir/metaphlan3_results/sample-level_profile

        rm *_maxEE* *bowtie2out.txt  #save storage
done

source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/biobakery3
merge_metaphlan_tables.py $WorkDir/metaphlan3_results/sample-level_profile/*_metaphlan3_profile.txt > $WorkDir/metaphlan3_results/merged_metaphlan3_profiles.txt

cd $WorkDir/metaphlan3_results/

#extract species profile only
grep -E "(s__)|(^clade_name)" merged_metaphlan3_profiles.txt | sed 's/clade_name/#OTU ID/g' | cut -f 1,3- > merged_metaphlan3_profiles_species.txt

biom convert -i merged_metaphlan3_profiles_species.txt -o merged_metaphlan3_profiles_species.biom --table-type="OTU table" --to-hdf5
beta_diversity.py -i merged_metaphlan3_profiles_species.biom  -o beta-div_abund-jaccard -m abund_jaccard
principal_coordinates.py -i ./beta-div_abund-jaccard/* -o ./beta-div_abund-jaccard/pcoa_abund-jaccard.txt

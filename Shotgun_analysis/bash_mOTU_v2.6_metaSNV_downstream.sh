#!/bin/bash

WorkDir=/mnt/home1/@Untargeted_metagenomics/WGS_CDI/WGS_BGI_batch2/QW_SubAnalysis_FMT

mkdir $WorkDir/mOTU_v2.6_metaSNV_results/temp

cd $WorkDir/mOTU_v2.6_metaSNV_results/SNV_output_custom/filtered-*

for file in $PWD/*.freq
do
        FileName=`basename $file .filtered.freq`
        sed '1d' $file > temp1  #remove first row
        sed -n '1p' $file | cut -f 2- | sed 's/.bam//g' > temp2   #extract first row
        echo "#OTU ID" > temp3
        paste temp3 temp2 > temp4
        cat temp4 temp1 > ''$FileName'_parsed.txt'
        biom convert -i ''$FileName'_parsed.txt' -o ''$FileName'_parsed.biom' --table-type="OTU table" --to-hdf5
        mv ''$FileName'_parsed.biom' $WorkDir/mOTU_v2.6_metaSNV_results/temp/
        rm temp* ''$FileName'_parsed.txt'
done

cd $WorkDir/mOTU_v2.6_metaSNV_results/temp/
BiomFileList=$(ls -m | sed -e 's/, /,/g' | tr -d '\n')
merge_otu_tables.py -i $BiomFileList -o ../mOTU_v2.6_metaSNV_filtered-freq_AllSamples.biom
rm -rf $WorkDir/mOTU_v2.6_metaSNV_results/temp


cd $WorkDir/mOTU_v2.6_metaSNV_results
beta_diversity.py -i mOTU_v2.6_metaSNV_filtered-freq_AllSamples.biom -o beta-div_binary-jaccard -m binary_jaccard
principal_coordinates.py -i ./beta-div_binary-jaccard/* -o ./beta-div_binary-jaccard/pcoa_binary-jaccard.txt


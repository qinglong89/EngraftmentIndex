#!/bin/bash

WorkDir=/mnt/home1/@Untargeted_metagenomics/WGS_CDI/WGS_BGI_batch2/FMT-panphlan3/panphlan3_results

#process absence/presence profile
cd $WorkDir

mkdir $WorkDir/FakeOTU_detection

for file in $WorkDir/*.tsv
do
        FileName=`basename $file .tsv`
        sed '1d' $file > temp1  #remove first row
        sed -n '1p' $file | cut -f 2- > temp2   #extract first row
        echo "#OTU ID" > temp3
        paste temp3 temp2 > temp4
        cat temp4 temp1 > ''$FileName'_parsed.txt'
        biom convert -i ''$FileName'_parsed.txt' -o ''$FileName'_parsed.biom' --table-type="OTU table" --to-hdf5
        mv ''$FileName'_parsed.biom' $WorkDir/FakeOTU_detection
        rm temp* ''$FileName'_parsed.txt'
done

cd $WorkDir/FakeOTU_detection
BiomFileList=$(ls -m | sed -e 's/, /,/g' | tr -d '\n')
merge_otu_tables.py -i $BiomFileList -o Panphlan3_AllSpecies.biom

beta_diversity.py -i Panphlan3_AllSpecies.biom -o beta-div_binary-jaccard -m binary_jaccard
principal_coordinates.py -i ./beta-div_binary-jaccard/* -o ./beta-div_binary-jaccard/pcoa_binary-jaccard.txt


#process coverage profile
cd $WorkDir

mkdir $WorkDir/FakeOTU_coverage

for file in $WorkDir/*.coverage
do
        FileName=`basename $file .coverage`
        sed '1d' $file > temp1  #remove first row
        sed -n '1p' $file | cut -f 2- > temp2   #extract first row
        echo "#OTU ID" > temp3
        paste temp3 temp2 > temp4
        cat temp4 temp1 > ''$FileName'_parsed.txt'
        biom convert -i ''$FileName'_parsed.txt' -o ''$FileName'_parsed.biom' --table-type="OTU table" --to-hdf5
        mv ''$FileName'_parsed.biom' $WorkDir/FakeOTU_coverage
        rm temp* ''$FileName'_parsed.txt'
done

cd $WorkDir/FakeOTU_coverage
BiomFileList=$(ls -m | sed -e 's/, /,/g' | tr -d '\n')
merge_otu_tables.py -i $BiomFileList -o Panphlan3_AllSpecies.biom

beta_diversity.py -i Panphlan3_AllSpecies.biom -o beta-div_abund-jaccard -m abund_jaccard
principal_coordinates.py -i ./beta-div_abund-jaccard/* -o ./beta-div_abund-jaccard/pcoa_abund-jaccard.txt

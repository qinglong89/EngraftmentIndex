#!/bin/bash

WorkDir=/home1/qinglong/phylophlan3_Cdiff

THREADS=66

cd $WorkDir

source /home1/qinglong/Programs/miniconda2/bin/activate phylophlan3

#one-time set-up
#phylophlan_setup_database -g s__Clostridioides_difficile -o phylophlan3_Cdiff-DB --verbose 2>&1 | tee logs/phylophlan_setup_database.log
#phylophlan_write_config_file -o phylophlan3_config.cfg -d a --force_nucleotides --db_aa diamond --map_aa diamond --map_dna diamond --msa mafft --trim trimal --tree1 fasttree --tree2 raxml

#add RefSeq C. difficile genomes
scp -r qinglong@10.26.24.62:/mnt/home1/qinglong/@Shotgun/@FMT_Kellermayer/Cdiff_isolates/CustomAlignment/Cdiff_CompleteGenomes_RefSeq $WorkDir
mkdir $WorkDir/input_genomes
cp ./Cdiff_CompleteGenomes_RefSeq/refseq//bacteria/*/*.fna.gz $WorkDir/input_genomes/
for file in $WorkDir/input_genomes/*; do FileName=`basename $file | cut -d '_' -f 1-2` && mv $file ./input_genomes/'RefSeq_'$FileName'.fna.gz';done

#add selected C. difficile genomes sharing the same alteration to tcdA by IS component
scp -r qinglong@10.26.24.62:/mnt/home1/qinglong/@Shotgun/@FMT_Kellermayer/Cdiff_isolates/CustomAlignment/Cdiff_SelectedGenomes_P6-1_tcdA-like $WorkDir
for genome in $WorkDir/Cdiff_SelectedGenomes_P6-1_tcdA-like/*; do GenomeName=`basename $genome` && cp $genome/*_genomic.fna.gz ./ && mv *_genomic.fna.gz $WorkDir/input_genomes/'P6-1_tcdA-like_'$GenomeName'.fna.gz';done

#add TCH and BTH C. difficile isolates

phylophlan -i $WorkDir/input_genomes -o phylophlan3_out_Cdiff \
	-d phylophlan3_Cdiff-DB -f phylophlan3_config.cfg \
	--trim greedy --remove_fragmentary_entries -t a \
	--fragmentary_threshold 0.67 --not_variant_threshold 0.99 --min_num_entries 408 \
	--diversity low --force_nucleotides \
	--nproc $THREADS --verbose 2>&1 | tee logs/phylophlan__output_isolates.log


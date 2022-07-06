#!/bin/bash

WorkDir=/mnt/home1/qinglong/@16S-reanalysis/Illumina_data/V4_C26_IMTCDIFF

THREADS=44

cd $WorkDir/DADA2_process

#mkdir DADA2_process

#cp $WorkDir/RawReadProcess/*/* $WorkDir/DADA2_process

######################################################################################################################################
#choose the type of processing, "you will need to change some parameters on the corresponding R file"

#for 16S V4 (paired-end reads)
#Rscript --vanilla /mnt/home1/qinglong/@16S-reanalysis/Illumina_data/V4_C26_IMTCDIFF/bash_DADA2+IDTAXA_workflow_Illumina-PairedEnd_V4.R $WorkDir/DADA2_process

#for 16S V1V3 (paired-end reads)
#Rscript --vanilla /home/qinglong/Programs/DADA2_pipeline/bash_DADA2+IDTAXA_workflow_Illumina-PairedEnd_V1V3.R $WorkDir/DADA2_process

#for single end reads
Rscript --vanilla /mnt/home1/qinglong/@16S-reanalysis/Illumina_data/V4_C26_IMTCDIFF/bash_DADA2+IDTAXA_workflow_single-end_fixed-length.R $WorkDir/DADA2_process
########################################################################################################################################

cd $WorkDir/DADA2_process

rm -rf filtered *.fastq* *.fq*

###parse ASV feature table
awk 'NR==1 {print; exit}' DADA2_ASV-feature-table.txt > first-line.txt
cat DADA2_ASV-feature-table.txt | sed 1d > rest-lines.txt

echo "#OTU ID" > file.txt
paste -d"\t" file.txt first-line.txt > first-line_parsed.txt

#make sure you have "datamash" installed
cat first-line_parsed.txt rest-lines.txt | sed 's/"//g' | datamash transpose > DADA2_ASV-feature-table_parsed.txt

rm first-line.txt rest-lines.txt file.txt first-line_parsed.txt

###parse taxonomy file
cat DADA2_ASV-seq_taxa.txt | sed 's/"//g' | sed 1d | cut -f2- | tr "\\t" ";" > taxon_parsed.txt
cat DADA2_ASV-seq_taxa.txt | sed 's/"//g' | sed 1d | cut -f 1 > ASV_ids.txt

paste -d"\t" ASV_ids.txt taxon_parsed.txt > DADA2_ASV-seq_taxa_parsed.txt

rm taxon_parsed.txt ASV_ids.txt

#extract ASV ids that are not classified at kingdom level (second column/field of the file)
cat DADA2_ASV-seq_taxa.txt | sed 's/"//g' | awk '$2 == "NA" { print $0 }' | cut -f 1 > ASV_ids_unclassified_kingdom.txt

###convert to biom file, ASV sequences as OTU identifiers
source /home/DataAnalysis/miniconda2/bin/activate qiime191

biom convert -i DADA2_ASV-feature-table_parsed.txt -o DADA2_ASV-feature-table_parsed.biom --to-hdf5 --table-type="OTU table"

biom add-metadata -i DADA2_ASV-feature-table_parsed.biom -o DADA2_ASV-feature-table_parsed_taxa.biom \
                  --observation-metadata-fp DADA2_ASV-seq_taxa_parsed.txt \
                  --observation-header OTUID,taxonomy
#call QIIME version 1.91 from conda
filter_otus_from_otu_table.py -i DADA2_ASV-feature-table_parsed_taxa.biom -o DADA2_ASV-feature-table_parsed_taxa_nohost.biom -e ASV_ids_unclassified_kingdom.txt
biom convert -i DADA2_ASV-feature-table_parsed_taxa_nohost.biom -o DADA2_ASV-feature-table_parsed_taxa_nohost.txt --to-tsv --header-key taxonomy

rm ASV_ids_unclassified_kingdom.txt


###parse ASV-based OTU table and prepare ASV sequences for phylogeny
cat DADA2_ASV-feature-table_parsed_taxa_nohost.txt | sed 1d | sed 1d | cut -f 2- > temp1	#remove first two lines and remove first column
line=`cat temp1 | wc -l` && seq $line > temp2 && sed -e 's/^/ASV/' temp2 > temp3 && paste temp3 temp1 > temp4	#create tab file containing ASV1, ASV2, ASV3 rather than DNA sequence as identifier
awk 'NR==1,NR==2' DADA2_ASV-feature-table_parsed_taxa_nohost.txt > temp5	#extract first two lines
cat temp5 temp4 | sed 's/; NA//g' > Final_ASV-table_DADA2.txt	#must remove NA in the taxonomy, othwise will cause problem in STAMP analysis
biom convert -i Final_ASV-table_DADA2.txt -o Final_ASV-table_DADA2.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
biom summarize-table -i Final_ASV-table_DADA2.biom -o Final_ASV-table_DADA2_stats.txt	#get sequencing depth per sample

cat DADA2_ASV-feature-table_parsed_taxa_nohost.txt | sed 1d | sed 1d | cut -f 1 > temp6	#extract the DNA sequence of each ASVs
paste -d"\t" temp3 temp6 | awk '{print ">"$1"\n"$2}' > Final_ASV-seq_DADA2.fasta	#create tab file containing the sequences for ASV1, ASV2, ASV3..., and then convert to fasta file

rm temp* DADA2_ASV*

#align ASV sequence with DECIPHIER
Rscript --vanilla /home/qinglong/Programs/DADA2_pipeline/DECIPHER_SeqAlignment.R $WorkDir/DADA2_process

#construct phylogenetic tree
make_phylogeny.py -i Final_ASV-seq_DADA2_DECIPHIER-aligned.fasta -t fasttree -o Final_ASV-seq_DADA2_DECIPHIER-aligned.tre

#prepare file for STAMP analysis
/home/qinglong/Programs/DADA2_pipeline/biom_to_stamp.py -m taxonomy Final_ASV-table_DADA2.biom > Final_ASV-table_DADA2_STAMP.spf

source /home/DataAnalysis/miniconda2/bin/deactivate qiime191

#####################################################################################################################################################################
#"optional" - perform picrust2 analysis

cd $WorkDir/DADA2_process

mkdir picrust2_out_pipeline && cd picrust2_out_pipeline

#do not use this for miniconda2 (won't be activated in sub-shell): conda activate picrust2
source /home/DataAnalysis/miniconda2/bin/activate picrust2

#Place reads into reference tree
place_seqs.py -s ../Final_ASV-seq_DADA2.fasta -o out.tre -p $THREADS --intermediate intermediate/place_seqs

#Hidden-state prediction of gene families
hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p $THREADS -n
hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p $THREADS

#Generate metagenome predictions
metagenome_pipeline.py -i ../Final_ASV-table_DADA2.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz \
                       -o EC_metagenome_out --strat_out

convert_table.py EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
                 -c contrib_to_legacy \
                 -o EC_metagenome_out/pred_metagenome_contrib.legacy.tsv.gz

#Pathway-level inference
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
                    -o pathways_out -p $THREADS

#Add functional descriptions
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz

#If you would like to try out STAMP for visualizing the PICRUSt2 output, see https://github.com/picrust/picrust2/wiki/STAMP-example

#do not use this for miniconda2 (won't be activated in sub-shell): conda deactivate

source /home/DataAnalysis/miniconda2/bin/deactivate picrust2

#end of process

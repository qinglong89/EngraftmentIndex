#!/bin/bash

WorkDir=/mnt/home1/qinglong/Tutorial_Cdiff_toxin_genes_detection/
THREADS=10

cd $WorkDir

##################################################################################################################################
####################################"constructing database of toxin genes of C. difficile"########################################
##################################################################################################################################

#download the genomes (genbank format) of C. difficile from NCBI GenBank database
#follow instructions in github: https://github.com/kblin/ncbi-genome-download
ncbi-genome-download --section genbank --format genbank --assembly-level all \
			--genus "Clostridioides difficile" --species-taxid 1496 \
			--parallel $THREADS --no-cache \
			-o $WorkDir/NCBI_CDiff \
			--metadata-table $WorkDir/NCBI_CDiff_metadata.txt \
			bacteria

#extract toxin genes and binary toxin genes from genbank file of each C. difficile strain
for StrainFolder in $WorkDir/NCBI_CDiff/genbank/bacteria/GCA*
do
	cd $StrainFolder
	#StrainName=`basename $StrainFolder`  #might need to consider the header identifier of each gene for building BT2 index
	pigz -d *_genomic.gbff.gz
	python $WorkDir/parseGenbanksExtractCDSsByGeneName.py -c tcdA -g *_genomic.gbff
	python $WorkDir/parseGenbanksExtractCDSsByGeneName.py -c tcdB -g *_genomic.gbff
	python $WorkDir/parseGenbanksExtractCDSsByGeneName.py -c tcdC -g *_genomic.gbff
	python $WorkDir/parseGenbanksExtractCDSsByGeneName.py -c tcdE -g *_genomic.gbff
	python $WorkDir/parseGenbanksExtractCDSsByGeneName.py -c tcdR -g *_genomic.gbff
	python $WorkDir/parseGenbanksExtractCDSsByGeneName.py -c cdtA -g *_genomic.gbff
	python $WorkDir/parseGenbanksExtractCDSsByGeneName.py -c cdtB -g *_genomic.gbff
	python $WorkDir/parseGenbanksExtractCDSsByGeneName.py -c cdtR -g *_genomic.gbff
done

#concatenate nucleotide sequences from all strains of C. difficile for each toxin gene
cd $WorkDir
mkdir Cdiff_toxin_genes
cat $WorkDir/NCBI_CDiff/genbank/bacteria/GCA*/tcdA.fasta > $WorkDir/Cdiff_toxin_genes/all_Cdiff_tcdA_DNA.fasta &
cat $WorkDir/NCBI_CDiff/genbank/bacteria/GCA*/tcdB.fasta > $WorkDir/Cdiff_toxin_genes/all_Cdiff_tcdB_DNA.fasta &
cat $WorkDir/NCBI_CDiff/genbank/bacteria/GCA*/tcdC.fasta > $WorkDir/Cdiff_toxin_genes/all_Cdiff_tcdC_DNA.fasta &
cat $WorkDir/NCBI_CDiff/genbank/bacteria/GCA*/tcdE.fasta > $WorkDir/Cdiff_toxin_genes/all_Cdiff_tcdE_DNA.fasta &
cat $WorkDir/NCBI_CDiff/genbank/bacteria/GCA*/tcdR.fasta > $WorkDir/Cdiff_toxin_genes/all_Cdiff_tcdR_DNA.fasta &
cat $WorkDir/NCBI_CDiff/genbank/bacteria/GCA*/cdtA.fasta > $WorkDir/Cdiff_toxin_genes/all_Cdiff_cdtA_DNA.fasta &
cat $WorkDir/NCBI_CDiff/genbank/bacteria/GCA*/cdtB.fasta > $WorkDir/Cdiff_toxin_genes/all_Cdiff_cdtB_DNA.fasta &
cat $WorkDir/NCBI_CDiff/genbank/bacteria/GCA*/cdtR.fasta > $WorkDir/Cdiff_toxin_genes/all_Cdiff_cdtR_DNA.fasta &
wait

##################################################################################################################################
################################"manually inspect and manipulate of toxin gene sequences"###########################################
##################################################################################################################################
cd $WorkDir/Cdiff_toxin_genes

#"do the same for TcdB gene"
##################################################################################################################################
#"manipulate the sequence identifier to replace space with underline"
cat all_Cdiff_tcdA_DNA.fasta | perl -pe 's/>*( )+/_/g' > all_Cdiff_tcdA_DNA_header-modified.fasta
grep '^>' all_Cdiff_tcdA_DNA_header-modified.fasta | sed 's/^>//g' > all_Cdiff_tcdA_DNA_header-modified_identifiers.txt

#"with nano, manually inspect the gene names are correct and save to _checked.txt file" and "filter gene names that are not correct"
#samtools faidx does not take tab-delimited file list, even convert to one line with space-delimited; "for small, do manually; for big, use seqtk subseq"
xargs samtools faidx all_Cdiff_tcdA_DNA_header-modified.fasta < all_Cdiff_tcdA_DNA_header-modified_identifiers_checked.txt \
			> all_Cdiff_tcdA_DNA_header-modified_filtered.fasta

#"make the header identifier clean and identical", and convert multiple-lined fasta sequence to single-line
/home/qinglong/softwares/fastx_toolkit_0.0.13_binaries/fasta_formatter -w 0 \
	-i all_Cdiff_tcdA_DNA_header-modified_filtered.fasta \
	| sed 's/^>/>TcdA-/g' \
	| cut -f1 -d '.' > TcdA.fasta
##################################################################################################################################

#"combine the sequences of toxin A and toxin B since TcdA and TcdB share 63% homology in their amino acid sequences"
cat TcdA.fasta TcdA.fasta > TcdAB.fasta

#make sure you have local database creadted before running
/home/qinglong/softwares/ncbi-blast-binary-2.7.1/bin/makeblastdb -in /mnt/home1/qinglong/Toxin_genes_detection-test/Cdiff_toxin_genes/TcdAB.fasta \
                                               -out /mnt/home1/qinglong/Toxin_genes_detection-test/Cdiff_toxin_genes/TcdAB-blastn-DB \
                                               -parse_seqids -dbtype nucl

#build bowtie2 index
/home/qinglong/softwares/bowtie2-binary-2.3.4.3/bowtie2-build --threads $THREADS TcdAB.fasta TcdAB_BT2-index

#build bowtie2 index for each gene; "make sure the header identifier is identical across strains of C. difficile for each gene"
cd $WorkDir/Cdiff_toxin_genes
for Gene in $WorkDir/Cdiff_toxin_genes/*.fasta
do
	GeneName=`basename $Gene .fasta`
	/home/qinglong/softwares/bowtie2-binary-2.3.4.3/bowtie2-build --threads $THREADS $Gene ''$GeneName'_BT2-index'
done

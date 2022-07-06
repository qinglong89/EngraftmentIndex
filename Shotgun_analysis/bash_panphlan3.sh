#!/bin/bash

WorkDir=/mnt/home1/@Untargeted_metagenomics/WGS_CDI/WGS_BGI_batch2/FMT-panphlan3

THREADS=80

cd $WorkDir

mkdir -p $WorkDir/panphlan3_results

for folder in $WorkDir/RawReadProcess/*
do
        cd $folder
        SampleName=`basename $folder`

        R1=`ls *_1.fq.gz`    #forward read
        R2=`ls *_2.fq.gz`    #reverse read

        #perform expected accumulated error-based filtering
        /home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --threads $THREADS \
                                                                  --fastq_filter $R1 \
                                                                  --fastqout ''$SampleName'_1_maxEE1.fq' \
                                                                  --fastq_truncee 1 \
                                                                  --fastq_stripleft 2 \
                                                                  --fastq_minlen 90 &

        /home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --threads $THREADS \
                                                                  --fastq_filter $R2 \
                                                                  --fastqout ''$SampleName'_2_maxEE1.fq' \
                                                                  --fastq_truncee 1 \
                                                                  --fastq_stripleft 2 \
                                                                  --fastq_minlen 90
        wait


	#run panphlan_map.py
	source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/biobakery3

#	for Species in Bacteroides_vulgatus Flavonifractor_plautii Blautia_wexlerae Dorea_longicatena Bacteroides_uniformis Bacteroides_ovatus Intestinibacter_bartlettii Ruminococcus_torques Anaerostipes_hadrus Blautia_obeum Faecalibacterium_prausnitzii
	for Species in Bacteroides_ovatus Blautia_obeum
	do
		mkdir -p $WorkDir/panphlan3_results/$Species

		panphlan_map.py -p /mnt/home1/qinglong/PanPhlan3_pangenomes/$Species/''$Species'_pangenome.tsv' \
				--indexes /mnt/home1/qinglong/PanPhlan3_pangenomes/$Species/$Species \
				-i ''$SampleName'_1_maxEE1.fq' \
				-o $WorkDir/panphlan3_results/$Species/''$SampleName'_forward.csv' \
				--nproc $THREADS &

		panphlan_map.py -p /mnt/home1/qinglong/PanPhlan3_pangenomes/$Species/''$Species'_pangenome.tsv' \
                                --indexes /mnt/home1/qinglong/PanPhlan3_pangenomes/$Species/$Species \
                                -i ''$SampleName'_2_maxEE1.fq' \
                                -o $WorkDir/panphlan3_results/$Species/''$SampleName'_reverse.csv' \
                                --nproc $THREADS
		wait
	done

	source /home/DataAnalysis/miniconda2/bin/deactivate /home/qinglong/Programs/biobakery3

	rm ''$SampleName'_1_maxEE1.fq' ''$SampleName'_2_maxEE1.fq' #save storage
done


#run panphlan_profiling.py
for map in $WorkDir/panphlan3_results/*
do
	mapid=`basename $map`

	#apply very sensitive mode for detection

	source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/biobakery3

	panphlan_profiling.py -i $map \
			      --min_coverage 1 --left_max 1.70 --right_min 0.30 \
			      -p /mnt/home1/qinglong/PanPhlan3_pangenomes/$mapid/''$mapid'_pangenome.tsv' \
			      --o_matrix $WorkDir/panphlan3_results/'panphlan3_'$mapid'.tsv' \
			      --o_covmat $WorkDir/panphlan3_results/'panphlan3_'$mapid'.coverage'

	source /home/DataAnalysis/miniconda2/bin/deactivate /home/qinglong/Programs/biobakery3
done



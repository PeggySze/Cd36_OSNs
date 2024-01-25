#!/bin/bash

## PBS configure
#PBS -N run_RNA_Velocity
#PBS -j oe
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system
START=$(date +%s.%N)

## Activate conda python3 environment
source /md01/shipy3/ori/miniconda3/bin/activate
conda activate python3

## Set variables
cellranger_dir=~/Cd36_OSNs/output/scRNA/cellranger/
output_dir=~/Cd36_OSNs/output/scRNA/Velocyto/
samples=(19d_rep1 19d_rep2 8w_rep1 8w_rep2)

for i in `seq 0 3`;
do
	if [ ! -d ${output_dir}${samples[$i]} ]; then
		mkdir ${output_dir}${samples[$i]}
	fi

	cp ${cellranger_dir}${samples[$i]}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ${output_dir}${samples[$i]}/
	gunzip -f ${output_dir}${samples[$i]}/barcodes.tsv.gz

	cp ${cellranger_dir}${samples[$i]}/outs/possorted_genome_bam.bam ${output_dir}${samples[$i]}/

	# Generate spliced and unspliced matrices by Velocyto
	velocyto run \
	-b ${output_dir}${samples[$i]}/barcodes.tsv \
	-m /md01/shipy3/ori/DB/mm10/annotation/mm10_rmsk.gtf \
	-o ${output_dir}${samples[$i]} \
	${output_dir}${samples[$i]}/possorted_genome_bam.bam \
	/md01/shipy3/ori/DB/10X_reference/cellranger-arc-mm10/genes.gtf
done

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration
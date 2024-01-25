#!/bin/bash

## PBS configure
# PBS -N mouse_OE_scATAC_cellranger
# PBS -j oe
# PSB -q batch
# PBS -S /bin/sh
# PBS -l nodes=1:ppn=20

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system

## Set variables
input_dir=~/Cd36_OSNs/input/data/scATAC/
output_dir=~/Cd36_OSNs/output/scATAC/cellranger/
reference_dir=~/DB/10X_reference/
samples=(19d_0818 8w_0720)

for sample in ${samples[@]};
do
	# create the output directory
	mkdir -p ${output_dir}${sample}
	cd ${output_dir}${sample}

	## generate single cell feature counts by cellranger-atac count
	/public/home/yangjw28/software/cellranger-atac-2.0.0/cellranger-atac count --id=${sample} \
	--fastqs=${input_dir}${sample} \
	--reference=/md01/yangjw28/ori/genomeDB/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
	--sample=${sample} \
	--localcores=20
done

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration


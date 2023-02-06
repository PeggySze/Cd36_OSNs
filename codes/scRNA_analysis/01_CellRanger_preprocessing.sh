#!/bin/bash

## PBS configure
# PBS -N mouse_OE_scRNA_cellranger
# PBS -j oe
# PSB -q batch
# PBS -S /bin/sh
# PBS -l nodes=1:ppn=20

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system

## Source module envrionment and load tools
source /etc/profile.d/modules.sh
module load cellranger/4.0.0

## Set variables
input_dir=~/Cd36_OSNs/input/data/scRNA/
output_dir=~/Cd36_OSNs/output/scRNA/cellranger/
reference_dir=~/DB/10X_reference/
samples=(19d_rep1 19d_rep2 8w_rep1 8w_rep2)

## generate single cell feature counts by cellranger count
cd ${output_dir}
for sample in ${samples[@]};
do
	cellranger count --id=${sample} \
	--fastqs=${input_dir}${sample}/ \
	--sample=${sample} \
	--transcriptome=${reference_dir}refdata-gex-mm10-2020-A/ 
done

## Unload tools
module unload cellranger/4.0.0

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration


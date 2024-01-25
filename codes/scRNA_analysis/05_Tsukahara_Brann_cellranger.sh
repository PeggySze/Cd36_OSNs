#!/bin/bash

## PBS configure
#PBS -N Tsukahara_Brann_cellranger
#PBS -j oe
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system
START=$(date +%s.%N)

## Source module envrionment and load tools
export MODULEPATH=$MODULEPATH:/public/home/shipy3/modulefiles
source /etc/profile.d/modules.sh
module load cellranger/4.0.0

## Set variables
input_dir=~/Cd36_OSNs/input/Tsukahara_Brann_OSN/data/raw/fastq/
output_dir=~/Cd36_OSNs/input/Tsukahara_Brann_OSN/data/cellranger/
reference_dir=/md01/yangjw28/db/mm10-scRNA/refdata-gex-mm10-2020-A/

## generate single cell feature counts by cellranger count
cd ${output_dir}
cellranger count --id=homecage_1 \
--fastqs=${input_dir}homecage_1/baseline-1_0_1_HM3C5DSXY/,${input_dir}homecage_1/baseline-1_0_1_HM533DSXY/ \
--sample=bamtofastq \
--transcriptome=${reference_dir} \
--localcores=10 \
--localmem=64

cellranger count --id=homecage_2 \
--fastqs=${input_dir}homecage_2/baseline-2_0_1_HM3C5DSXY/,${input_dir}homecage_2/baseline-2_0_1_HM533DSXY/ \
--sample=bamtofastq \
--transcriptome=${reference_dir} \
--localcores=10 \
--localmem=4

cellranger count --id=homecage_3 \
--fastqs=${input_dir}homecage_3/baseline-3_0_1_HM3C5DSXY/,${input_dir}homecage_3/baseline-3_0_1_HM533DSXY/ \
--sample=bamtofastq \
--transcriptome=${reference_dir} \
--localcores=10 \
--localmem=64

cellranger count --id=homecage_4 \
--fastqs=${input_dir}homecage_4/baseline-4_0_1_HM3C5DSXY/,${input_dir}homecage_4/baseline-4_0_1_HM533DSXY/ \
--sample=bamtofastq \
--transcriptome=${reference_dir} \
--localcores=10 \
--localmem=64

cellranger count --id=homecage_5 \
--fastqs=${input_dir}homecage_5/baseline-9_0_1_H22KHDRXY/,${input_dir}homecage_5/baseline-9_0_1_HVWHYDRXX/ \
--sample=bamtofastq \
--transcriptome=${reference_dir} \
--localcores=10 \
--localmem=64

cellranger count --id=homecage_6 \
--fastqs=${input_dir}homecage_6/baseline-10_0_1_H22KHDRXY/,${input_dir}homecage_6/baseline-10_0_1_HVWHYDRXX/ \
--sample=bamtofastq \
--transcriptome=${reference_dir} \
--localcores=10 \
--localmem=64

## Unload tools
module unload cellranger/4.0.0

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration
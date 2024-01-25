#!/bin/bash

## PBS configure
#PBS -N Tsukahara_Brann_bamtofastq
#PBS -j oe
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20

## Source module envrionment and load tools
export MODULEPATH=$MODULEPATH:/public/home/shipy3/modulefiles
source /etc/profile.d/modules.sh
module load bamtofastq/1.2.0

cd ~/Cd36_OSNs/input/Tsukahara_Brann_OSN/data/raw/
for i in `seq 1 6`;
do
        bamtofastq-1.2.0 --nthreads 10 homecage_${i}_possorted_genome_bam.bam ~/Cd36_OSNs/input/Tsukahara_Brann_OSN/data/raw/fastq/homecage_${i} ;
done
#!/bin/bash

#PBS -N binning
#PBS -P berger_common
#PBS -j oe
#PBS -q workq
#PBS -J 1-12
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -o /home/GMI/sean.montgomery/scripts/out-logs/binning#^array_index^.txt

# === begin ENVIRONMENT SETUP ===

#1. set directory with lists of bam and sample names
sample_list=/home/GMI/sean.montgomery/inputfiles/cutrun_merged.txt
#4. set directory containing sample folders
work_folder=$WORK/cutrun/Tak1v4/merged
#5. bin size
binsize=100000

#4. sample name
NAME=`sed -n "${PBS_ARRAY_INDEX} p" $sample_list | awk '{print $1}'`


# Load the required modules
module load deepTools/2.2.4-foss-2015a-Python-2.7.9
# ... and then change to  working directory:
cd $work_folder
# === end ENVIRONMENT SETUP ===

# bamCoverage -bs $binsize -of bedgraph -b $NAME/$NAME.sizednuc150.bam --normalizeTo1x 218000000  -o $NAME/$NAME.sizednuc150.100kb.bed 
# binsize=10000
# bamCoverage -bs $binsize -of bedgraph -b $NAME/$NAME.sizednuc150.bam --normalizeTo1x 218000000  -o $NAME/$NAME.sizednuc150.10kb.bed
binsize=50000
bamCoverage -bs $binsize -of bedgraph -b $NAME/$NAME.sizednuc150.bam --normalizeTo1x 218000000  -o $NAME/$NAME.sizednuc150.50kb.bed

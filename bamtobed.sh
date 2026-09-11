#!/bin/bash

#PBS -N ChIPseq_BED
#PBS -P berger_common
#PBS -j oe
#PBS -q workq
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -o /lustre/scratch/users/sean.montgomery/chipseqTak1v4.1/logs/bed.txt

#Script to trim adapters then align ChIP-seq paired end data with Bowtie2

# === begin ENVIRONMENT SETUP ===

#1. set directory with lists of bam and sample names
list_dir=/home/GMI/$USER/Scripts/input_files
#2. list of bam file names on separate lines (eg copied from excel)
bam_list=chipseq_bamfiles.txt
#3. list of sample names on separate lines (eg copied from excel)
sample_list=chipseq_samplenames.txt
#4. set directory containing sample folders
work_folder=$WORK/chipseq/Tak1v4.1

bamfiles=$list_dir/$bam_list #list of bam files on separate lines (eg copied from excel)
F=`sed -n "3 p" $bamfiles`

samplenames=$list_dir/$sample_list #list of sample names on separate lines (eg 23955)
NAME=`sed -n "3 p" $samplenames`



# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
module load BEDTools/v2.17.0-goolf-1.4.10
module load cutadapt/1.9-goolf-1.4.10-Python-2.7.5
module load Bowtie2/2.1.0-goolf-1.4.10
module load deepTools/2.2.4-foss-2015a-Python-2.7.9
module load Picard/1.141
module load R/3.3.1-foss-2016b
# ... and then change to  working directory:
cd $work_folder

export TMPDIR=$work_folder/tmp
mkdir -p $TMPDIR
samtools sort -n ${NAME}/${NAME}.sorted_uniq.bam | bamToBed -i stdin | sort -k1,1 -k2,2n -T $TMPDIR > ${NAME}/${NAME}.sorted_uniq.bed

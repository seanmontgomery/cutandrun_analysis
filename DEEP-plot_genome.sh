#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --output=/users/sean.montgomery/logs/deeptoolsgenomeplot.txt
# === end SBATCH directives ===


# === begin ENVIRONMENT SETUP ===

#4. set directory containing sample folders
work_folder=/scratch-cbe/users/sean.montgomery/cutrun/homeogypsy

# Load the required modules
module load deeptools/2.5.4-foss-2018b-python-2.7.15
# module load deeptools/3.1.2-foss-2018b-python-2.7.15

# ... and then change to  working directory:
cd $work_folder
# mkdir -p $work_folder/profiles

  computeMatrix reference-point \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  -R $WORK/lab/Tak1v4/annotations/Tak1v4.gene.bed \
  -S \
	H3K9me1_Tak1v4/H3K9me1_Tak1v4.sizednuc150.bw \
	H3K27me1_Tak1v4/H3K27me1_Tak1v4.sizednuc150.bw \
	H3K9ac_Tak1v4/H3K9ac_Tak1v4.sizednuc150.bw \
	H3K14ac_Tak1v4/H3K14ac_Tak1v4.sizednuc150.bw \
	H3K4me1_Tak1v4/H3K4me1_Tak1v4.sizednuc150.bw \
	H3K36me3_Tak1v4/H3K36me3_Tak1v4.sizednuc150.bw \
	H3K4me3_Tak1v4/H3K4me3_Tak1v4.sizednuc150.bw \
  	H2A.Z_Tak1v4/H2A.Z_Tak1v4.sizednuc150.bw \
	H3K27me3_Tak1v4/H3K27me3_Tak1v4.sizednuc150.bw \
	../H3-1_Tak1v4/H3-1_Tak1v4.sizednuc150.bw \
	--skipZeros \
	--blackListFileName profiles/blacklist.bed \
  -o profiles/selected.reference-point.TSS.gz \
  --outFileSortedRegions profiles/selected.reference-point.TSS.regions.bed
  plotHeatmap -m profiles/selected.reference-point.TSS.gz \
          -out profiles/selected.TSS.heatmap.k5.pdf \
          --colorMap RdBu \
          --whatToShow 'heatmap and colorbar' \
          --zMin 0 --zMax 3 \
          --outFileSortedRegions profiles/selected.TSS.profiledata.k5.heatmap.txt \
          --kmeans 5
 computeMatrix reference-point \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  -R profiles/selected.TSS.profiledata.k5.heatmap.txt \
  -S \
	H3K9me1_Tak1v4/H3K9me1_Tak1v4.sizednuc150.bw \
	H3K27me1_Tak1v4/H3K27me1_Tak1v4.sizednuc150.bw \
	H3K9ac_Tak1v4/H3K9ac_Tak1v4.sizednuc150.bw \
	H3K14ac_Tak1v4/H3K14ac_Tak1v4.sizednuc150.bw \
	H3K4me1_Tak1v4/H3K4me1_Tak1v4.sizednuc150.bw \
	H3K36me3_Tak1v4/H3K36me3_Tak1v4.sizednuc150.bw \
	H3K4me3_Tak1v4/H3K4me3_Tak1v4.sizednuc150.bw \
  	H2A.Z_Tak1v4/H2A.Z_Tak1v4.sizednuc150.bw \
	H3K27me3_Tak1v4/H3K27me3_Tak1v4.sizednuc150.bw \
	../H3-1_Tak1v4/H3-1_Tak1v4.sizednuc150.bw \
	--skipZeros \
  -o profiles/selected.reference-point.TSS.clusters.gz 

 plotProfile -m profiles/selected.reference-point.TSS.clusters.gz \
          --perGroup \
          --yMin 0 --yMax 9 --plotHeight 18 --plotWidth 13 \
          --refPointLabel TSS \
          --colors "#CB181D" "#FB6A4A" "#005A32" "#238B45" "#74C476" "#A1D99B" "#084594" "#4292C6" "#C6DBEF" "#cdcd00" \
          --samplesLabel H3K9me1 H3K27me1 H3K9ac H3K14ac H3K4me1 H3K36me3 H3K4me3 H2A.Z H3K27me3 H3 \
          --numPlotsPerRow 1 \
          -out profiles/selected.TSS.plotProfile.k5.clusters.pdf
 plotHeatmap -m profiles/selected.reference-point.TSS.clusters.gz \
          -out profiles/selected.TSS.heatmap.k5.clusters.pdf \
          --colorMap coolwarm \
          --whatToShow 'heatmap and colorbar' \
          --sortRegions no \
          --zMin 0 --zMax 4 



computeMatrix reference-point \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  -R $WORK/lab/Tak1v4/annotations/Tak1v4.repeat.bed  \
  -S \
	H3K9me1_Tak1v4/H3K9me1_Tak1v4.sizednuc150.bw \
	H3K27me1_Tak1v4/H3K27me1_Tak1v4.sizednuc150.bw \
	H3K9ac_Tak1v4/H3K9ac_Tak1v4.sizednuc150.bw \
	H3K14ac_Tak1v4/H3K14ac_Tak1v4.sizednuc150.bw \
	H3K4me1_Tak1v4/H3K4me1_Tak1v4.sizednuc150.bw \
	H3K36me3_Tak1v4/H3K36me3_Tak1v4.sizednuc150.bw \
	H3K4me3_Tak1v4/H3K4me3_Tak1v4.sizednuc150.bw \
  	H2A.Z_Tak1v4/H2A.Z_Tak1v4.sizednuc150.bw \
	H3K27me3_Tak1v4/H3K27me3_Tak1v4.sizednuc150.bw \
	../H3-1_Tak1v4/H3-1_Tak1v4.sizednuc150.bw \
	--skipZeros \
	--blackListFileName profiles/blacklist.bed \
  -o profiles/selected.reference-point.repeatTSS.gz \
  --outFileSortedRegions profiles/selected.reference-point.repeatTSS.regions.bed
  plotHeatmap -m profiles/selected.reference-point.repeatTSS.gz \
          -out profiles/selected.repeatTSS.heatmap.k5.pdf \
          --colorMap RdBu \
          --whatToShow 'heatmap and colorbar' \
          --zMin 0 --zMax 3 \
          --outFileSortedRegions profiles/selected.repeat.TSS.profiledata.k5.heatmap.txt \
          --kmeans 5
 computeMatrix reference-point \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  -R profiles/selected.repeat.TSS.profiledata.k5.heatmap.txt \
  -S \
	H3K9me1_Tak1v4/H3K9me1_Tak1v4.sizednuc150.bw \
	H3K27me1_Tak1v4/H3K27me1_Tak1v4.sizednuc150.bw \
	H3K9ac_Tak1v4/H3K9ac_Tak1v4.sizednuc150.bw \
	H3K14ac_Tak1v4/H3K14ac_Tak1v4.sizednuc150.bw \
	H3K4me1_Tak1v4/H3K4me1_Tak1v4.sizednuc150.bw \
	H3K36me3_Tak1v4/H3K36me3_Tak1v4.sizednuc150.bw \
	H3K4me3_Tak1v4/H3K4me3_Tak1v4.sizednuc150.bw \
  	H2A.Z_Tak1v4/H2A.Z_Tak1v4.sizednuc150.bw \
	H3K27me3_Tak1v4/H3K27me3_Tak1v4.sizednuc150.bw \
	../H3-1_Tak1v4/H3-1_Tak1v4.sizednuc150.bw \
	--skipZeros \
  -o profiles/selected.reference-point.repeatTSS.clusters.gz 

 plotProfile -m profiles/selected.reference-point.repeatTSS.clusters.gz  \
          --perGroup \
          --yMin 0 --yMax 9 --plotHeight 18 --plotWidth 13 \
          --refPointLabel TSS \
          --colors "#CB181D" "#FB6A4A" "#005A32" "#238B45" "#74C476" "#A1D99B" "#084594" "#4292C6" "#C6DBEF" "#cdcd00" \
          --samplesLabel H3K9me1 H3K27me1 H3K9ac H3K14ac H3K4me1 H3K36me3 H3K4me3 H2A.Z H3K27me3 H3 \
          --numPlotsPerRow 1 \
          -out profiles/selected.repeatTSS.plotProfile.k5.clusters.pdf
 plotHeatmap -m profiles/selected.reference-point.repeatTSS.clusters.gz  \
          -out profiles/selected.repeatTSS.heatmap.k5.clusters.pdf \
          --colorMap coolwarm \
          --whatToShow 'heatmap and colorbar' \
          --sortRegions no \
          --zMin 0 --zMax 4

 computeMatrix reference-point \
  --referencePoint TES \
  -b 1000 -a 1000 \
  -R profiles/selected.TSS.profiledata.k5.heatmap.txt \
  -S \
	H3K9me1_Tak1v4/H3K9me1_Tak1v4.sizednuc150.bw \
	H3K27me1_Tak1v4/H3K27me1_Tak1v4.sizednuc150.bw \
	H3K9ac_Tak1v4/H3K9ac_Tak1v4.sizednuc150.bw \
	H3K14ac_Tak1v4/H3K14ac_Tak1v4.sizednuc150.bw \
	H3K4me1_Tak1v4/H3K4me1_Tak1v4.sizednuc150.bw \
	H3K36me3_Tak1v4/H3K36me3_Tak1v4.sizednuc150.bw \
	H3K4me3_Tak1v4/H3K4me3_Tak1v4.sizednuc150.bw \
  	H2A.Z_Tak1v4/H2A.Z_Tak1v4.sizednuc150.bw \
	H3K27me3_Tak1v4/H3K27me3_Tak1v4.sizednuc150.bw \
	../H3-1_Tak1v4/H3-1_Tak1v4.sizednuc150.bw \
	--skipZeros \
  -o profiles/selected.reference-point.TES.clusters.gz 

 plotProfile -m profiles/selected.reference-point.TES.clusters.gz \
          --perGroup \
          --yMin 0 --yMax 9 --plotHeight 18 --plotWidth 13 \
          --refPointLabel TES \
          --colors "#CB181D" "#FB6A4A" "#005A32" "#238B45" "#74C476" "#A1D99B" "#084594" "#4292C6" "#C6DBEF" "#cdcd00" \
          --samplesLabel H3K9me1 H3K27me1 H3K9ac H3K14ac H3K4me1 H3K36me3 H3K4me3 H2A.Z H3K27me3 H3 \
          --numPlotsPerRow 1 \
          -out profiles/selected.TES.plotProfile.k5.clusters.pdf
 plotHeatmap -m profiles/selected.reference-point.TES.clusters.gz \
          -out profiles/selected.TES.heatmap.k5.clusters.pdf \
          --colorMap coolwarm \
          --whatToShow 'heatmap and colorbar' \
          --sortRegions no \
          --zMin 0 --zMax 4


 computeMatrix reference-point \
  --referencePoint TES \
  -b 1000 -a 1000 \
  -R profiles/selected.repeat.TSS.profiledata.k5.heatmap.txt \
  -S \
	H3K9me1_Tak1v4/H3K9me1_Tak1v4.sizednuc150.bw \
	H3K27me1_Tak1v4/H3K27me1_Tak1v4.sizednuc150.bw \
	H3K9ac_Tak1v4/H3K9ac_Tak1v4.sizednuc150.bw \
	H3K14ac_Tak1v4/H3K14ac_Tak1v4.sizednuc150.bw \
	H3K4me1_Tak1v4/H3K4me1_Tak1v4.sizednuc150.bw \
	H3K36me3_Tak1v4/H3K36me3_Tak1v4.sizednuc150.bw \
	H3K4me3_Tak1v4/H3K4me3_Tak1v4.sizednuc150.bw \
  	H2A.Z_Tak1v4/H2A.Z_Tak1v4.sizednuc150.bw \
	H3K27me3_Tak1v4/H3K27me3_Tak1v4.sizednuc150.bw \
	../H3-1_Tak1v4/H3-1_Tak1v4.sizednuc150.bw \
	--skipZeros \
  -o profiles/selected.reference-point.repeatTES.clusters.gz 

 plotProfile -m profiles/selected.reference-point.repeatTES.clusters.gz  \
          --perGroup \
          --yMin 0 --yMax 9 --plotHeight 18 --plotWidth 13 \
          --refPointLabel TES \
          --colors "#CB181D" "#FB6A4A" "#005A32" "#238B45" "#74C476" "#A1D99B" "#084594" "#4292C6" "#C6DBEF" "#cdcd00" \
          --samplesLabel H3K9me1 H3K27me1 H3K9ac H3K14ac H3K4me1 H3K36me3 H3K4me3 H2A.Z H3K27me3 H3 \
          --numPlotsPerRow 1 \
          -out profiles/selected.repeatTES.plotProfile.k5.clusters.pdf
 plotHeatmap -m profiles/selected.reference-point.repeatTES.clusters.gz  \
          -out profiles/selected.repeatTES.heatmap.k5.clusters.pdf \
          --colorMap coolwarm \
          --whatToShow 'heatmap and colorbar' \
          --sortRegions no \
          --zMin 0 --zMax 4



##Plot with TPM ordered genes
 computeMatrix reference-point \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  -R profiles/selected.TSS.profiledata.k5.heatmap.tpmorder.txt \
  -S \
	H3K9me1_Tak1v4/H3K9me1_Tak1v4.sizednuc150.bw \
	H3K27me1_Tak1v4/H3K27me1_Tak1v4.sizednuc150.bw \
	H3K9ac_Tak1v4/H3K9ac_Tak1v4.sizednuc150.bw \
	H3K14ac_Tak1v4/H3K14ac_Tak1v4.sizednuc150.bw \
	H3K4me1_Tak1v4/H3K4me1_Tak1v4.sizednuc150.bw \
	H3K36me3_Tak1v4/H3K36me3_Tak1v4.sizednuc150.bw \
	H3K4me3_Tak1v4/H3K4me3_Tak1v4.sizednuc150.bw \
  	H2A.Z_Tak1v4/H2A.Z_Tak1v4.sizednuc150.bw \
	H3K27me3_Tak1v4/H3K27me3_Tak1v4.sizednuc150.bw \
	H3-1_Tak1v4/H3-1_Tak1v4.sizednuc150.bw \
	--skipZeros \
  -o profiles/selected.reference-point.TSS.clusters.tpmorder.gz 

 plotProfile -m profiles/selected.reference-point.TSS.clusters.tpmorder.gz \
          --perGroup \
          --yMin 0 --yMax 9 --plotHeight 18 --plotWidth 13 \
          --refPointLabel TSS \
          --colors "#CB181D" "#FB6A4A" "#005A32" "#238B45" "#74C476" "#A1D99B" "#084594" "#4292C6" "#C6DBEF" "#cdcd00" \
          --samplesLabel H3K9me1 H3K27me1 H3K9ac H3K14ac H3K4me1 H3K36me3 H3K4me3 H2A.Z H3K27me3 H3 \
          --numPlotsPerRow 1 \
          -out profiles/selected.TSS.plotProfile.k5.clusters.tpmorder.pdf
 plotHeatmap -m profiles/selected.reference-point.TSS.clusters.tpmorder.gz \
          -out profiles/selected.TSS.heatmap.k5.clusters.tpmorder.pdf \
          --colorMap coolwarm \
          --whatToShow 'heatmap and colorbar' \
          --sortRegions no \
          --zMin 0 --zMax 4 




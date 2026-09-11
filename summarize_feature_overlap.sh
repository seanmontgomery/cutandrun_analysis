#!/bin/bash

#Summarize ratio of overlapping bp of peak and features
#Created 8.1.2018 by Sean Montgomery

##Sample usage:
## for f in *bedgraph; do /volumes/berger/user/sean.montgomery/Documents/Scripts/chipseq/summarize_feature_overlap.sh $f; done

# #1. Sample bedgraph file from homer
# file1=$1
# #2 gff3 file of genomic features
# file2=/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.merged.gff 

# # Summarize feature overlap
# # feature=$2 #CDS, five_prime_UTR, three_prime_UTR, exon, mRNA, intron

# # coverage_length=`bedtools intersect -a $file1 -b $file2 -wo | grep $feature | awk '{sum+=$14} END {print sum}'`
# # feature_length=`grep $feature $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# ip_length=`awk '{a=$3-$2;print $0,a;}' $file1 | awk '{sum+=$5} END {print sum}'`


# cov_mRNA=`bedtools intersect -a $file1 -b $file2 -wo | grep mRNA | awk '{sum+=$14} END {print sum}'`
# cov_three_prime=`bedtools intersect -a $file1 -b $file2 -wo | grep three_prime_UTR | grep '\.1' | awk '{sum+=$14} END {print sum}'`
# cov_five_prime=`bedtools intersect -a $file1 -b $file2 -wo | grep five_prime_UTR | grep '\.1' | awk '{sum+=$14} END {print sum}'`
# cov_cds=`bedtools intersect -a $file1 -b $file2 -wo | grep CDS | grep '\.1' | awk '{sum+=$14} END {print sum}'`
# cov_repeat=`bedtools intersect -a $file1 -b $file2 -wo | grep repeat_region | awk '{sum+=$14} END {print sum}'`
# cov_gene=`bedtools intersect -a $file1 -b $file2 -wo | grep gene | awk '{sum+=$14} END {print sum}'`
# cov_rRNA=`bedtools intersect -a $file1 -b $file2 -wo | grep rRNA | awk '{sum+=$14} END {print sum}'`
# cov_simple_repeat=`bedtools intersect -a $file1 -b $file2 -wo | grep Simple | awk '{sum+=$14} END {print sum}'`
# cov_satellite=`bedtools intersect -a $file1 -b $file2 -wo | grep Satellite | awk '{sum+=$14} END {print sum}'`
# cov_unknown=`bedtools intersect -a $file1 -b $file2 -wo | grep Unknown | awk '{sum+=$14} END {print sum}'`
# cov_other_repeat=`bedtools intersect -a $file1 -b $file2 -wo | grep repeat_region | grep -v Simple | grep -v Satellite | grep -v Unknown | awk '{sum+=$14} END {print sum}'`

# mRNA=`grep mRNA $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# three_prime=`grep three_prime_UTR $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# five_prime=`grep five_prime_UTR $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# cds=`grep CDS $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# repeat=`grep repeat_region $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# gene=`grep gene $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# rRNA=`grep rRNA $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# simple_repeat=`grep Simple $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# satellite=`grep Satellite $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# unknown=`grep Unknown $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
# other_repeat=`grep repeat_region $file2 | grep -v Simple | grep -v Satellite | grep -v Unknown | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`

# intron=`echo "scale=3; $mRNA - $three_prime - $five_prime - $cds" | bc -l`
# cov_intron=`echo "scale=3; $cov_mRNA - $cov_three_prime - $cov_five_prime - $cov_cds" | bc -l`
# gene_other=`echo "scale=3; $gene - $three_prime - $five_prime - $cds" | bc -l`
# cov_gene_other=`echo "scale=3; $cov_gene - $cov_three_prime - $cov_five_prime - $cov_cds" | bc -l`

# echo "Total length of" $file1 ":" $ip_length

# echo "Total coverage length for 5_prime_UTR in" $file1 ":" $cov_five_prime
# fiver=`echo "scale=3; $cov_five_prime / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for 5_prime_UTR:" $five_prime
# fiver2=`echo "scale=3; $cov_five_prime / $five_prime" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for CDS in" $file1 ":" $cov_cds
# threer=`echo "scale=3; $cov_cds / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for CDS:" $cds
# threer2=`echo "scale=3; $cov_cds / $cds" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for 3_prime_UTR in" $file1 ":" $cov_three_prime
# cdser=`echo "scale=3; $cov_three_prime / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for 3_prime_UTR:" $three_prime
# cdser2=`echo "scale=3; $cov_three_prime / $three_prime" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for intron in" $file1 ":" $cov_intron
# introner=`echo "scale=3; $cov_intron / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for intron:" $intron
# introner2=`echo "scale=3; $cov_intron / $intron" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for repeats in" $file1 ":" $cov_repeat
# repeater=`echo "scale=3; $cov_repeat / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for repeats:" $repeat
# repeater2=`echo "scale=3; $cov_repeat / $repeat" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for genes in" $file1 ":" $cov_gene
# gener=`echo "scale=3; $cov_gene / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for genes:" $gene
# gener2=`echo "scale=3; $cov_gene / $gene" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for rRNA in" $file1 ":" $cov_rRNA
# rRNAer=`echo "scale=3; $cov_rRNA / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for rRNA:" $rRNA
# rRNAer2=`echo "scale=3; $cov_rRNA / $rRNA" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for gene in" $file1 ":" $cov_gene_other
# gene_otherer=`echo "scale=3; $cov_gene_other / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for gene:" $gene_other
# gene_otherer2=`echo "scale=3; $cov_gene_other / $gene_other" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for repeats in" $file1 ":" $cov_simple_repeat
# simpler=`echo "scale=3; $cov_simple_repeat / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for repeats:" $simple_repeat
# simpler2=`echo "scale=3; $cov_simple_repeat / $simple_repeat" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for repeats in" $file1 ":" $cov_satellite
# satelliter=`echo "scale=3; $cov_satellite / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for repeats:" $satellite
# satelliter2=`echo "scale=3; $cov_satellite / $satellite" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for repeats in" $file1 ":" $cov_unknown
# unknowner=`echo "scale=3; $cov_unknown / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for repeats:" $unknown
# unknowner2=`echo "scale=3; $cov_unknown / $unknown" | bc -l` #Portion of feature that is sample
# echo ""

# echo "Total coverage length for repeats in" $file1 ":" $cov_other_repeat
# repeat_otherer=`echo "scale=3; $cov_other_repeat / $ip_length" | bc -l` #Portion of sample that is feature
# echo "Total feature length for repeats:" $other_repeat
# repeat_otherer2=`echo "scale=3; $cov_other_repeat / $other_repeat" | bc -l` #Portion of feature that is sample
# echo ""

# NAME=`basename $file1 _Tak1v4.regions.narrow.bedgraph`
# # echo -e "Sample \t 5primeUTR \t CDS \t 3primeUTR \t intron \t repeats \t gene \t rRNA \t gene \t simple \t satellite \t unknown \t repeat" > sample_overlap_summary.txt
# echo -e "$NAME \t $fiver \t $cdser \t $threer \t $introner \t $repeater \t $gener \t $rRNAer \t $gene_otherer \t $simpler \t $satelliter \t $unknowner \t $repeat_otherer"  >> sample_overlap_summary.txt
# # echo -e "Sample \t 5primeUTR \t CDS \t 3primeUTR \t intron \t repeats \t gene \t rRNA \t gene \t simple \t satellite \t unknown \t repeat" > feature_overlap_summary.txt
# echo -e "$NAME \t $fiver2 \t $cdser2 \t $threer2 \t $introner2 \t $repeater2 \t $gener2 \t $rRNAer2 \t $gene_otherer2 \t $simpler2 \t $satelliter2 \t $unknowner2 \t $repeat_otherer2" >> feature_overlap_summary.txt




# #By clusters
# #1. Determine genes in clusters
# # grep 'cluster 1' selected.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.TSS.profiledata.k5.heatmap.cluster1.txt
# # grep 'cluster 2' selected.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.TSS.profiledata.k5.heatmap.cluster2.txt
# # grep 'cluster 3' selected.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.TSS.profiledata.k5.heatmap.cluster3.txt
# # grep 'cluster 4' selected.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.TSS.profiledata.k5.heatmap.cluster4.txt
# # grep 'cluster 5' selected.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.TSS.profiledata.k5.heatmap.cluster5.txt
# #1. Sample bedgraph file from homer
# file1=$1
# # Summarize feature overlap
# ip_length=`awk '{a=$3-$2;print $0,a;}' $file1 | awk '{sum+=$5} END {print sum}'`
# for i in {1..5}; do
# 	file2=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.cluster$i.txt
# 	cov_gene=`bedtools intersect -a $file1 -b $file2 -wo | awk '{sum+=$13} END {print sum}'`
# 	gene=`awk '{a=$3-$2;print $0,a;}' $file2 | awk '{sum+=$9} END {print sum}'`
# 	echo "Total coverage length for genes in" $file1 ":" $cov_gene
# 	echo "Total feature length for genes:" $gene
# 	echo ""
# 	eval "samplecluster$i=`echo "scale=3; $cov_gene / $ip_length" | bc -l`"
# 	eval "featurecluster$i=`echo "scale=3; $cov_gene / $gene" | bc -l`"
# done
# echo "Total length of" $file1 ":" $ip_length
# NAME=`basename $file1 _Tak1v4.regions.narrow.bedgraph`
# # echo -e "Sample \t Cluster_1 \t Cluster_2 \t Cluster_3 \t Cluster_4 \t Cluster_5" > sample_overlap_summary.clusters.txt
# echo -e "$NAME \t $samplecluster1 \t $samplecluster2 \t $samplecluster3 \t $samplecluster4 \t $samplecluster5" >> sample_overlap_summary.clusters.txt
# # echo -e "Sample \t Cluster_1 \t Cluster_2 \t Cluster_3 \t Cluster_4 \t Cluster_5" > feature_overlap_summary.clusters.txt
# echo -e "$NAME \t $featurecluster1 \t $featurecluster2 \t $featurecluster3 \t $featurecluster4 \t $featurecluster5" >> feature_overlap_summary.clusters.txt





# ##By quintile
# ##1. Determine quintiles
# # quintiles <- quantile(gametop$mean, probs = seq(0,1,0.2))
# # bedlist <- read.table("/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.gene.bed")
# # bedlist$TPM <- gametop$mean
# # write.table(subset(bedlist,bedlist$TPM<quintiles[[2]]), file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile1.txt",sep='\t',col.names = F,row.names = F)
# # write.table(subset(bedlist,bedlist$TPM>=quintiles[[5]]), file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile5.txt",sep='\t',col.names = F,row.names = F)
# # write.table(subset(bedlist,bedlist$TPM>=quintiles[[4]] & bedlist$TPM<quintiles[[5]]), file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile4.txt",sep='\t',col.names = F,row.names = F)
# # write.table(subset(bedlist,bedlist$TPM>=quintiles[[3]] & bedlist$TPM<quintiles[[4]]), file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile3.txt",sep='\t',col.names = F,row.names = F)
# # write.table(subset(bedlist,bedlist$TPM>=quintiles[[2]] & bedlist$TPM<quintiles[[3]]), file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile2.txt",sep='\t',col.names = F,row.names = F)
# ##2. Sample bedgraph file from homer
# file1=$1
# # Summarize feature overlap
# ip_length=`awk '{a=$3-$2;print $0,a;}' $file1 | awk '{sum+=$5} END {print sum}'`
# for i in {1..5}; do
# 	file2=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile$i.txt
# 	sed -i -e 's/"//g' $file2
# 	cov_gene=`bedtools intersect -a $file1 -b $file2 -wo | awk '{sum+=$11} END {print sum}'`
# 	gene=`awk '{a=$3-$2;print $0,a;}' $file2 | awk '{sum+=$7} END {print sum}'`
# 	echo "Total coverage length for genes in" $file1 ":" $cov_gene
# 	echo "Total feature length for genes:" $gene
# 	echo ""
# 	eval "samplecluster$i=`echo "scale=3; $cov_gene / $ip_length" | bc -l`"
# 	eval "featurecluster$i=`echo "scale=3; $cov_gene / $gene" | bc -l`"
# done
# echo "Total length of" $file1 ":" $ip_length
# NAME=`basename $file1 _Tak1v4.regions.narrow.bedgraph`
# # echo -e "Sample \t Quintile_1 \t Quintile_2 \t Quintile_3 \t Quintile_4 \t Quintile_5" > sample_overlap_summary.quintiles.txt
# echo -e "$NAME \t $samplecluster1 \t $samplecluster2 \t $samplecluster3 \t $samplecluster4 \t $samplecluster5" >> sample_overlap_summary.quintiles.txt
# # echo -e "Sample \t Quintile_1 \t Quintile_2 \t Quintile_3 \t Quintile_4 \t Quintile_5" > feature_overlap_summary.quintiles.txt
# echo -e "$NAME \t $featurecluster1 \t $featurecluster2 \t $featurecluster3 \t $featurecluster4 \t $featurecluster5" >> feature_overlap_summary.quintiles.txt


# ##4. Summarize repeat cluster overlap
# echo -e "Cluster \t Quintile1 \t Quintile2 \t Quintile3 \t Quintile4 \t Quintile5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/sample_overlap_summary.atac.quintiles.cluster.repeat.txt
# echo -e "Cluster \t Quintile1 \t Quintile2 \t Quintile3 \t Quintile4 \t Quintile5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/feature_overlap_summary.atac.quintiles.cluster.repeat.txt
# for i in {1..5}; do
# 	file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.cluster$i.bed.txt
# 	ip_length=`awk '{a=$3-$2;print $0,a;}' $file1 | awk '{sum+=$9} END {print sum}'`
# 	for j in {1..5}; do
# 		file2=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/atac.quintile$j.txt
# 		cov_gene=`bedtools intersect -a $file1 -b $file2 -wo | awk '{sum+=$19} END {print sum}'`
# 		gene=`awk '{a=$3-$2;print $0,a;}' $file2 | awk '{sum+=$11} END {print sum}'`
# 		echo "Total coverage length for genes in" $file1 ":" $cov_gene
# 		echo "Total feature length for genes:" $gene
# 		echo ""
# 		eval "samplecluster$j=`echo "scale=3; $cov_gene / $ip_length" | bc -l`"
# 		eval "featurecluster$j=`echo "scale=3; $cov_gene / $gene" | bc -l`"
# 	done
# 	NAME=`echo "Repeat"$i`
# 	echo -e "$NAME \t $samplecluster1 \t $samplecluster2 \t $samplecluster3 \t $samplecluster4 \t $samplecluster5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/sample_overlap_summary.atac.quintiles.cluster.repeat.txt
# 	echo -e "$NAME \t $featurecluster1 \t $featurecluster2 \t $featurecluster3 \t $featurecluster4 \t $featurecluster5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/feature_overlap_summary.atac.quintiles.cluster.repeat.txt
# done


# ##Overlap between peaks of different marks
# ## Ratio of mark (rowname) in (colname) (overlap divided by total length of rowname)
# # Summarize feature overlap
# ls *bedgraph | sed 's/_Tak1v4.regions.narrow.bedgraph//g' | tr '\n' '\t' > sample_overlap_summary.peaks.txt
# echo "" >> sample_overlap_summary.peaks.txt
# for f in *bedgraph; do
# 	NAME=`basename $f _Tak1v4.regions.narrow.bedgraph`
# 	eval "output1=`echo $NAME`"
# 	for file in *bedgraph; do
# 		f_ip_length=`awk '{a=$3-$2;print $0,a;}' $f | awk '{sum+=$5} END {print sum}'`
# 		cov_gene=`bedtools intersect -a $f -b $file -wo | awk '{sum+=$9} END {print sum}'`
# 		echo "Total coverage length for " $f " in " $file ":" $cov_gene
# 		echo ""
# 		output1+="boop"
# 		eval "output1+=`echo "scale=3; $cov_gene / $f_ip_length" | bc -l`"
# 	done
# 	echo $output1 | sed 's/boop/	/g' >> sample_overlap_summary.peaks.txt
# done
# sed -i 's/1.000/0.000/g' sample_overlap_summary.peaks.txt


# # ##Distance between features in certain clusters
# # ##1. Determine repeat clusters
# # # grep 'cluster 1' selected.repeat.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.repeat.TSS.profiledata.k5.heatmap.cluster1.txt
# # # grep 'cluster 2' selected.repeat.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.repeat.TSS.profiledata.k5.heatmap.cluster2.txt
# # # grep 'cluster 3' selected.repeat.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.repeat.TSS.profiledata.k5.heatmap.cluster3.txt
# # # grep 'cluster 4' selected.repeat.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.repeat.TSS.profiledata.k5.heatmap.cluster4.txt
# # # grep 'cluster 5' selected.repeat.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.repeat.TSS.profiledata.k5.heatmap.cluster5.txt
# # # Find average distance of nearest gene in specified cluster to repeat in specified cluster
# # echo -e "Cluster \t Gene_1 \t Gene_2 \t Gene_3 \t Gene_4 \t Gene_5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.txt
# # echo -e "Cluster \t Gene_1 \t Gene_2 \t Gene_3 \t Gene_4 \t Gene_5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/gene_to_repeat_distance.clusters.txt
# # for i in {1..5}; do
# # 	for j in {1..5}; do
# # 		file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.cluster$i.bed.txt
# # 		file2=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.cluster$j.bed.txt
# # 		eval "distance$j=`bedtools closest -a $file1 -b $file2 -d | awk '{ total += $17 } END { print total/NR }'`"
# # 		eval "distancegene$j=`bedtools closest -a $file2 -b $file1 -d | awk '{ total += $17 } END { print total/NR }'`"
# # 	done
# # 	echo -e "Repeat_$i \t $distance1 \t $distance2 \t $distance3 \t $distance4 \t $distance5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.txt
# # 	echo -e "Repeat_$i \t $distancegene1 \t $distancegene2 \t $distancegene3 \t $distancegene4 \t $distancegene5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/gene_to_repeat_distance.clusters.txt
# # done
# echo -e "Repeat1.Gene1 \t Repeat1.Gene2 \t Repeat1.Gene3 \t Repeat1.Gene4 \t Repeat1.Gene5 \t Repeat2.Gene1 \t Repeat2.Gene2 \t Repeat2.Gene3 \t Repeat2.Gene4 \t Repeat2.Gene5 \t Repeat3.Gene1 \t Repeat3.Gene2 \t Repeat3.Gene3 \t Repeat3.Gene4 \t Repeat3.Gene5 \t Repeat4.Gene1 \t Repeat4.Gene2 \t Repeat4.Gene3 \t Repeat4.Gene4 \t Repeat4.Gene5 \t Repeat5.Gene1 \t Repeat5.Gene2 \t Repeat5.Gene3 \t Repeat5.Gene4 \t Repeat5.Gene5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.all.txt
# echo -e "Gene1.Repeat1 \t Gene1.Repeat2 \t Gene1.Repeat3 \t Gene1.Repeat4 \t Gene1.Repeat5 \t Gene2.Repeat1 \t Gene2.Repeat2 \t Gene2.Repeat3 \t Gene2.Repeat4 \t Gene2.Repeat5 \t Gene3.Repeat1 \t Gene3.Repeat2 \t Gene3.Repeat3 \t Gene3.Repeat4 \t Gene3.Repeat5 \t Gene4.Repeat1 \t Gene4.Repeat2 \t Gene4.Repeat3 \t Gene4.Repeat4 \t Gene4.Repeat5 \t Gene5.Repeat1 \t Gene5.Repeat2 \t Gene5.Repeat3 \t Gene5.Repeat4 \t Gene5.Repeat5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/gene_to_repeat_distance.clusters.all.txt
# for i in {1..5}; do
# 	for j in {1..5}; do
# 		file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.cluster$i.bed.txt
# 		file2=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.cluster$j.bed.txt
# 		bedtools closest -a $file1 -b $file2 -d | cut -f15 > repeat$i.gene$j
# 		bedtools closest -a $file2 -b $file1 -d | cut -f15 > gene$j.repeat$i
# 	done
# done 
# paste repeat1.gene1 repeat1.gene2 repeat1.gene3 repeat1.gene4 repeat1.gene5 repeat2.gene1 repeat2.gene2 repeat2.gene3 repeat2.gene4 repeat2.gene5 repeat3.gene1 repeat3.gene2 repeat3.gene3 repeat3.gene4 repeat3.gene5 repeat4.gene1 repeat4.gene2 repeat4.gene3 repeat4.gene4 repeat4.gene5 repeat5.gene1 repeat5.gene2 repeat5.gene3 repeat5.gene4 repeat5.gene5 >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.all.txt
# paste gene1.repeat1 gene1.repeat2 gene1.repeat3 gene1.repeat4 gene1.repeat5 gene2.repeat1 gene2.repeat2 gene2.repeat3 gene2.repeat4 gene2.repeat5 gene3.repeat1 gene3.repeat2 gene3.repeat3 gene3.repeat4 gene3.repeat5 gene4.repeat1 gene4.repeat2 gene4.repeat3 gene4.repeat4 gene4.repeat5 gene5.repeat1 gene5.repeat2 gene5.repeat3 gene5.repeat4 gene5.repeat5 >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/gene_to_repeat_distance.clusters.all.txt

# ##Number of repeat class in a cluster
# ##1. Determine repeat clusters
# echo -e "Cluster \t Unknown \t LTR \t DNA \t LINE \t Helitron" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratioperclass.txt
# echo -e "Cluster \t Unknown \t LTR \t DNA \t LINE \t Helitron" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratiopercluster.txt
# file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.ids.txt
# unknown_total=`grep "Unknown" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# ltr_total=`grep "LTR" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# dna_total=`grep "DNA" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# line_total=`grep "LINE" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# helitron_total=`grep "Helitron" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# for i in {1..5}; do
# 	cluster_total=`grep "cluster $i" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	unknown_length=`grep "Unknown" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	ltr_length=`grep "LTR" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	dna_length=`grep "DNA" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	line_length=`grep "LINE" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	helitron_length=`grep "Helitron" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	eval "unknown_ratio=`echo "scale=3; $unknown_length / $unknown_total" | bc -l`"
# 	eval "ltr_ratio=`echo "scale=3; $ltr_length / $ltr_total" | bc -l`"
# 	eval "dna_ratio=`echo "scale=3; $dna_length / $dna_total" | bc -l`"
# 	eval "line_ratio=`echo "scale=3; $line_length / $line_total" | bc -l`"
# 	eval "helitron_ratio=`echo "scale=3; $helitron_length / $helitron_total" | bc -l`"
# 	echo -e "Repeat_$i \t $unknown_ratio \t $ltr_ratio \t $dna_ratio \t $line_ratio \t $helitron_ratio" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratioperclass.txt
# 	eval "unknown_ratio=`echo "scale=3; $unknown_length / $cluster_total" | bc -l`"
# 	eval "ltr_ratio=`echo "scale=3; $ltr_length / $cluster_total" | bc -l`"
# 	eval "dna_ratio=`echo "scale=3; $dna_length / $cluster_total" | bc -l`"
# 	eval "line_ratio=`echo "scale=3; $line_length / $cluster_total" | bc -l`"
# 	eval "helitron_ratio=`echo "scale=3; $helitron_length / $cluster_total" | bc -l`"
# 	echo -e "Repeat_$i \t $unknown_ratio \t $ltr_ratio \t $dna_ratio \t $line_ratio \t $helitron_ratio" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratiopercluster.txt
# done

# ##TE superfamilies by repeat cluster overlap
# echo -e "Repeat \t Cluster1 \t Cluster2 \t Cluster3 \t Cluster4 \t Cluster5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratioperclass.all.txt
# echo -e "Repeat \t Cluster1 \t Cluster2 \t Cluster3 \t Cluster4 \t Cluster5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratiopercluster.all.txt
# file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.ids.txt
# while read p; do
# 	file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.ids.txt
# 	repeat_total=`grep "$p" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	for i in {1..5}; do
# 		cluster_total=`grep "cluster $i" $file1 | grep -v "Unknown" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 		repeat_length=`grep "$p" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 		eval "ratioperclass$i=`echo "scale=3; $repeat_length / $repeat_total" | bc -l`"
# 		eval "ratiopercluster$i=`echo "scale=3; $repeat_length / $cluster_total" | bc -l`"
# 	done
# 	echo -e "$p \t $ratioperclass1 \t $ratioperclass2 \t $ratioperclass3 \t $ratioperclass4 \t $ratioperclass5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratioperclass.all.txt
# 	echo -e "$p \t $ratiopercluster1 \t $ratiopercluster2 \t $ratiopercluster3 \t $ratiopercluster4 \t $ratiopercluster5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratiopercluster.all.txt	
# done </volumes/berger/lab/Marchantia_v4.4/repeat.superfamilies.txt 





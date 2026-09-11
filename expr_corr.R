#correlate gene expression with clusters defined in deepTools for C&R data

##Merge TPM from each file 
# tail -n +2 spore-25960/rsem/spore-25960.genes.results | cut -f1 | sed 's/_Mp.*//g' > TPM.txt
# awk -v OFS='\t' '{_[FNR]=(_[FNR] OFS $6)}END{for (i=1; i<=FNR; i++) {sub(/\t/,"",_[i]); print _[i]}}' */rsem/*genes.results | tail -n +2 | paste TPM.txt - > TPM2.txt
# mv TPM2.txt TPM.txt
# ls <copy result>
# vi TPM.txt <paste result>

setwd("/Volumes/berger/user/sean.montgomery/Documents/rna-seq/Tak1v4/higo2016/expression")

tpm <- read.table("TPM.txt", header = T)
gametop <- subset(tpm, select = c("g.1","g.2","g.3"))
rownames(gametop) <- tpm[,1]
gametop$mean <- (rowSums(gametop)/3)
antheridia <- subset(tpm, select = c("a.1","a.2"))
rownames(antheridia) <- tpm[,1]
antheridia$mean <- (rowSums(antheridia)/2)
antheridiophore <- subset(tpm, select = c("an.1","an.2","an.3"))
rownames(antheridiophore) <- tpm[,1]
antheridiophore$mean <- (rowSums(antheridiophore)/2)
archegoniophore <- subset(tpm, select = c("ar.1","ar.2","ar.3"))
rownames(archegoniophore) <- tpm[,1]
archegoniophore$mean <- (rowSums(archegoniophore)/2)

tpmfrank <- read.table("/Volumes/berger/user/sean.montgomery/Documents/rna-seq/Tak1v4/frank2014/expression/tpm.txt", header = T)
frankgametop <- subset(tpmfrank, select = c("g.1"))
rownames(frankgametop) <- tpmfrank[,1]
frankgametop$mean <- (rowSums(frankgametop)/1)
sporop <- subset(tpmfrank, select = c("sp.1","sp.2","sp.3"))
rownames(sporop) <- tpmfrank[,1]
sporop$mean <- (rowSums(sporop)/3)
tip <- subset(tpmfrank, select = c("tip.1","tip.2","tip.3"))
rownames(tip) <- tpmfrank[,1]
tip$mean <- (rowSums(tip)/3)

# clusters <- read.table("../../cutrun/Tak1v3/genes.TSS.k5.black.profiledata.txt")
# clust <- clusters[order(clusters$V4),]
# adder <- c()
# adder[which(clust$V4 %in% row.names(gametop))] <- as.numeric(gametop$X23955[which(row.names(gametop) %in% clust$V4)])
# adder1 <- c()
# adder1[which(clust$V4 %in% row.names(gametop))] <- as.numeric(gametop$X23956[which(row.names(gametop) %in% clust$V4)])
# adder2 <- c()
# adder2[which(clust$V4 %in% row.names(gametop))] <- as.numeric(gametop$X24193[which(row.names(gametop) %in% clust$V4)])
# clust0 <- cbind(clust, rep1=c(adder), rep2=c(adder1), rep3=c(adder2))
# clust0$V8 <- as.factor(clust0$V8)

# boxplot(rep1~V8,data=clust0)
# boxplot(rep1~V8,data=clust0,outline=FALSE)
# boxplot(rep2~V8,data=clust0,outline=FALSE)
# boxplot(rep3~V8,data=clust0,outline=FALSE)



#Plot gene expression for each mark; possibly by quintiles
#Creat +1 -1 nucleosome gene annotations
genes=/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.gene.bed
# newgenes=/volumes/berger/lab/NGS_annotations/backup_berger_common/MpTak1v4.4/annotations/MpTak1v4.4.-+1.bed
# sortednewgenes=/volumes/berger/lab/NGS_annotations/backup_berger_common/MpTak1v4.4/annotations/MpTak1v4.4.-+1.sorted.bed
# awk '{                                                                                                      
# if ($4 == "+")
# print $1,$2+200,$3,$4,$5;
# else if($4 == "-")
# print $1,$2,$3-200,$4,$5;
# }' OFS="\t" $genes > $newgenes
# awk '{                                                                                                      
# if ($2>$3)
# print $1,$3,$2,$4,$5;
# else if ($2<$3)
# print $1,$2,$3,$4,$5;
# else if ($2==$3)
# print $1,$2,$3,$4,$5;
# }' OFS="\t" $newgenes > $sortednewgenes
for f in *narrow.bedgraph; do bedtools intersect -wa -u -F 5E-1 -a $genes -b $f > genes/$f.5e-1.bed; done
for f in genes/*5e-1.bed; do cut -f4 $f > $f.tmp; done
paste genes/*.tmp > genes/genelist2.5e-1.txt 
rm genes/*.tmp
ll genes/ | tail -n +2 | tr -s ' ' | cut -d ' ' -f9 | sed 's/_Tak1v4.regions.narrow.bedgraph.5e-1.bed//g' | tr '\n' '\t'

genelist <- as.matrix(read.table("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/genes/genelist2.5e-1.txt",header = T, fill =T, sep='\t'))
a1 <- as.numeric(gametop$mean[pmatch(genelist[,1],row.names(gametop))]) #FLAG-DUO1
a2 <- as.numeric(gametop$mean[pmatch(genelist[,2],row.names(gametop))]) #FLAG-WT
a3 <- as.numeric(gametop$mean[pmatch(genelist[,3],row.names(gametop))]) #H2A
a4 <- as.numeric(gametop$mean[pmatch(genelist[,4],row.names(gametop))]) #H2A
a5 <- as.numeric(gametop$mean[pmatch(genelist[,5],row.names(gametop))]) #H2A
a6 <- as.numeric(gametop$mean[pmatch(genelist[,6],row.names(gametop))]) #H2A.W.2
a7 <- as.numeric(gametop$mean[pmatch(genelist[,7],row.names(gametop))]) #H2A.W.2
a8 <- as.numeric(gametop$mean[pmatch(genelist[,8],row.names(gametop))]) #H2A.W.2
a9 <- as.numeric(gametop$mean[pmatch(genelist[,9],row.names(gametop))]) #H2A.W.2
a10 <- as.numeric(gametop$mean[pmatch(genelist[,10],row.names(gametop))]) #H2A.W.2
a11 <- as.numeric(gametop$mean[pmatch(genelist[,11],row.names(gametop))]) #H2A.W.2
a12 <- as.numeric(gametop$mean[pmatch(genelist[,12],row.names(gametop))]) #H2A.X
a13 <- as.numeric(gametop$mean[pmatch(genelist[,13],row.names(gametop))]) #H2A.X
a14 <- as.numeric(gametop$mean[pmatch(genelist[,14],row.names(gametop))]) #H2A.X
a15 <- as.numeric(gametop$mean[pmatch(genelist[,15],row.names(gametop))]) #H2A.Z
a16 <- as.numeric(gametop$mean[pmatch(genelist[,16],row.names(gametop))]) #H2A.Z
a17 <- as.numeric(gametop$mean[pmatch(genelist[,17],row.names(gametop))]) #H3
a18 <- as.numeric(gametop$mean[pmatch(genelist[,18],row.names(gametop))]) #H3K14ac
a19 <- as.numeric(gametop$mean[pmatch(genelist[,19],row.names(gametop))]) #H3K14ac
a20 <- as.numeric(gametop$mean[pmatch(genelist[,20],row.names(gametop))]) #H3K27me1
a21 <- as.numeric(gametop$mean[pmatch(genelist[,21],row.names(gametop))]) #H3K27me1
a22 <- as.numeric(gametop$mean[pmatch(genelist[,22],row.names(gametop))]) #H3K27me3
a23 <- as.numeric(gametop$mean[pmatch(genelist[,23],row.names(gametop))]) #H3K27me3
a24 <- as.numeric(gametop$mean[pmatch(genelist[,24],row.names(gametop))]) #H3K27me3
a25 <- as.numeric(gametop$mean[pmatch(genelist[,25],row.names(gametop))]) #H3K36me3
a26 <- as.numeric(gametop$mean[pmatch(genelist[,26],row.names(gametop))]) #H3K36me3
a27 <- as.numeric(gametop$mean[pmatch(genelist[,27],row.names(gametop))]) #H3K36me3
a28 <- as.numeric(gametop$mean[pmatch(genelist[,28],row.names(gametop))]) #H3K36me3
a29 <- as.numeric(gametop$mean[pmatch(genelist[,29],row.names(gametop))]) #H3K36me3dig
a30 <- as.numeric(gametop$mean[pmatch(genelist[,30],row.names(gametop))]) #H3K36me3dig
a31 <- as.numeric(gametop$mean[pmatch(genelist[,31],row.names(gametop))]) #H3K36me3dig
a32 <- as.numeric(gametop$mean[pmatch(genelist[,32],row.names(gametop))]) #H3K36me3dig
a33 <- as.numeric(gametop$mean[pmatch(genelist[,33],row.names(gametop))]) #H3K36me3dig
a34 <- as.numeric(gametop$mean[pmatch(genelist[,34],row.names(gametop))]) #H3K4me1
a35 <- as.numeric(gametop$mean[pmatch(genelist[,35],row.names(gametop))]) #H3K4me1
a36 <- as.numeric(gametop$mean[pmatch(genelist[,36],row.names(gametop))]) #H3K4me3
a37 <- as.numeric(gametop$mean[pmatch(genelist[,37],row.names(gametop))]) #H3K4me3
a38 <- as.numeric(gametop$mean[pmatch(genelist[,38],row.names(gametop))]) #H3K9ac
a39 <- as.numeric(gametop$mean[pmatch(genelist[,39],row.names(gametop))]) #H3K9ac
a40 <- as.numeric(gametop$mean[pmatch(genelist[,40],row.names(gametop))]) #H3K9acK14ac
a41 <- as.numeric(gametop$mean[pmatch(genelist[,41],row.names(gametop))]) #H3K9me1
a42 <- as.numeric(gametop$mean[pmatch(genelist[,42],row.names(gametop))]) #H3K9me1
a43 <- as.numeric(gametop$mean[pmatch(genelist[,43],row.names(gametop))]) #H3K9me2
a44 <- as.numeric(gametop$mean[pmatch(genelist[,44],row.names(gametop))]) #H3K9me2
a45 <- as.numeric(gametop$mean[pmatch(genelist[,45],row.names(gametop))]) #H3K9me2
a46 <- as.numeric(gametop$mean[pmatch(genelist[,46],row.names(gametop))]) #H3K9me2
a47 <- as.numeric(gametop$mean[pmatch(genelist[,47],row.names(gametop))]) #H3K9me2
a48 <- as.numeric(gametop$mean[pmatch(genelist[,48],row.names(gametop))]) #H3K9me3
a49 <- as.numeric(gametop$mean[pmatch(genelist[,49],row.names(gametop))]) #H4K20me1
a50 <- as.numeric(gametop$mean[pmatch(genelist[,50],row.names(gametop))]) #H4K20me1
a51 <- as.numeric(gametop$mean[pmatch(genelist[,51],row.names(gametop))]) #IgG
a52 <- as.numeric(gametop$mean[pmatch(genelist[,52],row.names(gametop))]) #SpH2A.W.1
a53 <- as.numeric(gametop$mean[pmatch(genelist[,53],row.names(gametop))]) #SpH2A.X.1
a54 <- as.numeric(gametop$mean[pmatch(genelist[,54],row.names(gametop))]) #SpH2A.Z
a55 <- as.numeric(gametop$mean[pmatch(genelist[,55],row.names(gametop))]) #SpH3
a56 <- as.numeric(gametop$mean[pmatch(genelist[,56],row.names(gametop))]) #SpH3
a57 <- as.numeric(gametop$mean[pmatch(genelist[,57],row.names(gametop))]) #SpH3K27me1
a58 <- as.numeric(gametop$mean[pmatch(genelist[,58],row.names(gametop))]) #SpH3K27me3
a59 <- as.numeric(gametop$mean[pmatch(genelist[,59],row.names(gametop))]) #SpH3K27me3
a60 <- as.numeric(gametop$mean[pmatch(genelist[,60],row.names(gametop))]) #SpH3K27me3
a61 <- as.numeric(gametop$mean[pmatch(genelist[,61],row.names(gametop))]) #SpH3K36me3
a62 <- as.numeric(gametop$mean[pmatch(genelist[,62],row.names(gametop))]) #SpH3K36me3
a63 <- as.numeric(gametop$mean[pmatch(genelist[,63],row.names(gametop))]) #SpH3K36me3
a64 <- as.numeric(gametop$mean[pmatch(genelist[,64],row.names(gametop))]) #SpH3K36me3
a65 <- as.numeric(gametop$mean[pmatch(genelist[,65],row.names(gametop))]) #SpH3K4me1
a66 <- as.numeric(gametop$mean[pmatch(genelist[,66],row.names(gametop))]) #SpH3K9me2
a67 <- as.numeric(gametop$mean[pmatch(genelist[,67],row.names(gametop))]) #SpH3K9me3
a68 <- as.numeric(gametop$mean[pmatch(genelist[,68],row.names(gametop))]) #cenH3
b25 <- as.numeric(gametop$mean[pmatch(genelist[which(genelist[,25] %in% h3k27me3.1$V10),25],row.names(gametop))]) #H3K4me3 + H3K27me3
b26 <- as.numeric(gametop$mean[pmatch(genelist[which(genelist[,26] %in% h3k27me3.1$V10),26],row.names(gametop))]) #H3K4me3 + H3K27me3
c25 <- as.numeric(gametop$mean[pmatch(genelist[which(!genelist[,25] %in% h3k27me3.1$V10),25],row.names(gametop))]) #H3K4me3 - H3K27me3
c26 <- as.numeric(gametop$mean[pmatch(genelist[which(!genelist[,26] %in% h3k27me3.1$V10),26],row.names(gametop))]) #H3K4me3 - H3K27me3

# par(mar = c(7,4,2,4))
# boxplot(a3,a8,a13,a15,a22,a42,a20,a25,a18,a34,a36,a38, 
#         outline=F, names = c(colnames(genelist)[3],colnames(genelist)[8],colnames(genelist)[13],colnames(genelist)[15],
#         colnames(genelist)[22],colnames(genelist)[42],colnames(genelist)[20],
#         colnames(genelist)[25],colnames(genelist)[18],colnames(genelist)[34],colnames(genelist)[36],colnames(genelist)[38]), 
#         las=2, ylim = c(0,10), 
#         ylab="Gene expression (asinh(TPM))", 
#         col=c("light blue","light blue","light blue","light blue",
#         "red","red","red",
#         "light green","light green","light green","light green","light green"))
# boxplot(a1,a2, a3, a4, a5, a6, a7, a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,
# 	a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,
# 	a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68, 
# 	outline=F, names = colnames(genelist[,1:68]), las=2, ylab="Gene expression (TPM)", 
# 	col=c("black","black","black","green","green","green","green","green","grey","red","red","orange",
# 		"orange","green","green","green","green","green","green","green","green","green","green","green",
# 		"green","green","green","red","red","black","black","black","black","grey","black","black","black",
# 		"black","black","black","black","black","black","black","black","black","black","black","purple"))

## Used for figure
df <- data.frame(asinh(a9), asinh(a3), asinh(a8), asinh(a2), asinh(a6), asinh(a5), asinh(a7), asinh(a1), asinh(a4), asinh(a10))
colnames(df) <- colnames(genelist)[c(9,3,8,2,6,5,7,1,4,10)]
df.m <- reshape2::melt(df, id.vars = NULL)
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/expr.permark.pdf",height=6.3,width=10)
ggplot(df.m, aes(x = variable, y = value, color=variable)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Histone modification")+scale_y_continuous(name="asinh(TPM)")
dev.off()


##Expr per k means cluster
# bedtools intersect -wo -F 1 -f 1 -a $genes -b /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.txt | cut -f1,2,3,4,5,12 | uniq > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.ids.txt
##Create separate bed files of genes per cluster
## for i in {1..5}; do grep "cluster $i" selected.TSS.profiledata.k5.ids.txt | cut -f1,2,3,4,5 > selected.TSS.profiledata.k5.ids.cluster$i.bed; done
k4profile <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.ids.txt", sep = '\t')
k4profile$TPM <- gametop$mean[which(rownames(gametop) %in% k4profile$V4)]
k4profile$asinhTPM <- asinh(k4profile$TPM)
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/expr.clusters.k5.gametophyte.pdf",height=6.3,width=10)
ggplot(k4profile, aes(x = V6, y = asinhTPM, color=V6)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="asinh(TPM)",limits = c(0,11))
dev.off()
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/expr.clusters.k5.archegoniophore.pdf",height=6.3,width=10)
k4profile$TPM <- archegoniophore$mean[which(rownames(archegoniophore) %in% k4profile$V4)]
k4profile$asinhTPM <- asinh(k4profile$TPM)
ggplot(k4profile, aes(x = V6, y = asinhTPM, color=V6)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="asinh(TPM)",limits = c(0,11))
dev.off()
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/expr.clusters.k5.antheridiophore.pdf",height=6.3,width=10)
k4profile$TPM <- antheridiophore$mean[which(rownames(antheridiophore) %in% k4profile$V4)]
k4profile$asinhTPM <- asinh(k4profile$TPM)
ggplot(k4profile, aes(x = V6, y = asinhTPM, color=V6)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="asinh(TPM)",limits = c(0,11))
dev.off()
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/expr.clusters.k5.antheridia.pdf",height=6.3,width=10)
k4profile$TPM <- antheridia$mean[which(rownames(antheridia) %in% k4profile$V4)]
k4profile$asinhTPM <- asinh(k4profile$TPM)
ggplot(k4profile, aes(x = V6, y = asinhTPM, color=V6)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="asinh(TPM)",limits = c(0,11))
dev.off()
##For Frank2014 data
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/expr.clusters.k5.frankgametop.pdf",height=6.3,width=10)
k4profile$TPM <- frankgametop$mean[which(rownames(frankgametop) %in% k4profile$V4)]
k4profile$asinhTPM <- asinh(k4profile$TPM)
ggplot(k4profile, aes(x = V6, y = asinhTPM, color=V6)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="asinh(TPM)",limits = c(0,11))
dev.off()
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/expr.clusters.k5.sporop.pdf",height=6.3,width=10)
k4profile$TPM <- sporop$mean[which(rownames(sporop) %in% k4profile$V4)]
k4profile$asinhTPM <- asinh(k4profile$TPM)
ggplot(k4profile, aes(x = V6, y = asinhTPM, color=V6)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="asinh(TPM)",limits = c(0,11))
dev.off()
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/expr.clusters.k5.tip.pdf",height=6.3,width=10)
k4profile$TPM <- tip$mean[which(rownames(tip) %in% k4profile$V4)]
k4profile$asinhTPM <- asinh(k4profile$TPM)
ggplot(k4profile, aes(x = V6, y = asinhTPM, color=V6)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="asinh(TPM)",limits = c(0,11))
dev.off()

##Length of features in cluster
k4profile$length <- k4profile$V3 - k4profile$V2
k4profile$loglength <- log10(k4profile$length)
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/length.genes.clusters.k5.pdf",height=6.3,width=10)
ggplot(k4profile, aes(x = V6, y = loglength, color=V6)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Length (log(bp))")
dev.off()
## For repeats
# bedtools intersect -wo -F 1 -f 1 -a /volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.repeat.bed -b /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.txt | cut -f1,2,3,4,5,12 | uniq > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.ids.txt
k5repeat <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.ids.txt", sep = '\t')
k5repeat$length <- k5repeat$V3 - k5repeat$V2
k5repeat$loglength <- log10(k5repeat$length)
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/length.repeats.clusters.k5.pdf",height=6.3,width=10)
ggplot(k5repeat, aes(x = V6, y = loglength, color=V6)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Length (log(bp))")
dev.off()

##Expr per quintile
quintiles <- quantile(gametop$mean, probs = seq(0,1,0.2))
bedlist <- read.table("/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.gene.bed")
bedlist$TPM <- gametop$mean
bedlist$asinhTPM <- asinh(bedlist$TPM)
bedlist$quintile <- cut(bedlist$TPM, breaks = c(-Inf, quintiles[2], quintiles[3], quintiles[4],quintiles[5],Inf),labels = c(1,2,3,4,5))
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/expr.quintiles.pdf",height=6.3,width=10)
ggplot(bedlist, aes(x = quintile, y = asinhTPM, color=quintile)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Quintile")+scale_y_continuous(name="asinh(TPM)")
dev.off()
##Length of features in quintile
bedlist$length <- bedlist$V3 - bedlist$V2
bedlist$loglength <- log10(bedlist$length)
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/length.genes.quintiles.k5.pdf",height=6.3,width=10)
ggplot(bedlist, aes(x = quintile, y = loglength, color=quintile)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Quintile")+scale_y_continuous(name="Length (log(bp))")
dev.off()

##DNA methylation level of peaks from marks
##Table with columns of 1)Mark, 2)Methylation type, 3)Methylation % for each peak (all together; super long table...)
# for f in /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/markpeaks/Thallus-Male_*; do MET=$(echo $f | cut -f2 -d '_' | cut -f1 -d '.'); MARK=$(echo $f | cut -f2 -d '_' | cut -f2 -d '.' | sed 's/H2A/H2A.Z/g'); cut -f5 $f | while read line; do echo "$MET	$MARK	$line" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/markpeaks/SRR5314038_Thallus-Male_2.met.marks.txt; done; done
for f in /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/markpeaks/SRR5314038_Thallus-Male_2_*; do MET=$(echo $f | cut -f4 -d '_' | cut -f1 -d '.'); MARK=$(echo $f | cut -f4 -d '_' | cut -f2 -d '.' | sed 's/H2A/H2A.Z/g'); cut -f5 $f | while read line; do echo "$MET	$MARK	$line" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/markpeaks/SRR5314038_Thallus-Male_2.met.marks.txt; done; done
##Grouped violin plots per mark of CG, CHG, CHH
##New colour code for different methylation types
metmark <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/markpeaks/SRR5314038_Thallus-Male_2.met.marks.txt", sep= '\t')
cyl_table <- table(metmark$V2)
cyl_levels <- names(cyl_table)[order(cyl_table)]
metmark$cyl2 <- factor(metmark$V2, levels = c("H3K9me1","H3K27me1","H3K9ac","H3K14ac","H3K4me1","H3K36me3","H3K4me3","H2A.Z","H3K27me3","H3"))
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.mark.SRR5314038_thallus-male_2.pdf",height=6.3,width=10)
ggplot(data=subset(metmark, !is.na(cyl2)) ,aes(x=cyl2, y=V3, color=V1)) + geom_violin(scale = "count",width=1,position = position_dodge(width=0.7)) + stat_summary(aes(group=V1),fun.y=median, geom="point", size=2, color="red", position = position_dodge(width = 0.7)) + scale_color_manual(values=c("#6C756B","#C59B76","#ADB2D3"))+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Mark")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,10), ylim = c(0,100))
dev.off()
##Methylation level of repeats by cluster
##Write a shortcut method for this...
dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/SRR5314038_Thallus-Male_2_CpG.repeat.cluster.bed", sep = '\t')
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CG.SRR5314038_Thallus-Male_2.pdf",height=6.3,width=10)
ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
dev.off()
dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/SRR5314038_Thallus-Male_2_CHG.repeat.cluster.bed", sep = '\t')
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHG.SRR5314038_Thallus-Male_2.pdf",height=6.3,width=10)
ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
dev.off()
dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/SRR5314038_Thallus-Male_2_CHH.repeat.cluster.bed", sep = '\t')
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHH.SRR5314038_Thallus-Male_2.pdf",height=6.3,width=10)
ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/Thallus-Male_CpG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CG.thallus-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/Thallus-Male_CHG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHG.thallus-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/Thallus-Male_CHH.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHH.thallus-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/Thallus-Female_CHH.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHH.thallus-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/Thallus-Female_CHG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHG.thallus-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/Thallus-Female_CpG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CG.thallus-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/AntheridiaSperm_CpG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CG.antheridia-sperm.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/AntheridiaSperm_CHG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHG.antheridia-sperm.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/AntheridiaSperm_CHH.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHH.antheridia-sperm.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/ApicalNotch-Female_CHH.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHH.apicalnotch-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/ApicalNotch-Female_CHG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHG.apicalnotch-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/ApicalNotch-Female_CpG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CG.apicalnotch-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/ApicalNotch-Male_CpG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CG.apicalnotch-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/ApicalNotch-Male_CHG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHG.apicalnotch-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/ApicalNotch-Male_CHH.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHH.apicalnotch-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/Archegonia_CHH.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHH.archegonia.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/Archegonia_CHG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHG.archegonia.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/Archegonia_CpG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CG.archegonia.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/EarlySporophyte_CpG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CG.earlysporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/EarlySporophyte_CHG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHG.earlysporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/EarlySporophyte_CHH.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHH.earlysporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/LateSporophyte_CHH.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHH.latesporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/LateSporophyte_CHG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CHG.latesporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/repeatcluster/LateSporophyte_CpG.repeat.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.repeat.CG.latesporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=2) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
##Methylation level of genes by cluster
dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/SRR5314038_Thallus-Male_2_CpG.gene.cluster.bed", sep = '\t')
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CG.SRR5314038_Thallus-Male_2.pdf",height=6.3,width=10)
ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
dev.off()
dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/SRR5314038_Thallus-Male_2_CHG.gene.cluster.bed", sep = '\t')
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHG.SRR5314038_Thallus-Male_2.pdf",height=6.3,width=10)
ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
dev.off()
dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/SRR5314038_Thallus-Male_2_CHH.gene.cluster.bed", sep = '\t')
pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHH.SRR5314038_Thallus-Male_2.pdf",height=6.3,width=10)
ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/Thallus-Male_CpG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CG.thallus-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/Thallus-Male_CHG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHG.thallus-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/Thallus-Male_CHH.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHH.thallus-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/Thallus-Female_CHH.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHH.thallus-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/Thallus-Female_CHG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHG.thallus-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/Thallus-Female_CpG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CG.thallus-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/AntheridiaSperm_CpG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CG.antheridia-sperm.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/AntheridiaSperm_CHG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHG.antheridia-sperm.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/AntheridiaSperm_CHH.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHH.antheridia-sperm.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/ApicalNotch-Female_CHH.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHH.apicalnotch-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/ApicalNotch-Female_CHG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHG.apicalnotch-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/ApicalNotch-Female_CpG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CG.apicalnotch-female.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/ApicalNotch-Male_CpG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CG.apicalnotch-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/ApicalNotch-Male_CHG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHG.apicalnotch-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/ApicalNotch-Male_CHH.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHH.apicalnotch-male.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/Archegonia_CHH.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHH.archegonia.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/Archegonia_CHG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHG.archegonia.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/Archegonia_CpG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CG.archegonia.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/EarlySporophyte_CpG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CG.earlysporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/EarlySporophyte_CHG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHG.earlysporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/EarlySporophyte_CHH.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHH.earlysporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/LateSporophyte_CHH.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHH.latesporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/LateSporophyte_CHG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHG.latesporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genecluster/LateSporophyte_CpG.gene.cluster.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CG.latesporophyte.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V6, y = V7, color=V6)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genji/SRR5314038_Thallus-Male_2_CpG.TY_Genji.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CG.genji.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genji/SRR5314038_Thallus-Male_2_CHG.TY_Genji.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHG.genji.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/genji/SRR5314038_Thallus-Male_2_CHH.TY_Genji.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/dnamet.gene.CHH.genji.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/SRR5314038_Thallus-Male_2_CpG.Copia.heatmap.k5.clusters.sorted.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/dnamet.CG.copia.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#2171B5","#FB6A4A","#74C476","#FB6A4A","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/SRR5314038_Thallus-Male_2_CHG.Copia.heatmap.k5.clusters.sorted.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/dnamet.CHG.copia.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#2171B5","#FB6A4A","#74C476","#FB6A4A","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/SRR5314038_Thallus-Male_2_CHH.Copia.heatmap.k5.clusters.sorted.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/dnamet.CHH.copia.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#2171B5","#FB6A4A","#74C476","#FB6A4A","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,5), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/SRR5314038_Thallus-Male_2_CpG.Gypsy.heatmap.k3.clusters.sorted.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/dnamet.CG.gypsy.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/SRR5314038_Thallus-Male_2_CHG.Gypsy.heatmap.k3.clusters.sorted.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/dnamet.CHG.gypsy.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/SRR5314038_Thallus-Male_2_CHH.Gypsy.heatmap.k3.clusters.sorted.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/dnamet.CHH.gypsy.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/SRR5314038_Thallus-Male_2_CpG.HomeoGypsy.heatmap.k3.clusters.sorted.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/dnamet.CG.homeogypsy.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/SRR5314038_Thallus-Male_2_CHG.HomeoGypsy.heatmap.k3.clusters.sorted.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/dnamet.CHG.homeogypsy.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
# dev.off()
# dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/SRR5314038_Thallus-Male_2_CHH.HomeoGypsy.heatmap.k3.clusters.sorted.bed", sep = '\t')
# pdf("/Volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/dnamet.CHH.homeogypsy.pdf",height=6.3,width=10)
# ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
# dev.off()

library(reshape2)
dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/Hi-C/CL/Tak2.G2_G2UA.matrix.50kb.foci.LateSporophyte_CpG.bed", sep = '\t')
dnamet2 <- read.table("/volumes/berger/user/sean.montgomery/Documents/Hi-C/CL/Tak2.G2_G2UA.matrix.50kb.foci.LateSporophyte_CHG.bed", sep = '\t')
dnamet3 <- read.table("/volumes/berger/user/sean.montgomery/Documents/Hi-C/CL/Tak2.G2_G2UA.matrix.50kb.foci.LateSporophyte_CHH.bed", sep = '\t')
dnamet4 <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/Tak1v4.totalchr.LateSporophyte_CpG.bed", sep = '\t')
dnamet5 <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/Tak1v4.totalchr.LateSporophyte_CHG.bed", sep = '\t')
dnamet6 <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/dnamet/Tak1v4.totalchr.LateSporophyte_CHH.bed", sep = '\t')
dnamet$V6 <- dnamet2$V5
dnamet$V7 <- dnamet3$V5
dnamet4$V8 <- dnamet4$V5
dnamet4$V9 <- dnamet5$V5
dnamet4$V10 <- dnamet6$V5
pdf("/volumes/berger/user/sean.montgomery/Documents/Hi-C/CL/dnamet.Tak2.G2_G2UA.matrix.50kb.foci.LateSporophyte.pdf",height=6.3,width=10)
ggplot(data=rbind(melt(dnamet[,5:7]),melt(dnamet4[,6:8])), aes(x = variable, y = value, color=variable)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,6), ylim = c(0,100))
dev.off()

dnamet <- read.table("/volumes/berger/user/sean.montgomery/Documents/Hi-C/CL/Tak2.G2_G2UA.matrix.50kb.foci.SRR5314038_Thallus-Male_2_CpG.bed", sep = '\t')
pdf("/volumes/berger/user/sean.montgomery/Documents/Hi-C/CL/dnamet.Tak2.G2_G2UA.matrix.50kb.foci.CG.pdf",height=6.3,width=10)
ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
dev.off()
dnamet2 <- read.table("/volumes/berger/user/sean.montgomery/Documents/Hi-C/CL/Tak2.G2_G2UA.matrix.50kb.foci.SRR5314038_Thallus-Male_2_CHG.bed", sep = '\t')
pdf("/volumes/berger/user/sean.montgomery/Documents/Hi-C/CL/dnamet.Tak2.G2_G2UA.matrix.50kb.foci.CHG.pdf",height=6.3,width=10)
ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
dev.off()
dnamet3 <- read.table("/volumes/berger/user/sean.montgomery/Documents/Hi-C/CL/Tak2.G2_G2UA.matrix.50kb.foci.SRR5314038_Thallus-Male_2_CHH.bed", sep = '\t')
pdf("/volumes/berger/user/sean.montgomery/Documents/Hi-C/CL/dnamet.Tak2.G2_G2UA.matrix.50kb.foci.CHH.pdf",height=6.3,width=10)
ggplot(dnamet, aes(x = V4, y = V5, color=V4)) + geom_violin(scale = "count",width=1) + stat_summary(fun.y=median, geom="point", size=2, color="red") + scale_color_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"),guide=F)+theme(axis.text.x = element_text(color = "black",size = 10, angle = 45,hjust = 1))+scale_x_discrete(name="Cluster")+scale_y_continuous(name="Methylation (%)")+coord_cartesian(xlim = c(1,3), ylim = c(0,100))
dev.off()

##Horizontal scatterplot of TPM for heatmap
library(tidyr)
data_summary <- function(data, varname, groupnames){
	require(plyr)
	summary_func <- function(x, col){
		c(mean = mean(x[[col]], na.rm=TRUE),
		sd = sd(x[[col]], na.rm=TRUE))
	}
	data_sum<-ddply(data, groupnames, .fun=summary_func,
	varname)
	data_sum <- rename(data_sum, c("mean" = varname))
	return(data_sum)
}
geneorder <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.txt", sep = '\t')
clusterids <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.ids.txt", sep = '\t')
tpm <- read.table("/volumes/berger/user/sean.montgomery/Documents/rna-seq/Tak1v4/higo2016/expression/TPM.txt", header = T)
clusterids$g1 <- tpm$g.1[match(clusterids$V4, tpm$gene_id)]
clusterids$g2 <- tpm$g.2[match(clusterids$V4, tpm$gene_id)]
clusterids$g3 <- tpm$g.3[match(clusterids$V4, tpm$gene_id)]
# clusterids$tpm <- (clusterids$g1+clusterids$g2+clusterids$g3)/3
# geneorder$gene_id <- clusterids$V4[match(geneorder$V3, clusterids$V3)]
# geneorder$g1 <- tpm$g.1[match(geneorder$gene_id, tpm$gene_id)]
# geneorder$g2 <- tpm$g.2[match(geneorder$gene_id, tpm$gene_id)]
# geneorder$g3 <- tpm$g.3[match(geneorder$gene_id, tpm$gene_id)]
geneorder2 <- clusterids[order(-rank(clusterids$V6), clusterids$tpm, decreasing = T),]
geneorder.output <- geneorder2[,1:7]
write.table(geneorder.output, file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.tpmorder.txt", row.names = F, col.names = F, quote = F, sep = '\t')
geneorder2$V4 <- factor(geneorder2$V4, levels=unique(geneorder2$V4))
# df <- gather(geneorder2, key = replicate, tpm, g1:g3, factor_key = TRUE)
# df2 <- data_summary(df, varname = "tpm", groupnames = c("gene_id", "V7"))
ggplot(subset(geneorder2, geneorder2$V6 == "cluster 1"), aes(x=asinh(tpm), y=factor(V4, levels = rev(levels(factor(V4)))))) + geom_point(color="#C6DBEF") + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(limits = c(0,11)) + scale_y_discrete(name= element_blank())
ggplot(subset(geneorder2, geneorder2$V6 == "cluster 2"), aes(x=asinh(tpm), y=factor(V4, levels = rev(levels(factor(V4)))))) + geom_point(color="#A1D99B") + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(limits = c(0,11)) + scale_y_discrete(name= element_blank())
ggplot(subset(geneorder2, geneorder2$V6 == "cluster 3"), aes(x=asinh(tpm), y=factor(V4, levels = rev(levels(factor(V4)))))) + geom_point(color="#005A32") + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(limits = c(0,11)) + scale_y_discrete(name= element_blank())
ggplot(subset(geneorder2, geneorder2$V6 == "cluster 4"), aes(x=asinh(tpm), y=factor(V4, levels = rev(levels(factor(V4)))))) + geom_point(color="#4292C6") + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(limits = c(0,11)) + scale_y_discrete(name= element_blank())
ggplot(subset(geneorder2, geneorder2$V6 == "cluster 5"), aes(x=asinh(tpm), y=factor(V4, levels = rev(levels(factor(V4)))))) + geom_point(color="#CB181D") + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(limits = c(0,11)) + scale_y_discrete(name= element_blank())
# ggplot(subset(geneorder2, geneorder2$V7 == "cluster 1"), aes(x=asinh(tpm), y=factor(gene_id, levels = rev(levels(factor(gene_id)))))) + geom_point(color="#C6DBEF") + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(limits = c(-4,11)) + geom_errorbarh(aes(xmin=asinh(tpm-sd), xmax=asinh(tpm+sd)))
# ggplot(subset(geneorder2, geneorder2$V7 == "cluster 2"), aes(x=asinh(tpm), y=factor(gene_id, levels = rev(levels(factor(gene_id)))))) + geom_point(color="#A1D99B") + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(limits = c(-4,11)) + geom_errorbarh(aes(xmin=asinh(tpm-sd), xmax=asinh(tpm+sd)))
# ggplot(subset(geneorder2, geneorder2$V7 == "cluster 3"), aes(x=asinh(tpm), y=factor(gene_id, levels = rev(levels(factor(gene_id)))))) + geom_point(color="#005A32") + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(limits = c(-4,11)) + geom_errorbarh(aes(xmin=asinh(tpm-sd), xmax=asinh(tpm+sd)))
# ggplot(subset(geneorder2, geneorder2$V7 == "cluster 4"), aes(x=asinh(tpm), y=factor(gene_id, levels = rev(levels(factor(gene_id)))))) + geom_point(color="#4292C6") + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(limits = c(-4,11)) + geom_errorbarh(aes(xmin=asinh(tpm-sd), xmax=asinh(tpm+sd)))
# ggplot(subset(geneorder2, geneorder2$V7 == "cluster 5"), aes(x=asinh(tpm), y=factor(gene_id, levels = rev(levels(factor(gene_id)))))) + geom_point(color="#CB181D") + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_x_continuous(limits = c(-4,11)) + geom_errorbarh(aes(xmin=asinh(tpm-sd), xmax=asinh(tpm+sd)))


library(circlize)
library("rtracklayer")
# df <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4.4/circos/circos.Tak1v4.4.txt", sep = "\t", header = F)


# cenH3narrow <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4.3/homer/cenH3_Tak1v4.3.regions.narrow.bedgraph"))[,c(1:3,6)]
# cenH3broad <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4.3/homer/cenH3_Tak1v4.3.regions.broad.bedgraph"))[,c(1:3,6)]
# repeats <- as.data.frame(import.bedGraph("/volumes/berger/lab/Marchantia_v4.3/repeat_hits.100.shuffle.bed"))[,c(1:4)]
# H3broad <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4.3/homer/H3-1_Tak1v4.3.regions.broad.bedgraph"))[,c(1:3,6)]
# genes <- as.data.frame(import.bedGraph("/volumes/berger/lab/Marchantia_v4.3/gene_models/transfered_genes_merged.gff.bed"))[,c(1:4)]

# # H3K27me3 <- as.data.frame(import.bedGraph("/Volumes/berger/user/michael.borg/ChIP-seq/ULIChIP-seq/SNchip_all_v2/merged/merged_bigwigs/H3K27me3.100kb.log2r.bedgraph"))[,c(1:3,6)]
# # H3K4me3 <- as.data.frame(import.bedGraph("/Volumes/berger/user/michael.borg/ChIP-seq/ULIChIP-seq/SNchip_all_v2/merged/merged_bigwigs/H3K4me3.100kb.log2r.bedgraph"))[,c(1:3,6)]
# # HTR10 <- as.data.frame(import.bedGraph("/Volumes/berger/user/michael.borg/ChIP-seq/ULIChIP-seq/SNchip_all_v2/merged/merged_bigwigs/HTR10.100kb.log2r.bedgraph"))[,c(1:3,6)]

# bed.tracks <- list(genes, cenH3broad, cenH3narrow, H3broad, repeats)#, H3K27ac, H3K4me3, HTR10)
# # rm(H3K9me2, H3K27me1, H3K27me3, H3K27ac, H3K4me3, HTR10)
# clr <- c("#252525", "#01665e", "#a50026", "#f46d43", "#1a9850")#, "#4575b4")
# names(clr) <- c("Genes", "cenH3-broad", "cenH3-narrow", "H3-broad", "Repeat")


# # library(gridBase)
# # plot.new()
# # circle_size = unit(1, "snpc") # snpc unit gives you a square region

# # pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
# #     just = c("left", "center")))
# # par(omi = gridOMI(), new = TRUE)


# circos.clear()
# circos.par(gap.degree = 3, start.degree = 90, track.height = 0.10)
# circos.genomicInitialize(df)

# for (i in 1:length(bed.tracks)){
# circos.genomicDensity(bed.tracks[[i]], 
# 	col = clr[i]
# 	)
# # circos.genomicLabels(clr, labels.column = 4, side = "inside")
# }
# legend(-1.2,1.2, c("Genes", "cenH3-broad", "cenH3-narrow", "H3-broad", "Repeat"), lty=c(1,1), lwd=c(2.5,2.5),col=c("#252525", "#01665e", "#a50026", "#f46d43", "#1a9850"), bty="n") 

# # lgd_points = Legend(at = c("label1", "label2", "label3", "label4", "label5"), type = "points", 
# #                     legend_gp = gpar(col = 1:5), title_position = "topleft", 
# #                     title = "Track1")
# # # discrete
# # # lgd_links = Legend(at = c(-2, -1, 0, 1, 2), col_fun = col_fun, 
# # #                    title_position = "topleft", title = "Links")

# # lgd_list_vertical = packLegend(lgd_points)


# # upViewport()
# # pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical), 
# #     height = grobHeight(lgd_list_vertical), just = c("left", "center")))
# # grid.draw(lgd_list_vertical)
# # upViewport()


# circos.genomicDensity(repeats, baseline=0)
# circos.genomicDensity(cenH3broad, baseline=0)
# circos.genomicDensity(cenH3narrow, baseline=0)
# circos.genomicDensity(H3broad, baseline=0)
# circos.genomicDensity(genes, baseline=0)

# circos.clear()

# circos.par(gap.degree = 3, start.degree = 90, track.height = 0.10)
# circos.genomicInitialize(df)
# for (i in 1:length(bed.tracks)){
# circos.genomicTrack(bed.tracks[[i]], 
#                     panel.fun = function(region, value, ...) {
#                       circos.genomicLines(region, value, ytop.column = 1, ybottom = 0, col = clr [i], ...)
#                       circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 1, col = "#00000040")
#                     })
# }

df <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/circos/circos.Tak1v4.txt", sep = "\t", header = F)

H2A.Z <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/H2A.Z_Tak1v4.regions.narrow.bedgraph"))[,c(1:3,6)]
H3K9ac <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/H3K9ac_Tak1v4.regions.narrow.bedgraph"))[,c(1:3,6)]
H3K14ac <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/H3K14ac_Tak1v4.regions.narrow.bedgraph"))[,c(1:3,6)]
H3K4me1 <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/H3K4me1_Tak1v4.regions.narrow.bedgraph"))[,c(1:3,6)]
H3K4me3 <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/H3K4me3_Tak1v4.regions.narrow.bedgraph"))[,c(1:3,6)]
H3K36me3 <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/H3K36me3_Tak1v4.regions.narrow.bedgraph"))[,c(1:3,6)]
H3K9me1 <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/H3K9me1_Tak1v4.regions.narrow.bedgraph"))[,c(1:3,6)]
H3K27me1 <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/H3K27me1_Tak1v4.regions.narrow.bedgraph"))[,c(1:3,6)]
H3K27me3 <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/H3K27me3_Tak1v4.regions.narrow.bedgraph"))[,c(1:3,6)]
H3 <- as.data.frame(import.bedGraph("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/H3_Tak1v4.regions.narrow.bedgraph"))[,c(1:3,6)]
genes <- read.table("/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.gene.bed",header=F,col.names = c("seqnames","start","end","geneid","strand"))[,c(1:3)] 
genes$score <- 1
repeats <- read.table("/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.repeat.bed",header=F,col.names = c("seqnames","start","end","geneid","strand"))
ltr <- repeats[grep("LTR",repeats$geneid),c(1:3)]
ltr$score <- 1
colnames(ltr) <- c("seqnames","start","end","score")
line <- repeats[grep("LINE",repeats$geneid),c(1:3)]
line$score <- 1
colnames(line) <- c("seqnames","start","end","score")
helitron <- repeats[grep("RC/Helitron",repeats$geneid),c(1:3)]
helitron$score <- 1
colnames(helitron) <- c("seqnames","start","end","score")
# dna <- repeats[grep("DNA",repeats$geneid),c(1:3)]
# dna$score <- 1
# colnames(dna) <- c("seqnames","start","end","score")
# unknown <- repeats[grep("Unknown",repeats$geneid),c(1:3)]
# unknown$score <- 1
copia <- repeats[grep("LTR/Copia",repeats$geneid),c(1:3)]
copia$score <- 1
colnames(copia) <- c("seqnames","start","end","score")
gypsy <- repeats[grep("LTR/Gypsy",repeats$geneid),c(1:3)]
gypsy$score <- 1
colnames(gypsy) <- c("seqnames","start","end","score")
rtex <- repeats[grep("LINE/RTE-X",repeats$geneid),c(1:3)]
rtex$score <- 1
colnames(rtex) <- c("seqnames","start","end","score")
charlie <- repeats[grep("DNA/hAT-Charlie",repeats$geneid),c(1:3)]
charlie$score <- 1
colnames(charlie) <- c("seqnames","start","end","score")
repeats <- read.table("/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.repeat.bed",header=F,col.names = c("seqnames","start","end","geneid","strand"))
repeats$score <- 1
centromere <- read.table("/volumes/berger/user/sean.montgomery/Documents/Centromere/cenH3_repeat.e30.bed",header=F,col.names = c("seqnames","start","end","escore"))
centromere$score <- 1
gene1 <- subset(k4profile, k4profile$V6 == "cluster 1")
gene2 <- subset(k4profile, k4profile$V6 == "cluster 2")
gene3 <- subset(k4profile, k4profile$V6 == "cluster 3")
gene4 <- subset(k4profile, k4profile$V6 == "cluster 4")
gene5 <- subset(k4profile, k4profile$V6 == "cluster 5")
repeat1 <- subset(k5repeat, k5repeat$V6 == "cluster 1")
repeat2 <- subset(k5repeat, k5repeat$V6 == "cluster 2")
repeat3 <- subset(k5repeat, k5repeat$V6 == "cluster 3")
repeat4 <- subset(k5repeat, k5repeat$V6 == "cluster 4")
repeat5 <- subset(k5repeat, k5repeat$V6 == "cluster 5")
rtexcluster5 <- repeat5[grep("LINE/RTE-X",repeat5$V4),c(1:3)]

bed.tracks.genes <- list(genes, H3K9ac, H3K14ac, H3K4me1, H3K36me3, H3K4me3, H2A.Z, H3K27me3)
bed.tracks.repeats <- list(repeats, H3K9me1, H3K27me1, ltr, line, centromere)
# bed.tracks.repeats.types <- list(repeats, ltr, line, helitron, dna, unknown)
bed.tracks.cluster.gene <- list(gene1, gene2, gene3, gene4, gene5)
bed.tracks.cluster.repeat <- list(repeat1, repeat2, repeat3, repeat4, repeat5)
bed.tracks.repeats.2 <- list(repeats, H3K9me1, H3K27me1, copia, gypsy, helitron, charlie)
bed.tracks.repeats.centromere <- list(rtex, centromere, repeat5, rtexcluster5)

clr.genes <- c("#252525", "#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF")
names(clr.genes) <- c("Genes", "H3K9ac", "H3K14ac", "H3K4me1", "H3K36me3", "H3K4me3","H2A.Z", "H3K27me3")
clr.repeats <- c("#252525", "#CB181D","#FB6A4A", "#01665e", "#a50026","#000000")
names(clr.repeats) <- c("Repeats", "H3K9me1", "H3K27me1", "LTR", "LINE","Centromeric repeat")
clr.repeats.types <- c("#252525", "#01665e", "#a50026", "#f46d43", "#1a9850", "#4575b4")
names(clr.repeats.types) <- c("LTR/Copia", "LTR/Gypsy", "LINE/RTE-X","RC/Helitron","DNA/hAT-Charlie","Centromeric repeat")
clr.cluster.gene <- c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D")
names(clr.cluster.gene) <- c("Gene 1", "Gene 2", "Gene 3", "Gene 4", "Gene 5")
clr.cluster.repeat <- c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00")
names(clr.cluster.repeat) <- c("Repeat 1","Repeat 2","Repeat 3","Repeat 4","Repeat 5")
clr.repeats.2 <- c("#252525","#CB181D","#FB6A4A","#252525", "#01665e", "#f46d43", "#1a9850")
names(clr.repeats.2) <- c("Repeats","H3K9me1", "H3K27me1","LTR/Copia", "LTR/Gypsy","RC/Helitron","DNA/hAT-Charlie")
clr.repeats.centromere <- c("#a50026","#000000","#cdcd00","#D5A021")
names(clr.repeats.centromere) <- c("LINE/RTE-X","Centromeric repeat","Repeat 5","LINE/RTE-X Repeat 5")

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/circos.genes.pdf")
circos.clear()
circos.par(gap.degree = 3, start.degree = 90, track.height = 0.07)
circos.genomicInitialize(df)
for (i in 1:length(bed.tracks.genes)){
    circos.genomicDensity(bed.tracks.genes[[i]], 
                          col = clr.genes[i]
    )
}
legend(0.6,1.1, c("Genes", "H3K9ac", "H3K14ac", "H3K4me1", "H3K36me3", "H3K4me3","H2A.Z", "H3K27me3"), lty=c(1,1), lwd=c(2.5,2.5),col=c("#252525", "#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF"), bty="n") 
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/circos.repeats.pdf")
circos.clear()
circos.par(gap.degree = 3, start.degree = 90, track.height = 0.10)
circos.genomicInitialize(df)
for (i in 1:length(bed.tracks.repeats)){
circos.genomicDensity(bed.tracks.repeats[[i]], 
	col = clr.repeats[i]
	)
}
legend(0.6,1.1, c("Repeats", "H3K9me1", "H3K27me1", "LTR", "LINE","Centromeric repeat"), lty=c(1,1), lwd=c(2.5,2.5),col=c("#252525", "#CB181D","#FB6A4A", "#01665e", "#a50026","#000000"), bty="n") 
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/circos.cluster.gene.pdf")
circos.clear()
circos.par(gap.degree = 3, start.degree = 90, track.height = 0.10)
circos.genomicInitialize(df)
for (i in 1:length(bed.tracks.cluster.gene)){
circos.genomicDensity(bed.tracks.cluster.gene[[i]], 
	col = clr.cluster.gene[i]
	)
}
legend(0.6,1.1, c("Gene 1", "Gene 2", "Gene 3", "Gene 4", "Gene 5"), lty=c(1,1), lwd=c(2.5,2.5),col=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"), bty="n") 
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/circos.cluster.repeat.pdf")
circos.clear()
circos.par(gap.degree = 3, start.degree = 90, track.height = 0.10)
circos.genomicInitialize(df)
for (i in 1:length(bed.tracks.cluster.repeat)){
circos.genomicDensity(bed.tracks.cluster.repeat[[i]], 
	col = clr.cluster.repeat[i]
	)
}
legend(0.6,1.1, c("Repeat 1","Repeat 2","Repeat 3","Repeat 4","Repeat 5"), lty=c(1,1), lwd=c(2.5,2.5),col=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"), bty="n") 
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/circos.repeats.2.pdf")
circos.clear()
circos.par(gap.degree = 3, start.degree = 90, track.height = 0.10)
circos.genomicInitialize(df)
for (i in 1:length(bed.tracks.repeats.2)){
circos.genomicDensity(bed.tracks.repeats.2[[i]], 
	col = clr.repeats.2[i]
	)
}
legend(0.6,1.1, c("Repeats","H3K9me1", "H3K27me1","LTR/Copia", "LTR/Gypsy","RC/Helitron","DNA/hAT-Charlie"), lty=c(1,1), lwd=c(2.5,2.5),col=c("#252525","#CB181D","#FB6A4A","#252525", "#01665e", "#f46d43", "#1a9850"), bty="n") 
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/circos.repeats.centromere.pdf")
circos.clear()
circos.par(gap.degree = 3, start.degree = 90, track.height = 0.10)
circos.genomicInitialize(df)
for (i in 1:length(bed.tracks.repeats.centromere)){
circos.genomicDensity(bed.tracks.repeats.centromere[[i]], 
	col = clr.repeats.centromere[i]
	)
}
legend(0.6,1.1, c("LINE/RTE-X","Centromeric repeat","Repeat 5","LINE/RTE-X Repeat 5"), lty=c(1,1), lwd=c(2.5,2.5),col=c("#a50026","#000000","#cdcd00","#D5A021"), bty="n") 
dev.off()

# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/circos.repeats.types.pdf")
# circos.clear()
# circos.par(gap.degree = 3, start.degree = 90, track.height = 0.10)
# circos.genomicInitialize(df)
# for (i in 1:length(bed.tracks.repeats.types)){
#     circos.genomicDensity(bed.tracks.repeats.types[[i]], 
#                           col = clr.repeats.types[i]
#     )
# }
# legend(0.6,1.1, c("Repeats", "LTR", "LINE","Helitron","DNA","Unknown"), lty=c(1,1), lwd=c(2.5,2.5),col=c("#252525", "#01665e", "#a50026", "#f46d43", "#1a9850", "#4575b4"), bty="n")
# dev.off()

#!/usr/bin/env Rscript

# CIPHER DGE Script
# Input from CIPHER is several folders under a single directory located
# in results/rnaseq/salmon ... folders are named using ID from CIPHER
# ARGUMENTS: 1 -> main working directory
# 2 -> file of sample names and conditions
# 3 -> First condition
# 4 -> Second condition

# Import libraries
library()

# Set up variable to control command line arguments
args <- commandArgs(TRUE)

# Set working directory
# work_dir <- args[1]

## Using SNPsplit data
NAME <- args[1]

#Set/input sample names?
#Set/input window size?

#Loop over samples?
scaffold <- "chr1"
# sample <- read.table(NAME)
# #Loop over chromosomes?
# samplespline <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
# #Better way to plot blank?
# plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", ylim=c(0,2000))
# lines(samplespline, col="green")

pdf(paste("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/",scaffold,".H3K36me3.pdf",sep=""))
sample <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/bed/H3K36me3_Tak1v4.sizednuc150.100kb.bed")
plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", ylim=c(0,3000))
samplespline <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
lines(samplespline, col="#A1D99B")
dev.off()
pdf(paste("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/",scaffold,".H3K27me1.pdf",sep=""))
plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", ylim=c(0,3000))
sample <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/bed/H3K27me1_Tak1v4.sizednuc150.100kb.bed")
samplespline <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
lines(samplespline, col="#FB6A4A")
dev.off()
pdf(paste("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/",scaffold,".H3K27me3.pdf",sep=""))
plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", ylim=c(0,3000))
sample <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/bed/H3K27me3_Tak1v4.sizednuc150.100kb.bed")
samplespline <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
lines(samplespline, col="#C6DBEF")
dev.off()
pdf(paste("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/",scaffold,".H3K9me1.pdf",sep=""))
plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", ylim=c(0,3000))
sample <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/bed/H3K9me1_Tak1v4.sizednuc150.100kb.bed")
samplespline <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
lines(samplespline, col="#CB181D")
dev.off()
pdf(paste("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/",scaffold,".H3K4me1.pdf",sep=""))
plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", ylim=c(0,3000))
sample <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/bed/H3K4me1_Tak1v4.sizednuc150.100kb.bed")
samplespline <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
lines(samplespline, col="#74C476")
dev.off()
pdf(paste("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/",scaffold,".H3K4me3.pdf",sep=""))
plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", ylim=c(0,3000))
sample <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/bed/H3K4me3_Tak1v4.sizednuc150.100kb.bed")
samplespline <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
lines(samplespline, col="#084594")
dev.off()
pdf(paste("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/",scaffold,".H3K9ac.pdf",sep=""))
plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", ylim=c(0,3000))
sample <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/bed/H3K9ac_Tak1v4.sizednuc150.100kb.bed")
samplespline <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
lines(samplespline, col="#005A32")
dev.off()
pdf(paste("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/",scaffold,".H3K14ac.pdf",sep=""))
plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", ylim=c(0,3000))
sample <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/bed/H3K14ac_Tak1v4.sizednuc150.100kb.bed")
samplespline <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
lines(samplespline, col="#238B45")
dev.off()
pdf(paste("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/",scaffold,".H2A.Z.pdf",sep=""))
plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", ylim=c(0,3000))
sample <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/bed/H2A.Z_Tak1v4.sizednuc150.100kb.bed")
samplespline <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
lines(samplespline, col="#4292C6")
dev.off()



# "#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"






#Repeats
scalefactor <- 20
repeats <- read.table("/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.repeat.bed")
repeatschr1 <- subset(repeats, repeats$V1 == "chr1")[grep("LTR",repeats$V4),]
repeatschr1bin <- as.data.frame(table(cut(as.numeric(as.character(unlist(repeatschr1[2]))), seq(from=0, to =30825580, by = 100000), include.lowest = TRUE)))
repeatschr1bin[1] <- seq(from=100000, to =30825580, by = 100000)
repeatschr1bin[2] <- repeatschr1bin[2]*scalefactor
repeatschr1sp <- smooth.spline(repeatschr1bin)
lines(repeatschr1sp, col="black")


# repeatsscaff1 <- subset(repeats, repeats[1]==scaffold)
# repeatsscaff1bin <- as.data.frame(table(cut(as.numeric(as.character(unlist(repeatsscaff1[2]))), seq(from=0, to =30825580, by = 10000), include.lowest = TRUE)))
# repeatsscaff1bin[1] <- seq(from=10000, to =30825580, by = 10000)
# repeatsscaff1bin[2] <- repeatsscaff1bin[2]*scalefactor
# repeatsscaff1sp <- smooth.spline(repeatsscaff1bin)
# lines(repeatsscaff1sp, col="black")
# #Subset repeat type
# repeattype <- "LINE/RTE-X"
# repeatsscaff1sub <- subset(repeatsscaff1, repeatsscaff1[4]==repeattype)
# repeatsscaff1bin <- as.data.frame(table(cut(as.numeric(as.character(unlist(repeatsscaff1sub[2]))), seq(from=0, to =30825580, by = 10000), include.lowest = TRUE)))
# repeatsscaff1bin[1] <- seq(from=10000, to =30825580, by = 10000)
# repeatsscaff1bin[2] <- repeatsscaff1bin[2]*scalefactor
# repeatsscaff1sp <- smooth.spline(repeatsscaff1bin)
# lines(repeatsscaff1sp, col="cyan")

# chr1	30.584173
# chr2	29.643427
# chr3	27.142341
# chr4	26.988051
# chr5	26.794015
# chr6	23.861560
# chr7	21.963529
# chr8	21.314552
# chrY	7.543715

# scaffold <- "chr2"
# sample <- bhagyk27me3
# sample[2] <- sample[2]/1000000
# sample[3] <- sample[3]/1000000
# samplespline7 <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
# sample <- bhagyk27me1
# sample[2] <- sample[2]/1000000
# sample[3] <- sample[3]/1000000
# samplespline8 <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
# sample <- bhagyk36me3
# sample[2] <- sample[2]/1000000
# sample[3] <- sample[3]/1000000
# samplespline9 <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
# plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", xlab="Mbp", ylab="Reads per 10kb", main="Chr1", ylim=c(0,20000))
# lines(samplespline9, col="green")
# lines(samplespline7, col="blue")
# lines(samplespline8, col="orange")

# pdf("chr6.pdf")
# scaffold <- "chr6"
# sample <- h3k27me3.1.100
# sample[2] <- sample[2]/1000000
# sample[3] <- sample[3]/1000000
# samplespline7 <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
# sample <- h3k27me1.1.100
# sample[2] <- sample[2]/1000000
# sample[3] <- sample[3]/1000000
# samplespline8 <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
# sample <- h3k36me3.2.100
# sample[2] <- sample[2]/1000000
# sample[3] <- sample[3]/1000000
# samplespline9 <- smooth.spline(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]))
# plot(data.frame(x=subset(sample, sample[1]==scaffold)[2], y=subset(sample, sample[1]==scaffold)[4]), col="white", xlab="Mbp", ylab="Reads per 100kb", main="Chr6", ylim=c(0,12000),xlim=c(1,22.761782))
# lines(samplespline9, col="green")
# lines(samplespline7, col="blue")
# lines(samplespline8, col="orange")
# legend(x=15,y=12000,title="Histone features",legend=c("H3K27me3","H3K36me3","H3K27me1"),col=c("blue","green","orange"),lty=1)
# dev.off()

# # pdf("test2.pdf")
# # plot(c(0,33.138655),ylim=c(0,2000),col="white", ylab="Reads per 10kb",xlab="Mbp")
# # legend(x=22,y=4,title="Histone features",legend=c("H3K27me3","H3K36me3","H3K27me1"),col=c("blue","green","orange"),lty=1)
# # lines(samplespline9, col="green")
# # lines(samplespline7, col="blue")
# # lines(samplespline8, col="orange")
# # dev.off()

# scalefactor <- 20
# repeatsscaff1 <- subset(repeats, repeats[1]==scaffold)
# repeatsscaff1bin <- as.data.frame(table(cut(as.numeric(as.character(unlist(repeatsscaff1[2]))), seq(from=0, to =30825580, by = 10000), include.lowest = TRUE)))
# repeatsscaff1bin[1] <- seq(from=10000, to =30825580, by = 10000)
# repeatsscaff1bin[2] <- repeatsscaff1bin[2]*scalefactor
# repeatsscaff1sp <- smooth.spline(repeatsscaff1bin)
# lines(repeatsscaff1sp, col="black")
# repeattype <- "LINE/RTE-X"
# repeatsscaff1sub <- subset(repeatsscaff1, repeatsscaff1[4]==repeattype)
# repeatsscaff1bin <- as.data.frame(table(cut(as.numeric(as.character(unlist(repeatsscaff1sub[2]))), seq(from=0, to =30825580, by = 10000), include.lowest = TRUE)))
# repeatsscaff1bin[1] <- seq(from=10000, to =30825580, by = 10000)
# repeatsscaff1bin[2] <- repeatsscaff1bin[2]*scalefactor
# repeatsscaff1sp <- smooth.spline(repeatsscaff1bin)
# lines(repeatsscaff1sp, col="cyan")



breaks.corn <- seq(from = 0, to = 301476924, by = 100000)
bintable.gene.corn <- data.frame(breaks.corn)
bintable.repeat.corn <- data.frame(breaks.corn)
repeats.corn <- read.table("/volumes/berger/lab/cluster_files/sean/ZmaysE18/ZmaysE18.repeat.bed",header=F,col.names = c("seqnames","start","end","geneid","strand"))[,c(1:3)]
repeats.corn$bp <- (repeats.corn$end - repeats.corn$start)
genes.corn <- read.table("/volumes/berger/lab/cluster_files/sean/ZmaysE18/ZmaysE18.gene.bed",header=F,col.names = c("seqnames","start","end","geneid","strand"))[,c(1:3)]
genes.corn$bp <- (genes.corn$end - genes.corn$start)
bintable <- bintable.repeat.corn; feature <- repeats.corn
bintable <- bintable.gene.corn; feature <- genes.corn
breaks <- breaks.corn

breaks.mp <- seq(from = 0, to = 30584173, by = 100000)
bintable.gene.mp <- data.frame(breaks.mp)
bintable.repeat.mp <- data.frame(breaks.mp)
repeats.mp <- read.table("/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.repeat.bed",header=F,sep = '\t',col.names = c("seqnames","start","end","geneid","strand"))[,c(1:3)]
repeats.mp$bp <- (repeats.mp$end - repeats.mp$start)
genes.mp <- read.table("/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.gene.bed",header=F,col.names = c("seqnames","start","end","geneid","strand"))[,c(1:3)]
genes.mp$bp <- (genes.mp$end - genes.mp$start)
bintable <- bintable.repeat.mp; feature <- repeats.mp
bintable <- bintable.gene.mp; feature <- genes.mp
breaks <- breaks.mp

breaks.hg <- seq(from = 0, to = 249250621, by = 100000)
bintable.gene.hg <- data.frame(breaks.hg)
bintable.repeat.hg <- data.frame(breaks.hg)
repeats.hg <- read.table("/volumes/berger/lab/cluster_files/sean/hg19/hg19.repeat.bed",header=F,sep = '\t',col.names = c("seqnames","start","end","geneid","strand"))[,c(1:3)]
repeats.hg$bp <- (repeats.hg$end - repeats.hg$start)
genes.hg <- read.table("/volumes/berger/lab/cluster_files/sean/hg19/hg19.gene.bed",header=F,col.names = c("seqnames","start","end","geneid","strand"))[,c(1:3)]
genes.hg$bp <- (genes.hg$end - genes.hg$start)

bintable <- bintable.repeat.hg; feature <- repeats.hg
bintable <- bintable.gene.hg; feature <- genes.hg

i <- 1
j <- 1
u <- 0
summer <- 0
extrasum <- 0
while(u<=bintable[nrow(bintable),1]){
	if (feature[i,2]<u){
		if (feature[i,3]>u){
			summer <- (summer + (u - feature[i,2]))
			extrasum <- (extrasum + (feature[i,3] - u))
		}
		else {
			summer <- (summer + feature[i,4])
		}
		i <- (i+1)
	}
	else {
		bintable[j,2] <- summer
		summer <- extrasum
		extrasum <- 0
		j <- (j+1)
		u <- bintable[j,1]
	}
}

ggplot(data=bintable, aes(x=breaks.hg,y=V2))+geom_bar(stat="identity")+ylim(c(0,100000))

plot(data.frame(data.frame(x=bintable[1], y=bintable[2])), col="white", ylim=c(0,100000))
samplespline <- smooth.spline(data.frame(x=bintable[1], y=bintable[2]))
lines(samplespline, col="#74C476")

samplespline <- smooth.spline(data.frame(x=bintable[1], y=bintable[2]))
lines(samplespline, col="#CB181D")
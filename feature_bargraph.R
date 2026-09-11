
#Manual creation of bar plots of feature overlaps with peaks from HOMER
library(RColorBrewer)
sequential <- brewer.pal(6, "BuGn")

cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/feature_overlap_summary.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.bysample.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(9,3,8,2,6,5,7,1,4,10),]), 
        col=c("#005100","#006a00", "#008400", "#009d00", "#E69F00","#00d000") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,90))
mtext(side = 1, "Feature", line = 5)
legend(x = c(60, 90), y = c(0.6, 1), #adjust position of legend
       legend = colnames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.byfeature.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov[c(9,3,8,2,6,5,7,1,4,10),])), 
        col=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,70))
mtext(side = 1, "Feature", line = 5)
legend(x = c(40, 70), y = c(0.4, 1), #adjust position of legend
       legend = rownames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/sample_overlap_summary.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.bysample.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(9,3,8,2,6,5,7,1,4,10),]), 
        col=c("#005100","#006a00", "#008400", "#009d00", "#E69F00","#00d000") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,90))
mtext(side = 1, "Feature", line = 5)
legend(x = c(60, 90), y = c(0.6, 1), #adjust position of legend
       legend = colnames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.byfeature.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov[c(9,3,8,2,6,5,7,1,4,10),])), 
        col=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,70))
mtext(side = 1, "Feature", line = 5)
legend(x = c(40, 70), y = c(0.4, 1), #adjust position of legend
       legend = rownames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

## By gene clusters
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/feature_overlap_summary.clusters.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.bysample.clusters.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(9,3,8,2,6,5,7,1,4,10),]), 
        col=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,90))
mtext(side = 1, "Feature", line = 5)
legend(x = c(60, 90), y = c(0.6, 1), #adjust position of legend
       legend = colnames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.byfeature.clusters.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov[c(9,3,8,2,6,5,7,1,4,10),])), 
        col=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,70))
mtext(side = 1, "Feature", line = 5)
legend(x = c(40, 70), y = c(0.4, 1), #adjust position of legend
       legend = rownames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/sample_overlap_summary.clusters.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.bysample.clusters.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(9,3,8,2,6,5,7,1,4,10),]), 
        col=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,90))
mtext(side = 1, "Feature", line = 5)
legend(x = c(60, 90), y = c(0.6, 1), #adjust position of legend
       legend = colnames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.byfeature.clusters.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov[c(9,3,8,2,6,5,7,1,4,10),])), 
        col=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,70))
mtext(side = 1, "Feature", line = 5)
legend(x = c(40, 70), y = c(0.4, 1), #adjust position of legend
       legend = rownames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()


##By expression quintiles
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/feature_overlap_summary.quintiles.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.bysample.quintiles.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(9,3,8,2,6,5,7,1,4,10),]), 
        col=c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,90))
mtext(side = 1, "Feature", line = 5)
legend(x = c(60, 90), y = c(0.6, 1), #adjust position of legend
       legend = colnames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.byfeature.quintiles.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov[c(9,3,8,2,6,5,7,1,4,10),])), 
        col=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,79))
mtext(side = 1, "Feature", line = 5)
legend(x = c(50, 79), y = c(0.4, 1), #adjust position of legend
       legend = rownames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/sample_overlap_summary.quintiles.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.bysample.quintiles.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(9,3,8,2,6,5,7,1,4,10),]), 
        col=c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,90))
mtext(side = 1, "Feature", line = 5)
legend(x = c(60, 90), y = c(0.6, 1), #adjust position of legend
       legend = colnames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.byfeature.quintiles.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov[c(9,3,8,2,6,5,7,1,4,10),])), 
        col=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,79))
mtext(side = 1, "Feature", line = 5)
legend(x = c(50, 79), y = c(0.4, 1), #adjust position of legend
       legend = rownames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()


##By quintiles for ATAC-seq data
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/feature_overlap_summary.atac.quintiles.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.bysample.atac.quintiles.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(9,3,8,2,6,5,7,1,4,10),]), 
        col=c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,90))
mtext(side = 1, "Feature", line = 5)
legend(x = c(60, 90), y = c(0.6, 1), #adjust position of legend
       legend = colnames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.byfeature.atac.quintiles.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov[c(9,3,8,2,6,5,7,1,4,10),])), 
        col=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,79))
mtext(side = 1, "Feature", line = 5)
legend(x = c(50, 79), y = c(0.4, 1), #adjust position of legend
       legend = rownames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/sample_overlap_summary.atac.quintiles.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.bysample.atac.quintiles.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(9,3,8,2,6,5,7,1,4,10),]), 
        col=c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,90))
mtext(side = 1, "Feature", line = 5)
legend(x = c(60, 90), y = c(0.6, 1), #adjust position of legend
       legend = colnames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.byfeature.atac.quintiles.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov[c(9,3,8,2,6,5,7,1,4,10),])), 
        col=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,79))
mtext(side = 1, "Feature", line = 5)
legend(x = c(50, 79), y = c(0.4, 1), #adjust position of legend
       legend = rownames(cov[c(9,3,8,2,6,5,7,1,4,10),]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

#ATAC (feature) by gene cluster (sample)
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/feature_overlap_summary.atac.quintiles.cluster.gene.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.bysample.atac.quintiles.cluster.gene.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov), 
        col=c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,30))
mtext(side = 1, "Cluster", line = 5)
legend(x = c(20, 30), y = c(0.7, 1), #adjust position of legend
       legend = colnames(cov),  
       fill = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B"), # reorders so legend order matches graph
       title = "Coverage of quintiles")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.byfeature.atac.quintiles.cluster.gene.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov)), 
        col=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,30))
mtext(side = 1, "Cluster", line = 5)
legend(x = c(20, 30), y = c(0.7, 1), #adjust position of legend
       legend = rownames(cov),  
       fill = c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/sample_overlap_summary.atac.quintiles.cluster.gene.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.bysample.atac.quintiles.cluster.gene.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov), 
        col=c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,30))
mtext(side = 1, "Cluster", line = 5)
legend(x = c(20, 30), y = c(0.7, 1), #adjust position of legend
       legend = colnames(cov),  
       fill = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.byfeature.atac.quintiles.cluster.gene.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov)), 
        col=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,30))
mtext(side = 1, "Cluster", line = 5)
legend(x = c(20, 30), y = c(0.7, 1), #adjust position of legend
       legend = rownames(cov),  
       fill = c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()


#ATAC (feature) by repeat cluster (sample)
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/feature_overlap_summary.atac.quintiles.cluster.repeat.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.bysample.atac.quintiles.cluster.repeat.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov), 
        col=c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,30))
mtext(side = 1, "Cluster", line = 5)
legend(x = c(20, 30), y = c(0.7, 1), #adjust position of legend
       legend = colnames(cov),  
       fill = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B"), # reorders so legend order matches graph
       title = "Coverage of quintiles")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/feature_summary.byfeature.atac.quintiles.cluster.repeat.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov)), 
        col=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,30))
mtext(side = 1, "Cluster", line = 5)
legend(x = c(20, 30), y = c(0.7, 1), #adjust position of legend
       legend = rownames(cov),  
       fill = c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/sample_overlap_summary.atac.quintiles.cluster.repeat.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.bysample.atac.quintiles.cluster.repeat.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov), 
        col=c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,30))
mtext(side = 1, "Cluster", line = 5)
legend(x = c(20, 30), y = c(0.7, 1), #adjust position of legend
       legend = colnames(cov),  
       fill = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.byfeature.atac.quintiles.cluster.repeat.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov)), 
        col=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,30))
mtext(side = 1, "Cluster", line = 5)
legend(x = c(20, 30), y = c(0.7, 1), #adjust position of legend
       legend = rownames(cov),  
       fill = c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()


##By profiled mark
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/sample_overlap_summary.peaks.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.byrow.peaks.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(9,3,8,2,6,5,7,1,4,10),c(9,3,8,2,6,5,7,1,4,10)]), 
        col=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,140))
mtext(side = 1, "Feature", line = 5)
legend(x = c(95, 140), y = c(0.4, 1), #adjust position of legend
       legend = colnames(cov[c(9,3,8,2,6,5,7,1,4,10),c(9,3,8,2,6,5,7,1,4,10)]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.bycolumn.peaks.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov[c(9,3,8,2,6,5,7,1,4,10),c(9,3,8,2,6,5,7,1,4,10)])), 
        col=c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,140))
mtext(side = 1, "Feature", line = 5)
legend(x = c(95, 140), y = c(0.4, 1), #adjust position of legend
       legend = colnames(cov[c(9,3,8,2,6,5,7,1,4,10),c(9,3,8,2,6,5,7,1,4,10)]),  
       fill = c("#CB181D","#FB6A4A","#005A32","#238B45","#74C476","#A1D99B","#084594","#4292C6","#C6DBEF","yellow3"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()

##Distance of nearest gene in specified cluster to repeat in specified cluster
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.txt", header = T, row.names=1)
cov <- cov/1000
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/repeat_to_gene_distance.clusters.byrow.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov), 
        col=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Distance (kbp)", 
        font.lab=2,
        ylim=c(0,500),xlim=c(0,40))
mtext(side = 1, "Repeat cluster", line = 5)
legend(x = c(30, 40), y = c(300, 499), #adjust position of legend
       legend = colnames(cov),  
       fill = c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"), # reorders so legend order matches graph
       title = "Gene cluster")
dev.off()
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/repeat_to_gene_distance.clusters.bycolumn.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov)), 
        col=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Distance (kbp)", 
        font.lab=2,
        ylim=c(0,500),xlim=c(0,40))
mtext(side = 1, "Gene cluster", line = 5)
legend(x = c(30, 40), y = c(300, 499), #adjust position of legend
       legend = colnames(t(cov)),  
       fill = c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"), # reorders so legend order matches graph
       title = "Repeat cluster")
dev.off()
##Distance of gene in specified cluster to its nearest repeat in specified cluster
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/gene_to_repeat_distance.clusters.txt", header = T, row.names=1)
cov <- cov/1000
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/gene_to_repeat_distance.clusters.byrow.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov), 
        col=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Distance (kbp)", 
        font.lab=2,
        ylim=c(0,500),xlim=c(0,40))
mtext(side = 1, "Repeat cluster", line = 5)
legend(x = c(30, 40), y = c(300, 499), #adjust position of legend
       legend = colnames(cov),  
       fill = c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"), # reorders so legend order matches graph
       title = "Gene cluster")
dev.off()
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/gene_to_repeat_distance.clusters.bycolumn.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(t(cov)), 
        col=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Distance (kbp)", 
        font.lab=2,
        ylim=c(0,500),xlim=c(0,40))
mtext(side = 1, "Gene cluster", line = 5)
legend(x = c(30, 40), y = c(300, 499), #adjust position of legend
       legend = colnames(t(cov)),  
       fill = c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"), # reorders so legend order matches graph
       title = "Repeat cluster")
dev.off()

##Distances of repeats to their nearest gene neighbours of specific clusters
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
cov <- cov/1000
r1.g1 <- c(cov[,1])
r1.g2 <- c(cov[,2])
r1.g3 <- c(cov[,3])
r1.g4 <- c(cov[,4])
r1.g5 <- c(cov[,5])
r2.g1 <- c(cov[,6])
r2.g2 <- c(cov[,7])
r2.g3 <- c(cov[,8])
r2.g4 <- c(cov[,9])
r2.g5 <- c(cov[,10])
r3.g1 <- c(cov[,11])
r3.g2 <- c(cov[,12])
r3.g3 <- c(cov[,13])
r3.g4 <- c(cov[,14])
r3.g5 <- c(cov[,15])
r4.g1 <- c(cov[,16])
r4.g2 <- c(cov[,17])
r4.g3 <- c(cov[,18])
r4.g4 <- c(cov[,19])
r4.g5 <- c(cov[,20])
r5.g1 <- c(cov[,21])
r5.g2 <- c(cov[,22])
r5.g3 <- c(cov[,23])
r5.g4 <- c(cov[,24])
r5.g5 <- c(cov[,25])
df <- data.frame(
    x = c(r1.g1,r1.g2,r1.g3,r1.g4,r1.g5,r2.g1,r2.g2,r2.g3,r2.g4,r2.g5,r3.g1,r3.g2,r3.g3,r3.g4,r3.g5,r4.g1,r4.g2,r4.g3,r4.g4,r4.g5,r5.g1,r5.g2,r5.g3,r5.g4,r5.g5),
    y = rep(c("Repeat1","Repeat2","Repeat3","Repeat4","Repeat5"),each=221695),
    z = rep(rep(c("Gene1","Gene2","Gene3","Gene4","Gene5"),each=44339),5),
    stringsAsFactors = FALSE
)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/repeat_to_gene_distance.clusters.all.byrepeat.pdf",height=6.3,width=10)
ggplot(df, aes(y, x, fill=factor(z))) +
    geom_boxplot(outlier.colour=NA) +
    scale_fill_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"))+
    scale_x_discrete(name="Repeat cluster")+
    scale_y_continuous(name="Distance (kbp)") +
    labs(fill = "Gene cluster") +
    coord_cartesian(ylim = c(0, 1000))
dev.off()
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/repeat_to_gene_distance.clusters.all.bygene.pdf",height=6.3,width=10)
ggplot(df, aes(z, x, fill=factor(y))) +
    geom_boxplot(outlier.colour=NA) +
    scale_fill_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"))+
    scale_x_discrete(name="Gene cluster")+
    scale_y_continuous(name="Distance (kbp)") +
    labs(fill = "Repeat cluster") +
    coord_cartesian(ylim = c(0, 1000))
dev.off()
##Distance of genes to the nearest repeats of specified clusters
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/gene_to_repeat_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
cov <- cov/1000
g1.r1 <- c(cov[,1])
g1.r2 <- c(cov[,2])
g1.r3 <- c(cov[,3])
g1.r4 <- c(cov[,4])
g1.r5 <- c(cov[,5])
g2.r1 <- c(cov[,6])
g2.r2 <- c(cov[,7])
g2.r3 <- c(cov[,8])
g2.r4 <- c(cov[,9])
g2.r5 <- c(cov[,10])
g3.r1 <- c(cov[,11])
g3.r2 <- c(cov[,12])
g3.r3 <- c(cov[,13])
g3.r4 <- c(cov[,14])
g3.r5 <- c(cov[,15])
g4.r1 <- c(cov[,16])
g4.r2 <- c(cov[,17])
g4.r3 <- c(cov[,18])
g4.r4 <- c(cov[,19])
g4.r5 <- c(cov[,20])
g5.r1 <- c(cov[,21])
g5.r2 <- c(cov[,22])
g5.r3 <- c(cov[,23])
g5.r4 <- c(cov[,24])
g5.r5 <- c(cov[,25])
df <- data.frame(
    x = c(g1.r1,g1.r2,g1.r3,g1.r4,g1.r5,g2.r1,g2.r2,g2.r3,g2.r4,g2.r5,g3.r1,g3.r2,g3.r3,g3.r4,g3.r5,g4.r1,g4.r2,g4.r3,g4.r4,g4.r5,g5.r1,g5.r2,g5.r3,g5.r4,g5.r5),
    y = rep(c("Gene1","Gene2","Gene3","Gene4","Gene5"),each=40270),
    z = rep(rep(c("Repeat1","Repeat2","Repeat3","Repeat4","Repeat5"),each=8054),5),
    stringsAsFactors = FALSE
)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/gene_to_repeat_distance.clusters.all.byrepeat.pdf",height=6.3,width=10)
ggplot(df, aes(z, x, fill=factor(y))) +
    geom_boxplot(outlier.colour=NA) +
    scale_fill_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"))+
    scale_x_discrete(name="Repeat cluster")+
    scale_y_continuous(name="Distance (kbp)") +
    labs(fill = "Gene cluster") +
    coord_cartesian(ylim = c(0, 1000))
dev.off()
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/gene_to_repeat_distance.clusters.all.bygene.pdf",height=6.3,width=10)
ggplot(df, aes(y, x, fill=factor(z))) +
    geom_boxplot(outlier.colour=NA) +
    scale_fill_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"))+
    scale_x_discrete(name="Gene cluster")+
    scale_y_continuous(name="Distance (kbp)") +
    labs(fill = "Repeat cluster") +
    coord_cartesian(ylim = c(0, 1000))
dev.off()

# ##Distances of Genji elements to their nearest gene neighbours of specific clusters
# cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/genji/Genji_to_gene_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
# cov <- cov/1000
# r1.g1 <- c(cov[,1])
# r1.g2 <- c(cov[,2])
# r1.g3 <- c(cov[,3])
# r1.g4 <- c(cov[,4])
# r1.g5 <- c(cov[,5])
# r2.g1 <- c(cov[,6])
# r2.g2 <- c(cov[,7])
# r2.g3 <- c(cov[,8])
# r2.g4 <- c(cov[,9])
# r2.g5 <- c(cov[,10])
# r3.g1 <- c(cov[,11])
# r3.g2 <- c(cov[,12])
# r3.g3 <- c(cov[,13])
# r3.g4 <- c(cov[,14])
# r3.g5 <- c(cov[,15])
# df <- data.frame(
#     x = c(r1.g1,r1.g2,r1.g3,r1.g4,r1.g5,r2.g1,r2.g2,r2.g3,r2.g4,r2.g5,r3.g1,r3.g2,r3.g3,r3.g4,r3.g5),
#     y = rep(c("Genji1","Genji2","Genji3"),each=280),
#     z = rep(rep(c("Gene1","Gene2","Gene3","Gene4","Gene5"),each=56),3),
#     stringsAsFactors = FALSE
# )
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/genji/Genji_to_gene_distance.clusters.all.byGenji.pdf",height=6.3,width=10)
# ggplot(df, aes(y, x, fill=factor(z))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"))+
#     scale_x_discrete(name="Genji cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Gene cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/genji/Genji_to_gene_distance.clusters.all.bygene.pdf",height=6.3,width=10)
# ggplot(df, aes(z, x, fill=factor(y))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"))+
#     scale_x_discrete(name="Gene cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Genji cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()

# ##Distance of Genji elements to the nearest repeats of specified clusters
# cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/genji/Genji_to_repeat_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
# cov <- cov/1000
# g1.r1 <- c(cov[,1])
# g1.r2 <- c(cov[,2])
# g1.r3 <- c(cov[,3])
# g1.r4 <- c(cov[,4])
# g1.r5 <- c(cov[,5])
# g2.r1 <- c(cov[,6])
# g2.r2 <- c(cov[,7])
# g2.r3 <- c(cov[,8])
# g2.r4 <- c(cov[,9])
# g2.r5 <- c(cov[,10])
# g3.r1 <- c(cov[,11])
# g3.r2 <- c(cov[,12])
# g3.r3 <- c(cov[,13])
# g3.r4 <- c(cov[,14])
# g3.r5 <- c(cov[,15])
# df <- data.frame(
#     x = c(g1.r1,g1.r2,g1.r3,g1.r4,g1.r5,g2.r1,g2.r2,g2.r3,g2.r4,g2.r5,g3.r1,g3.r2,g3.r3,g3.r4,g3.r5),
#     y = rep(c("Genji1","Genji2","Genji3"),each=280),
#     z = rep(rep(c("Repeat1","Repeat2","Repeat3","Repeat4","Repeat5"),each=56),3),
#     stringsAsFactors = FALSE
# )
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/genji/Genji_to_repeat_distance.clusters.all.byrepeat.pdf",height=6.3,width=10)
# ggplot(df, aes(z, x, fill=factor(y))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"))+
#     scale_x_discrete(name="Repeat cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Genji cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/genji/Genji_to_repeat_distance.clusters.all.byGenji.pdf",height=6.3,width=10)
# ggplot(df, aes(y, x, fill=factor(z))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"))+
#     scale_x_discrete(name="Genji cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Repeat cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()

# ##Distances of HomeoGypsy elements to their nearest gene neighbours of specific clusters
# cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/HomeoGypsy_to_gene_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
# cov <- cov/1000
# r1.g1 <- c(cov[,1])
# r1.g2 <- c(cov[,2])
# r1.g3 <- c(cov[,3])
# r1.g4 <- c(cov[,4])
# r1.g5 <- c(cov[,5])
# r2.g1 <- c(cov[,6])
# r2.g2 <- c(cov[,7])
# r2.g3 <- c(cov[,8])
# r2.g4 <- c(cov[,9])
# r2.g5 <- c(cov[,10])
# r3.g1 <- c(cov[,11])
# r3.g2 <- c(cov[,12])
# r3.g3 <- c(cov[,13])
# r3.g4 <- c(cov[,14])
# r3.g5 <- c(cov[,15])
# df <- data.frame(
#     x = c(r1.g1,r1.g2,r1.g3,r1.g4,r1.g5,r2.g1,r2.g2,r2.g3,r2.g4,r2.g5,r3.g1,r3.g2,r3.g3,r3.g4,r3.g5),
#     y = rep(c("HomeoGypsy1","HomeoGypsy2","HomeoGypsy3"),each=2065), #Each is the number of observations *5
#     z = rep(rep(c("Gene1","Gene2","Gene3","Gene4","Gene5"),each=413),3), #Each is the number of observations
#     stringsAsFactors = FALSE
# )
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/HomeoGypsy_to_gene_distance.clusters.all.byHomeoGypsy.pdf",height=6.3,width=10)
# ggplot(df, aes(y, x, fill=factor(z))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"))+
#     scale_x_discrete(name="HomeoGypsy cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Gene cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/HomeoGypsy_to_gene_distance.clusters.all.bygene.pdf",height=6.3,width=10)
# ggplot(df, aes(z, x, fill=factor(y))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#FB6A4A","#2171B5","#cdcd00"))+
#     scale_x_discrete(name="Gene cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "HomeoGypsy cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# ##Distance of HomeoGypsy elements to the nearest repeats of specified clusters
# cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/HomeoGypsy_to_repeat_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
# cov <- cov/1000
# g1.r1 <- c(cov[,1])
# g1.r2 <- c(cov[,2])
# g1.r3 <- c(cov[,3])
# g1.r4 <- c(cov[,4])
# g1.r5 <- c(cov[,5])
# g2.r1 <- c(cov[,6])
# g2.r2 <- c(cov[,7])
# g2.r3 <- c(cov[,8])
# g2.r4 <- c(cov[,9])
# g2.r5 <- c(cov[,10])
# g3.r1 <- c(cov[,11])
# g3.r2 <- c(cov[,12])
# g3.r3 <- c(cov[,13])
# g3.r4 <- c(cov[,14])
# g3.r5 <- c(cov[,15])
# df <- data.frame(
#     x = c(g1.r1,g1.r2,g1.r3,g1.r4,g1.r5,g2.r1,g2.r2,g2.r3,g2.r4,g2.r5,g3.r1,g3.r2,g3.r3,g3.r4,g3.r5),
#     y = rep(c("HomeoGypsy1","HomeoGypsy2","HomeoGypsy3"),each=2065),
#     z = rep(rep(c("Repeat1","Repeat2","Repeat3","Repeat4","Repeat5"),each=413),3),
#     stringsAsFactors = FALSE
# )
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/HomeoGypsy_to_repeat_distance.clusters.all.byrepeat.pdf",height=6.3,width=10)
# ggplot(df, aes(z, x, fill=factor(y))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#FB6A4A","#2171B5","#cdcd00"))+
#     scale_x_discrete(name="Repeat cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "HomeoGypsy cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/HomeoGypsy_to_repeat_distance.clusters.all.byHomeoGypsy.pdf",height=6.3,width=10)
# ggplot(df, aes(y, x, fill=factor(z))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"))+
#     scale_x_discrete(name="HomeoGypsy cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Repeat cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# ##Distances of Gypsy elements to their nearest gene neighbours of specific clusters
# cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Gypsy_to_gene_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
# cov <- cov/1000
# r1.g1 <- c(cov[,1])
# r1.g2 <- c(cov[,2])
# r1.g3 <- c(cov[,3])
# r1.g4 <- c(cov[,4])
# r1.g5 <- c(cov[,5])
# r2.g1 <- c(cov[,6])
# r2.g2 <- c(cov[,7])
# r2.g3 <- c(cov[,8])
# r2.g4 <- c(cov[,9])
# r2.g5 <- c(cov[,10])
# r3.g1 <- c(cov[,11])
# r3.g2 <- c(cov[,12])
# r3.g3 <- c(cov[,13])
# r3.g4 <- c(cov[,14])
# r3.g5 <- c(cov[,15])
# df <- data.frame(
#     x = c(r1.g1,r1.g2,r1.g3,r1.g4,r1.g5,r2.g1,r2.g2,r2.g3,r2.g4,r2.g5,r3.g1,r3.g2,r3.g3,r3.g4,r3.g5),
#     y = rep(c("Gypsy1","Gypsy2","Gypsy3"),each=9030), #Each is the number of observations *5
#     z = rep(rep(c("Gene1","Gene2","Gene3","Gene4","Gene5"),each=1806),3), #Each is the number of observations
#     stringsAsFactors = FALSE
# )
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Gypsy_to_gene_distance.clusters.all.byGypsy.pdf",height=6.3,width=10)
# ggplot(df, aes(y, x, fill=factor(z))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"))+
#     scale_x_discrete(name="Gypsy cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Gene cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Gypsy_to_gene_distance.clusters.all.bygene.pdf",height=6.3,width=10)
# ggplot(df, aes(z, x, fill=factor(y))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#FB6A4A","#2171B5","#cdcd00"))+
#     scale_x_discrete(name="Gene cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Gypsy cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# ##Distance of Gypsy elements to the nearest repeats of specified clusters
# cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Gypsy_to_repeat_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
# cov <- cov/1000
# g1.r1 <- c(cov[,1])
# g1.r2 <- c(cov[,2])
# g1.r3 <- c(cov[,3])
# g1.r4 <- c(cov[,4])
# g1.r5 <- c(cov[,5])
# g2.r1 <- c(cov[,6])
# g2.r2 <- c(cov[,7])
# g2.r3 <- c(cov[,8])
# g2.r4 <- c(cov[,9])
# g2.r5 <- c(cov[,10])
# g3.r1 <- c(cov[,11])
# g3.r2 <- c(cov[,12])
# g3.r3 <- c(cov[,13])
# g3.r4 <- c(cov[,14])
# g3.r5 <- c(cov[,15])
# df <- data.frame(
#     x = c(g1.r1,g1.r2,g1.r3,g1.r4,g1.r5,g2.r1,g2.r2,g2.r3,g2.r4,g2.r5,g3.r1,g3.r2,g3.r3,g3.r4,g3.r5),
#     y = rep(c("Gypsy1","Gypsy2","Gypsy3"),each=12810),
#     z = rep(rep(c("Repeat1","Repeat2","Repeat3","Repeat4","Repeat5"),each=2562),3),
#     stringsAsFactors = FALSE
# )
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Gypsy_to_repeat_distance.clusters.all.byrepeat.pdf",height=6.3,width=10)
# ggplot(df, aes(z, x, fill=factor(y))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#FB6A4A","#2171B5","#cdcd00"))+
#     scale_x_discrete(name="Repeat cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Gypsy cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Gypsy_to_repeat_distance.clusters.all.byGypsy.pdf",height=6.3,width=10)
# ggplot(df, aes(y, x, fill=factor(z))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"))+
#     scale_x_discrete(name="Gypsy cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Repeat cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# ##Combine Gypsy and HomeoGypsy plots and do a Wilcoxon test
# cov2 <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Gypsy_to_gene_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
# cov2 <- cov2/1000
# r1.g12 <- c(cov2[,1])
# r1.g22 <- c(cov2[,2])
# r1.g32 <- c(cov2[,3])
# r1.g42 <- c(cov2[,4])
# r1.g52 <- c(cov2[,5])
# r2.g12 <- c(cov2[,6])
# r2.g22 <- c(cov2[,7])
# r2.g32 <- c(cov2[,8])
# r2.g42 <- c(cov2[,9])
# r2.g52 <- c(cov2[,10])
# r3.g12 <- c(cov2[,11])
# r3.g22 <- c(cov2[,12])
# r3.g32 <- c(cov2[,13])
# r3.g42 <- c(cov2[,14])
# r3.g52 <- c(cov2[,15])
# df2 <- data.frame(
#   x = c(r1.g12,r1.g22,r1.g32,r1.g42,r1.g52,r2.g12,r2.g22,r2.g32,r2.g42,r2.g52,r3.g12,r3.g22,r3.g32,r3.g42,r3.g52),
#   y = rep(c("Gypsy1","Gypsy2","Gypsy3"),each=9030), #Each is the number of observations *5
#   z = rep(rep(c("Gene1","Gene2","Gene3","Gene4","Gene5"),each=1806),3), #Each is the number of observations
#   stringsAsFactors = FALSE
# )
# df3 <- rbind(df, df2)
# df3$y <- factor(df3$y,levels = c("Gypsy1","HomeoGypsy1","Gypsy2","HomeoGypsy2","Gypsy3","HomeoGypsy3"))
# ggplot(df3, aes(z, x, fill=factor(y))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#FB6A4A","#963f2c","#C6DBEF","#2171B5","#cdcd00","#8f8f00"))+
#     scale_x_discrete(name="Gene cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Gypsy cluster") +
#     coord_cartesian(ylim = c(0, 2000))
# library(ggpubr)
# print(compare_means(x ~ y, data = df3, group.by = "z"), n=75)
# ##Distances of Copia elements to their nearest gene neighbours of specific clusters
# cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Copia_to_gene_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
# cov <- cov/1000
# r1.g1 <- c(cov[,1])
# r1.g2 <- c(cov[,2])
# r1.g3 <- c(cov[,3])
# r1.g4 <- c(cov[,4])
# r1.g5 <- c(cov[,5])
# r2.g1 <- c(cov[,6])
# r2.g2 <- c(cov[,7])
# r2.g3 <- c(cov[,8])
# r2.g4 <- c(cov[,9])
# r2.g5 <- c(cov[,10])
# r3.g1 <- c(cov[,11])
# r3.g2 <- c(cov[,12])
# r3.g3 <- c(cov[,13])
# r3.g4 <- c(cov[,14])
# r3.g5 <- c(cov[,15])
# r4.g1 <- c(cov[,16])
# r4.g2 <- c(cov[,17])
# r4.g3 <- c(cov[,18])
# r4.g4 <- c(cov[,19])
# r4.g5 <- c(cov[,20])
# r5.g1 <- c(cov[,21])
# r5.g2 <- c(cov[,22])
# r5.g3 <- c(cov[,23])
# r5.g4 <- c(cov[,24])
# r5.g5 <- c(cov[,25])
# df <- data.frame(
#     x = c(r1.g1,r1.g2,r1.g3,r1.g4,r1.g5,r2.g1,r2.g2,r2.g3,r2.g4,r2.g5,r3.g1,r3.g2,r3.g3,r3.g4,r3.g5,r4.g1,r4.g2,r4.g3,r4.g4,r4.g5,r5.g1,r5.g2,r5.g3,r5.g4,r5.g5),
#     y = rep(c("Copia1","Copia2","Copia3","Copia4","Copia5"),each=4085), #Each is the number of observations *5
#     z = rep(rep(c("Gene1","Gene2","Gene3","Gene4","Gene5"),each=817),5), #Each is the number of observations
#     stringsAsFactors = FALSE
# )
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Copia_to_gene_distance.clusters.all.byCopia.pdf",height=6.3,width=10)
# ggplot(df, aes(y, x, fill=factor(z))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#C6DBEF","#A1D99B","#005A32","#4292C6","#CB181D"))+
#     scale_x_discrete(name="Copia cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Gene cluster") +
#     coord_cartesian(ylim = c(0, 2000))
# dev.off()
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Copia_to_gene_distance.clusters.all.bygene.pdf",height=6.3,width=10)
# ggplot(df, aes(z, x, fill=factor(y))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#2171B5","#FB6A4A","#74C476","#FB6A4A","#cdcd00"))+
#     scale_x_discrete(name="Gene cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Copia cluster") +
#     coord_cartesian(ylim = c(0, 2000))
# dev.off()
# ##Distance of Copia elements to the nearest repeats of specified clusters
# cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Copia_to_repeat_distance.clusters.all.txt", header = T, fill = T,sep = '\t')
# cov <- cov/1000
# g1.r1 <- c(cov[,1])
# g1.r2 <- c(cov[,2])
# g1.r3 <- c(cov[,3])
# g1.r4 <- c(cov[,4])
# g1.r5 <- c(cov[,5])
# g2.r1 <- c(cov[,6])
# g2.r2 <- c(cov[,7])
# g2.r3 <- c(cov[,8])
# g2.r4 <- c(cov[,9])
# g2.r5 <- c(cov[,10])
# g3.r1 <- c(cov[,11])
# g3.r2 <- c(cov[,12])
# g3.r3 <- c(cov[,13])
# g3.r4 <- c(cov[,14])
# g3.r5 <- c(cov[,15])
# g4.r1 <- c(cov[,16])
# g4.r2 <- c(cov[,17])
# g4.r3 <- c(cov[,18])
# g4.r4 <- c(cov[,19])
# g4.r5 <- c(cov[,20])
# g5.r1 <- c(cov[,21])
# g5.r2 <- c(cov[,22])
# g5.r3 <- c(cov[,23])
# g5.r4 <- c(cov[,24])
# g5.r5 <- c(cov[,25])
# df <- data.frame(
#     x = c(g1.r1,g1.r2,g1.r3,g1.r4,g1.r5,g2.r1,g2.r2,g2.r3,g2.r4,g2.r5,g3.r1,g3.r2,g3.r3,g3.r4,g3.r5,g4.r1,g4.r2,g4.r3,g4.r4,g4.r5,g5.r1,g5.r2,g5.r3,g5.r4,g5.r5),
#     y = rep(c("Copia1","Copia2","Copia3","Copia4","Copia5"),each=8100),
#     z = rep(rep(c("Repeat1","Repeat2","Repeat3","Repeat4","Repeat5"),each=1620),5),
#     stringsAsFactors = FALSE
# )
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Copia_to_repeat_distance.clusters.all.byrepeat.pdf",height=6.3,width=10)
# ggplot(df, aes(z, x, fill=factor(y))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#2171B5","#FB6A4A","#74C476","#FB6A4A","#cdcd00"))+
#     scale_x_discrete(name="Repeat cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Copia cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()
# pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/homeogypsy/Copia_to_repeat_distance.clusters.all.byCopia.pdf",height=6.3,width=10)
# ggplot(df, aes(y, x, fill=factor(z))) +
#     geom_boxplot(outlier.colour=NA) +
#     scale_fill_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"))+
#     scale_x_discrete(name="Copia cluster")+
#     scale_y_continuous(name="Distance (kbp)") +
#     labs(fill = "Repeat cluster") +
#     coord_cartesian(ylim = c(0, 1000))
# dev.off()

##Distance of repeat clusters to any gene
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.repeatclusteronly.txt", header = T, fill = T,sep = '\t')
cov <- cov/1000
r1.g1 <- c(cov[,1])
r1.g2 <- c(cov[,2])
r1.g3 <- c(cov[,3])
r1.g4 <- c(cov[,4])
r1.g5 <- c(cov[,5])
df <- data.frame(
  x = c(r1.g1,r1.g2,r1.g3,r1.g4,r1.g5),
  y = rep(c("Repeat1","Repeat2","Repeat3","Repeat4","Repeat5"),each=43613),
  z = c("Gene"),
  stringsAsFactors = FALSE
)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/repeat_to_gene_distance.clusters.repeatclusteronly.pdf",height=6.3,width=10)
ggplot(df, aes(y, x, fill=factor(y))) +
    geom_boxplot(outlier.colour=NA) +
    scale_fill_manual(values=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"))+
    scale_x_discrete(name="Repeat cluster")+
    scale_y_continuous(name="Distance (kbp)") +
    labs(fill = "Gene cluster") +
    theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    coord_cartesian(ylim = c(0, 50))
dev.off()

##Number of repeats of a class in a repeat cluster
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratioperclass.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/repeat_classes.clusters.ratioperclass.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov), 
        col=c("#000000","#01665e","#c97f7f", "#a50026","#4C3100") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Ratio (Repeat per class/Total repeats per class", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,40))
mtext(side = 1, "Repeat cluster", line = 5)
legend(x = c(30, 40), y = c(0.7, 1), #adjust position of legend
       legend = colnames(cov),  
       fill = c("#000000","#01665e","#c97f7f", "#a50026","#4C3100"), # reorders so legend order matches graph
       title = "Repeat cluster")
dev.off()

cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratiopercluster.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/repeat_classes.clusters.ratiopercluster.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov), 
        col=c("#000000","#01665e","#c97f7f", "#a50026","#4C3100") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Ratio (Repeat per class/Total repeats in cluster", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,40))
mtext(side = 1, "Repeat cluster", line = 5)
legend(x = c(30, 40), y = c(0.7, 1), #adjust position of legend
       legend = colnames(cov),  
       fill = c("#000000","#01665e","#c97f7f", "#a50026","#4C3100"), # reorders so legend order matches graph
       title = "Repeat cluster")
dev.off()

##Number of repeats of a superfamily in a repeat cluster
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratioperclass.all.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/repeat_classes.clusters.ratioperclass.all.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(27,28,25,30,14),]), 
        col=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Ratio (Repeat per class/Total repeats per class)", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,40))
mtext(side = 1, "TE superfamily", line = 5)
legend(x = c(30, 40), y = c(0.7, 1), #adjust position of legend
       legend = colnames(cov[c(27,28,25,30,14),]),  
       fill = c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"), # reorders so legend order matches graph
       title = "Repeat cluster")
dev.off()

cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratiopercluster.all.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/repeat_classes.clusters.ratiopercluster.all.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(cov[c(27,28,25,30,14),]), 
        col=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Ratio (Repeat per class/Total repeats in cluster)", 
        font.lab=2,
        ylim=c(0,1),xlim=c(0,40))
mtext(side = 1, "TE superfamily", line = 5)
legend(x = c(30, 40), y = c(0.7, 1), #adjust position of legend
       legend = colnames(cov[c(27,28,25,30,14),]),  
       fill = c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"), # reorders so legend order matches graph
       title = "Repeat cluster")
dev.off()
covratios <- cov/rowSums(cov)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/repeat_classes.clusters.ratiopercluster.scaledperclass.pdf",height=6.3,width=10)
par(mai=c(1.52,0.82,0.82,0.42))
barplot(t(covratios[c(27,28,25,30,14),]), 
        col=c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00") , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=NULL, 
        las=2,
        ylab = "Proportion (TE per repeat cluster)",
        font.lab=2,
        ylim=c(0,1),xlim=c(0,40))
mtext(side = 1, "TE superfamily", line = 5)
legend(x = c(30, 40), y = c(0.7, 1), #adjust position of legend
       legend = colnames(covratios[c(27,28,25,30,14),]),  
       fill = c("#6BAED6","#74C476","#FB6A4A","#2171B5","#cdcd00"), # reorders so legend order matches graph
       title = "Repeat cluster")
dev.off()

#Stacked bar plot - doesn't look so good or work well when >100%
cov <- read.table("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/sample_overlap_summary.2.txt", header = T, row.names=1)
pdf("/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/plots/sample_summary.bysample.2.pdf",height=6.3,width=10)
ppar(mai=c(1.52,0.82,0.82,0.42)) #size of margins to fit axis labels
barplot(t(cov[c(9,3,8,2,6,5,7,1,4,10),c(1,2,3,4,5,6,7,8)]),
        names.arg = cov$Sample,
        cex.names = 0.6, # size of x axis labels
        col = c("#238B45", "#41AB5D", "#74C476", "#A1D99B","#CB181D","#EF3B2C","#FB6A4A","#FC9272"), # colors
        border="white", 
        font.axis=2, 
        legend=NULL, 
        ylab = "Proportion of mark overlapping genomic feature",
        font.lab=2,
        las=2, #x axis labels orientation to graph
        ylim=c(0,1),) #y axis limits
mtext(side = 1, "Sample", line = 5)
legend(x = c(7.5, 12), y = c(0.6,1), #adjust position of legend
       legend = c("5'UTR","CDS","3'UTR","Intron","rRNA","Satellite","Unknown repeat","Transposon"), #in order from top to bottom
       fill = c("#238B45", "#41AB5D", "#74C476", "#A1D99B","#CB181D","#EF3B2C","#FB6A4A","#FC9272"), # reorders so legend order matches graph
       title = "Coverage of features")
dev.off()



# library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# files <- list.files(pattern = ".*broad.bedgraph")
# names(files) <- gsub("_Tak1v4.4.regions.broad.bedgraph", "", files)
# peak <- readPeakFile(filelist[[14]])
# txdb <- makeTxDbFromGFF("/volumes/berger/lab/Marchantia_v4.4/annotation_MpTak1r1/MpTak1.r1.nochrname.gff3")
# promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
# tagMatrix <- getTagMatrix(peak, windows=promoter)
# peakAnno <- annotatePeak(files[[4]], tssRegion=c(-1000, 1000), TxDb=txdb)
# #Can also compare other features
# pdf("peakAnno.pdf",height=6.3,width=10)
# for (f in files){
#     peakAnno <- annotatePeak(f, tssRegion=c(-1000, 1000), TxDb=txdb)
#     plotAnnoPie(peakAnno)
#     title(f, line = -10)
# }
# dev.off()
###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)
library("gplots")
library(edgeR)

getcoverGnum <- function(indata){
	return(length(which(indata > 0)))
}
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
fpkm_1 <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/WQN_DUXno_mmESCs_GeneSymbol_fpkm.txt",header=T,row.names=1,sep="\t")
all_samples_1 <- colnames(fpkm_1)
all_genes_1 <- row.names(fpkm_1)
logfpkm_1 <- log2(fpkm_1+1)

fpkm_2 <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/WQN_duxpu_mmESCs_GeneSymbol_fpkm.txt",header=T,row.names=1,sep="\t")
all_samples_2 <- colnames(fpkm_2)
all_genes_2 <- row.names(fpkm_2)
logfpkm_2 <- log2(fpkm_2+1)

all_genes <- intersect(all_genes_1,all_genes_2)
fpkm <- cbind(fpkm_1[all_genes,],fpkm_2[all_genes,])
all_samples <- colnames(fpkm)
logfpkm <- log2(fpkm+1)

GSE66582_fpkm <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/data/GSE66582_mmEmbryo_RNAseq/GSE66582_GeneSymbol_fpkm.txt",header=T,row.names=1,sep="\t")
GSE66582_samples <- colnames(GSE66582_fpkm)
GSE66582_genes <- row.names(GSE66582_fpkm)
GSE66582_logfpkm <- log2(GSE66582_fpkm+1)

GSE85632_fpkm <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/data/GSE85632_mmESC_RNAseq/GSE85632_GeneSymbol_fpkm.txt",header=T,row.names=1,sep="\t")
GSE85632_samples <- colnames(GSE85632_fpkm)
GSE85632_genes <- row.names(GSE85632_fpkm)
GSE85632_logfpkm <- log2(GSE85632_fpkm+1)

GSE94325_fpkm <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/data/GSE94325_mmESC_RNAseq/GSE94325_GeneSymbol_fpkm.txt",header = T,row.names = 1)
GSE94325_samples <- colnames(GSE94325_fpkm)
GSE94325_genes <- row.names(GSE94325_fpkm)
GSE94325_logfpkm <- log2(GSE94325_fpkm+1)

# batch effect correction
common_genes <- intersect(intersect(intersect(intersect(all_genes_1,GSE66582_genes),GSE85632_genes),all_genes_2),GSE94325_genes)
data <- cbind(logfpkm_1[common_genes,],GSE66582_logfpkm[common_genes,],GSE85632_logfpkm[common_genes,],logfpkm_2[common_genes,],GSE94325_logfpkm[common_genes,])
experiment_batch <- as.factor(c(rep(1,ncol(logfpkm_1)),rep(2,ncol(GSE66582_logfpkm)),rep(3,ncol(GSE85632_logfpkm)),rep(4,ncol(logfpkm_2)),rep(5,ncol(GSE94325_logfpkm))))
bec_data <- removeBatchEffect(data,batch=experiment_batch)

###################################################################################################
################################                PLOT               ################################
###################################################################################################
# inhouse data correlation and PCA
timepoints <- c("dox_C","dox_E","pu_C","pu_E","pup5_E","duxpu_C","duxpu_E")

dox_C <- all_samples[grep("^dox_C",all_samples)]
dox_E <- all_samples[grep("^dox_E",all_samples)]
pu_C <- all_samples[grep("^pu_C",all_samples)]
pu_E <- all_samples[grep("^pu_E",all_samples)]
pup5_E <- all_samples[grep("^pup5_E",all_samples)]
duxpu_C <- all_samples[grep("^duxpu_C",all_samples)]
duxpu_E <- all_samples[grep("^duxpu_E",all_samples)]

inhouse_path <- c(dox_C,dox_E,pu_C,pu_E,pup5_E,duxpu_C,duxpu_E)

COLS <- rep("black",length(inhouse_path))
names(COLS) <- inhouse_path
COLS[dox_C] <- cccol[1]
COLS[dox_E] <- cccol[2]
COLS[pu_C] <- cccol[3]
COLS[pu_E] <- cccol[4]
COLS[pup5_E] <- cccol[5]
COLS[duxpu_C] <- cccol[6]
COLS[duxpu_E] <- cccol[7]
PCHS <- rep(1,length(inhouse_path))
names(PCHS) <- inhouse_path
PCHS[inhouse_path] <- 19

# Covered gene number
pdf("inhousedata.logfpkm.CGN.pdf",width=6,height=6)
par(mar=c(9,4,2,2))
covered_gene_number <- apply(fpkm[,inhouse_path],2,getcoverGnum)
barplot(covered_gene_number,las=2,col=COLS,border=COLS)
dev.off()

# box plot of all the expression values
pdf("inhousedata.logfpkm.box.pdf",width=8,height=8)
boxplot(logfpkm[,inhouse_path],las=2,outline=F)
dev.off()

# PCA of all the expression values
pdf("inhousedata.logfpkm.pca.pdf",width=10,height=10)
pca <- prcomp(t(logfpkm[,inhouse_path]))
plot(pca$x[,"PC1"],pca$x[,"PC2"],pch=PCHS,col=COLS,main="PCA",xlab="PC1",ylab="PC2")
text(pca$x[,"PC1"],pca$x[,"PC2"]+1,inhouse_path,col=COLS,cex=0.3)
dev.off()

# MDS of all the expression values
d <- dist(t(logfpkm[,inhouse_path])) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
pdf("inhousedata.logfpkm.MDS.pdf",width=6,height=6)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS",col=COLS,pch=PCHS,lwd=4)
legend("center", timepoints, col=cccol,pch=19,bty="n",ncol = 2)
dev.off()

# correlation of all the samples
sample_cor <- cor(fpkm[,inhouse_path])
pdf("inhousedata.fpkm.cor.pdf",width=8,height=8)
ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"), bias=1)(100)
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",margins=c(10,10), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
cor_cell <- matrix(as.character(round(sample_cor, 2)), ncol=dim(sample_cor)[2])
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",cellnote=cor_cell, notecol="black", notecex=0.5, margins=c(10,10), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
dev.off()
pdf("inhousedata.fpkm.cor.legend.pdf",width=3,height=1.8)
par(mar=c(4,2,4,2))
plotMatrix <- sample_cor
all_exp <- as.matrix(plotMatrix) # using same scale bar
zmax <- max(na.omit(all_exp))
zmin <- min(na.omit(all_exp))
ColorLevels <- seq(to=zmax,from=zmin, length=1000)   #number sequence
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="Correlation cofficient",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,1),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()


# all data (before bec) correlation and PCA
all_samples <- colnames(data)
timepoints <- c("dox_C","dox_E","pu_C","pu_E","pup5_E","duxpu_C","duxpu_E","GSE66582_2C","GSE66582_E2C","GSE66582_ESCs","GSE85632_2CLC","GSE85632_ESCs","GSE94325_ESCs")
dox_C <- all_samples[grep("^dox_C",all_samples)]
dox_E <- all_samples[grep("^dox_E",all_samples)]
pu_C <- all_samples[grep("^pu_C",all_samples)]
pu_E <- all_samples[grep("^pu_E",all_samples)]
pup5_E <- all_samples[grep("^pup5_E",all_samples)]
duxpu_C <- all_samples[grep("^duxpu_C",all_samples)]
duxpu_E <- all_samples[grep("^duxpu_E",all_samples)]
GSE66582_2C <- all_samples[grep("^GSE66582_2C",all_samples)]
GSE66582_E2C <- all_samples[grep("^GSE66582_E2C",all_samples)]
GSE66582_ESCs <- all_samples[grep("^GSE66582_ESCs",all_samples)]
GSE85632_2CLC <- all_samples[grep("^GSE85632_2CLC",all_samples)]
GSE85632_ESCs <- all_samples[grep("^GSE85632_ESCs",all_samples)]
GSE94325_ESCs <- all_samples[grep("^GSE94325_ESCs",all_samples)]

# all_path <- c(dox_C,dox_E,pu_C,pu_E,pup5_E,duxpu_C,duxpu_E,GSE66582_2C,GSE66582_E2C,GSE66582_ESCs,GSE85632_2CLC,GSE85632_ESCs,GSE94325_ESCs)
all_path <- c(dox_C,dox_E,pu_C,pu_E,pup5_E,duxpu_C,duxpu_E,GSE66582_2C,GSE85632_2CLC,GSE85632_ESCs)

COLS <- rep("black",length(all_path))
names(COLS) <- all_path
COLS[dox_C] <- cccol[1]
COLS[dox_E] <- cccol[2]
COLS[pu_C] <- cccol[3]
COLS[pu_E] <- cccol[4]
COLS[pup5_E] <- cccol[5]
COLS[duxpu_C] <- cccol[6]
COLS[duxpu_E] <- cccol[7]
COLS[GSE66582_2C] <- cccol[8]
# COLS[GSE66582_E2C] <- cccol[9]
# COLS[GSE66582_ESCs] <- cccol[10]
COLS[GSE85632_2CLC] <- cccol[11]
COLS[GSE85632_ESCs] <- cccol[12]
# COLS[GSE94325_ESCs] <- cccol[13]
PCHS <- rep(1,length(all_path))
names(PCHS) <- all_path
PCHS[all_path] <- 19

# Covered gene number
pdf("alldataBeforebec.logfpkm.CGN.pdf",width=6,height=6)
par(mar=c(9,4,2,2))
covered_gene_number <- apply(data[,all_path],2,getcoverGnum)
barplot(covered_gene_number,las=2,col=COLS,border=COLS)
dev.off()

# box plot of all the expression values
pdf("alldataBeforebec.logfpkm.box.pdf",width=8,height=8)
boxplot(data[,all_path],las=2,outline=F)
dev.off()

# PCA of all the expression values
pdf("alldataBeforebec.logfpkm.pca.pdf",width=10,height=10)
pca <- prcomp(t(data[,all_path]))
plot(pca$x[,"PC1"],pca$x[,"PC2"],pch=PCHS,col=COLS,main="PCA",xlab="PC1",ylab="PC2")
text(pca$x[,"PC1"],pca$x[,"PC2"]+1,all_path,col=COLS,cex=0.3)
dev.off()

# MDS of all the expression values
d <- dist(t(data[,all_path])) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
pdf("alldataBeforebec.logfpkm.MDS.pdf",width=6,height=6)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS",col=COLS,pch=PCHS,lwd=4)
legend("center", timepoints, col=cccol,pch=19,bty="n",ncol = 2)
dev.off()

# correlation of all the samples
sample_cor <- cor(data[,all_path])
pdf("alldataBeforebec.fpkm.cor.pdf",width=8,height=8)
ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"), bias=1)(100)
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",margins=c(10,10), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
cor_cell <- matrix(as.character(round(sample_cor, 2)), ncol=dim(sample_cor)[2])
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",cellnote=cor_cell, notecol="black", notecex=0.5, margins=c(10,10), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
dev.off()
pdf("alldataBeforebec.fpkm.cor.legend.pdf",width=3,height=1.8)
par(mar=c(4,2,4,2))
plotMatrix <- sample_cor
all_exp <- as.matrix(plotMatrix) # using same scale bar
zmax <- max(na.omit(all_exp))
zmin <- min(na.omit(all_exp))
ColorLevels <- seq(to=zmax,from=zmin, length=1000)   #number sequence
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="Correlation cofficient",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,1),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()

# all data (after bec) correlation and PCA
# Covered gene number
pdf("alldataAfterbec.logfpkm.CGN.pdf",width=6,height=6)
par(mar=c(9,4,2,2))
covered_gene_number <- apply(bec_data[,all_path],2,getcoverGnum)
barplot(covered_gene_number,las=2,col=COLS,border=COLS)
dev.off()

# box plot of all the expression values
pdf("alldataAfterbec.logfpkm.box.pdf",width=8,height=8)
boxplot(bec_data[,all_path],las=2,outline=F)
dev.off()

# PCA of all the expression values
pdf("alldataAfterbec.logfpkm.pca.pdf",width=10,height=10)
pca <- prcomp(t(bec_data[,all_path]))
plot(pca$x[,"PC1"],pca$x[,"PC2"],pch=PCHS,col=COLS,main="PCA",xlab="PC1",ylab="PC2")
text(pca$x[,"PC1"],pca$x[,"PC2"]+1,all_path,col=COLS,cex=0.3)
dev.off()

# MDS of all the expression values
d <- dist(t(bec_data[,all_path])) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
pdf("alldataAfterbec.logfpkm.MDS.pdf",width=6,height=6)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS",col=COLS,pch=PCHS,lwd=4)
legend("center", timepoints, col=cccol,pch=19,bty="n",ncol = 2)
dev.off()

# correlation of all the samples
sample_cor <- cor(bec_data[,all_path])
pdf("alldataAfterbec.fpkm.cor.pdf",width=8,height=8)
ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"), bias=1)(100)
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",margins=c(10,10), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
cor_cell <- matrix(as.character(round(sample_cor, 2)), ncol=dim(sample_cor)[2])
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",cellnote=cor_cell, notecol="black", notecex=0.5, margins=c(10,10), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
dev.off()
pdf("alldataAfterbec.fpkm.cor.legend.pdf",width=3,height=1.8)
par(mar=c(4,2,4,2))
plotMatrix <- sample_cor
all_exp <- as.matrix(plotMatrix) # using same scale bar
zmax <- max(na.omit(all_exp))
zmin <- min(na.omit(all_exp))
ColorLevels <- seq(to=zmax,from=zmin, length=1000)   #number sequence
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="Correlation cofficient",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,1),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()
###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################

############ Expression matrix
fpkm_0 <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/WQN_DUXno_mmESCs_GeneSymbol_fpkm.txt",header=T,row.names=1,sep="\t")
fpkm_1 <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/WQN_duxpu_mmESCs_GeneSymbol_fpkm.txt",header=T,row.names=1)
fpkm_2 <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/data/GSE94325_mmESC_RNAseq/GSE94325_GeneSymbol_fpkm.txt",header = T,row.names = 1)
fpkm_3 <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/data/GSE66582_mmEmbryo_RNAseq/GSE66582_GeneSymbol_fpkm.txt",header = T,row.names=1)
fpkm_4 <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/data/GSE85632_mmESC_RNAseq/GSE85632_GeneSymbol_fpkm.txt",header = T,row.names = 1)
common_genes <- intersect(intersect(intersect(intersect(rownames(fpkm_1),rownames(fpkm_2)),rownames(fpkm_3)),rownames(fpkm_4)),rownames(fpkm_0))
fpkm <- cbind(fpkm_0[common_genes,],fpkm_1[common_genes,],fpkm_2[common_genes,],fpkm_3[common_genes,],fpkm_4[common_genes,])
all_samples <- colnames(fpkm)
all_genes <- rownames(fpkm)
logfpkm <- log2(fpkm+1)

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
GSE85632_ESCs <- all_samples[grep("^GSE85632_ESCs",all_samples)]
GSE85632_2CLC <- all_samples[grep("^GSE85632_2CLC",all_samples)]
GSE94325_ESCs <- all_samples[grep("^GSE94325_ESCs",all_samples)]
# path <- c(dox_C,dox_E,pu_C,pu_E,pup5_E,duxpu_C,duxpu_E,GSE66582_2C,GSE66582_E2C,GSE66582_ESCs,GSE85632_ESCs,GSE85632_2CLC,GSE94325_ESCs)
path <- c(dox_C,dox_E,pu_C,pu_E,pup5_E,duxpu_C,duxpu_E,GSE66582_2C,GSE66582_E2C,GSE85632_2CLC,GSE85632_ESCs)
timepoint_ESC <- c("dox_C","dox_E","pu_C","pu_E","pup5_E","duxpu_C","duxpu_E","GSE66582_2C","GSE66582_E2C","GSE85632_2CLC","GSE85632_ESCs")


############ development 
fpkm_dev <- read.table("/data/zhaochengchen/Work/1.CSwitch/results/0.AverageSignal/GSE136714_RNAseqFPKM_Combined.txt",header=T,row.names=1)
all_cells <- colnames(fpkm_dev)
logfpkm_dev <- log2(fpkm_dev+1)

Zygote <- all_cells[grep("^Zygote",all_cells)]
C2 <- all_cells[grep("^C2([^0-9]|$)",all_cells)]
C4 <- all_cells[grep("^C4([^0-9]|$)",all_cells)]
LC4 <- all_cells[grep("^LC4([^0-9]|$)",all_cells)]
C8 <- all_cells[grep("^C8([^0-9]|$)",all_cells)]
C16 <- all_cells[grep("^C16([^0-9]|$)",all_cells)]
C32 <- all_cells[grep("^C32([^0-9]|$)",all_cells)]

avgfpkm_dev <- cbind(apply(fpkm_dev[,Zygote],1,mean),apply(fpkm_dev[,C2],1,mean),apply(fpkm_dev[,C4],1,mean),apply(fpkm_dev[,LC4],1,mean),apply(fpkm_dev[,C8],1,mean),apply(fpkm_dev[,C16],1,mean),apply(fpkm_dev[,C32],1,mean))
timepoint_dev <- c("Zygote","C2","C4","LC4","C8","C16","C32")
dev_labels <- c("Zygote","C2","C4","LC4","C8","C16","C32")
colnames(avgfpkm_dev) <- timepoint_dev
development_path <- timepoint_dev
logavgfpkm_dev <- log2(avgfpkm_dev+1)

all_path <- c(path,Zygote,C2,C4,LC4,C8,C16,C32)
common_genes <- intersect(common_genes, rownames(fpkm_dev))
data <- cbind(logfpkm[common_genes,path],logfpkm_dev[common_genes,])

# 2021- genes enriched in 2-cell embryos and the pluripotent marker genes used for gene enrichment analysis.xlsx
Totipotency_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/Totipotency_genes.txt")[,1])
Pluripotent_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/Pluripotent_genes.txt")[,1])

common_genes <- intersect(c(Totipotency_genes, Pluripotent_genes),common_genes)

###################################################################################################
################################                PLOT               ################################
###################################################################################################

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
COLS[GSE66582_E2C] <- cccol[9]
COLS[GSE85632_2CLC] <- cccol[10]
COLS[GSE85632_ESCs] <- cccol[11]
# COLS[GSE94325_ESCs] <- cccol[12]
COLS[Zygote] <- cccol[13]
COLS[C2] <- cccol[14]
COLS[C4] <- cccol[15]
COLS[LC4] <- cccol[16]
COLS[C8] <- cccol[17]
COLS[C16] <- cccol[18]
COLS[C32] <- cccol[19]
PCHS <- rep(1,length(all_path))
names(PCHS) <- all_path
PCHS[all_path] <- 19

# PCA of all the expression values
pdf("AlldataDevelop.logfpkm.pca.TPG.pdf",width=10,height=10)
pca <- prcomp(t(data[common_genes,all_path]))
plot(pca$x[,"PC1"],pca$x[,"PC2"],pch=PCHS,col=COLS,main="PCA",xlab="PC1",ylab="PC2")
text(pca$x[,"PC1"],pca$x[,"PC2"]+1,all_path,col=COLS,cex=0.3)
dev.off()

# MDS of all the expression values
d <- dist(t(data[common_genes,all_path])) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
pdf("AlldataDevelop.logfpkm.MDS.TPG.pdf",width=6,height=6)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS",col=COLS,pch=PCHS,lwd=4)
plot(1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("center", c(timepoint_ESC,timepoint_dev), col=cccol,pch=19,bty="n",ncol = 2)
dev.off()

############################################################################
#################    development.correlation.with.ESC   ####################
############################################################################

n_cor <- c()
for (each_sample in path){
	tmp_cor <- c()
	for (each in rev(timepoint_dev)){
		tmp_cor <- c(tmp_cor,cor(logavgfpkm_dev[common_genes,each],logfpkm[common_genes,each_sample]))
	}
	n_cor <- cbind(n_cor,tmp_cor)
}
rownames(n_cor) <- rev(timepoint_dev)
colnames(n_cor) <- path
plot_matrix <- n_cor

pdf("Fig.DevelopmentCorrelation.TPG.pdf",width=14,height=9)
all_exp <- c(as.matrix(plot_matrix))
zmax <- max(all_exp)
zmin <- min(all_exp)
par(mar=c(8,5,2,2))
layout(matrix(c(1,1,1,1,1,1,1,1,1,2),nrow=10,ncol=1,byrow=F))
ColorRamp <- colorRampPalette(c("white",cccol[2]), bias=1)(100)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=100)   #number sequence
image(1:ncol(plot_matrix), 1:nrow(plot_matrix), t(plot_matrix), axes=F, col=ColorRamp, xlab="",ylab="");box(lwd=2)
axis(side=1,1:ncol(plot_matrix),labels=path,las=2)
axis(side=2,1:nrow(plot_matrix),labels=rev(dev_labels),las=2)
text(matrix(rep(1:ncol(plot_matrix),nrow(plot_matrix)),ncol=ncol(plot_matrix),nrow=nrow(plot_matrix),byrow = T), matrix(rep(1:nrow(plot_matrix),ncol(plot_matrix)),ncol=ncol(plot_matrix),nrow=nrow(plot_matrix)),as.matrix(round(plot_matrix,2)))
par(mar=c(4,10,2,10))
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="log2(FPKM+1)",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,2),round(zmax,1)))
dev.off()


# ############################################################################
# #################                normalise              ####################
# ############################################################################

# # normalise
# library(edgeR)
# experiment_batch <- as.factor(c(rep(1,length(c(dox_C,dox_E,pu_C,pu_E,pup5_E,duxpu_C,duxpu_E))),rep(2,length(c(GSE66582_2C,GSE66582_E2C))),rep(3,length(c(GSE85632_2CLC,GSE85632_ESCs))),rep(4,length(c(GSE94325_ESCs))),rep(5,length(c(Zygote,C2,C4,LC4,C8,C16,C32)))))
# rmbatch_data <- removeBatchEffect(data,batch=experiment_batch)

# avgfpkm_dev <- cbind(apply(rmbatch_data[,Zygote],1,mean),apply(rmbatch_data[,C2],1,mean),apply(rmbatch_data[,C4],1,mean),apply(rmbatch_data[,LC4],1,mean),apply(rmbatch_data[,C8],1,mean),apply(rmbatch_data[,C16],1,mean),apply(rmbatch_data[,C32],1,mean))
# timepoint_dev <- c("Zygote","C2","C4","LC4","C8","C16","C32")
# dev_labels <- c("Zygote","C2","C4","LC4","C8","C16","C32")
# colnames(avgfpkm_dev) <- timepoint_dev

# logfpkm <- rmbatch_data[,path]

# COLS <- rep("black",length(all_path))
# names(COLS) <- all_path
# COLS[dox_C] <- cccol[1]
# COLS[dox_E] <- cccol[2]
# COLS[pu_C] <- cccol[3]
# COLS[pu_E] <- cccol[4]
# COLS[pup5_E] <- cccol[5]
# COLS[duxpu_C] <- cccol[6]
# COLS[duxpu_E] <- cccol[7]
# COLS[GSE66582_2C] <- cccol[8]
# COLS[GSE66582_E2C] <- cccol[9]
# COLS[GSE85632_2CLC] <- cccol[10]
# COLS[GSE85632_ESCs] <- cccol[11]
# COLS[GSE94325_ESCs] <- cccol[12]
# COLS[Zygote] <- cccol[13]
# COLS[C2] <- cccol[14]
# COLS[C4] <- cccol[15]
# COLS[LC4] <- cccol[16]
# COLS[C8] <- cccol[17]
# COLS[C16] <- cccol[18]
# COLS[C32] <- cccol[19]
# PCHS <- rep(1,length(all_path))
# names(PCHS) <- all_path
# PCHS[all_path] <- 19

# # PCA of all the expression values
# pdf("Allrmbatch_dataDevelop.logfpkm.pca.TPG.pdf",width=10,height=10)
# pca <- prcomp(t(rmbatch_data[,all_path]))
# plot(pca$x[,"PC1"],pca$x[,"PC2"],pch=PCHS,col=COLS,main="PCA",xlab="PC1",ylab="PC2")
# text(pca$x[,"PC1"],pca$x[,"PC2"]+1,all_path,col=COLS,cex=0.3)
# dev.off()

# # MDS of all the expression values
# d <- dist(t(rmbatch_data[,all_path])) # euclidean distances between the rows
# fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# # plot solution 
# x <- fit$points[,1]
# y <- fit$points[,2]
# pdf("Allrmbatch_dataDevelop.logfpkm.MDS.TPG.pdf",width=6,height=6)
# plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS",col=COLS,pch=PCHS,lwd=4)
# plot(1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
# legend("center", c(timepoint_ESC,timepoint_dev), col=cccol,pch=19,bty="n",ncol = 2)
# dev.off()

# ############################################################################
# #################    development.correlation.with.ESC   ####################
# ############################################################################

# n_cor <- c()
# for (each_sample in path){
# 	tmp_cor <- c()
# 	for (each in rev(timepoint_dev)){
# 		tmp_cor <- c(tmp_cor,cor(logavgfpkm_dev[common_genes,each],logfpkm[common_genes,each_sample]))
# 	}
# 	n_cor <- cbind(n_cor,tmp_cor)
# }
# rownames(n_cor) <- rev(timepoint_dev)
# colnames(n_cor) <- path
# plot_matrix <- n_cor

# pdf("Fig.DevelopmentCorrelation.rmbatch_data.TPG.pdf",width=14,height=9)
# all_exp <- c(as.matrix(plot_matrix))
# zmax <- max(all_exp)
# zmin <- min(all_exp)
# par(mar=c(8,5,2,2))
# layout(matrix(c(1,1,1,1,1,1,1,1,1,2),nrow=10,ncol=1,byrow=F))
# ColorRamp <- colorRampPalette(c("white",cccol[2]), bias=1)(100)   #color list
# ColorLevels <- seq(to=zmax,from=zmin, length=100)   #number sequence
# image(1:ncol(plot_matrix), 1:nrow(plot_matrix), t(plot_matrix), axes=F, col=ColorRamp, xlab="",ylab="");box(lwd=2)
# axis(side=1,1:ncol(plot_matrix),labels=path,las=2)
# axis(side=2,1:nrow(plot_matrix),labels=rev(dev_labels),las=2)
# text(matrix(rep(1:ncol(plot_matrix),nrow(plot_matrix)),ncol=ncol(plot_matrix),nrow=nrow(plot_matrix),byrow = T), matrix(rep(1:nrow(plot_matrix),ncol(plot_matrix)),ncol=ncol(plot_matrix),nrow=nrow(plot_matrix)),as.matrix(round(plot_matrix,2)))
# par(mar=c(4,10,2,10))
# image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="log2(FPKM+1)",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
# axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,2),round(zmax,1)))
# dev.off()



############################################################################
#################            selected samples           ####################
############################################################################

selected_path <- c(pu_C,pu_E,duxpu_C,dox_C,Zygote,C2,C4,LC4,C8,C16,C32)

# PCA of all the expression values
pdf("SelecteddataDevelop.logfpkm.pca.TPG.pdf",width=10,height=10)
pca <- prcomp(t(data[common_genes,selected_path]))
plot(pca$x[,"PC1"],pca$x[,"PC2"],pch=PCHS[selected_path],col=COLS[selected_path],main="PCA",xlab="PC1",ylab="PC2")
text(pca$x[,"PC1"],pca$x[,"PC2"]+1,selected_path,col=COLS[selected_path],cex=0.3)
dev.off()

# MDS of all the expression values
d <- dist(t(data[common_genes,selected_path])) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
pdf("SelecteddataDevelop.logfpkm.MDS.TPG.pdf",width=6,height=6)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS",col=COLS[selected_path],pch=PCHS[selected_path],lwd=4)
plot(1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("center", c("pu_C","pu_E","duxpu_C","dox_C","Zygote","C2","C4","LC4","C8","C16","C32"), col=cccol[c(3,4,6,1,13,14,15,16,17,18,19)],pch=19,bty="n",ncol = 2)
dev.off()

############################################################################
#################    development.correlation.with.ESC   ####################
############################################################################
selected_path <- c(pu_C,pu_E,duxpu_C,dox_C)

n_cor <- c()
for (each_sample in selected_path){
	tmp_cor <- c()
	for (each in rev(timepoint_dev)){
		tmp_cor <- c(tmp_cor,cor(logavgfpkm_dev[common_genes,each],logfpkm[common_genes,each_sample]))
	}
	n_cor <- cbind(n_cor,tmp_cor)
}
rownames(n_cor) <- rev(timepoint_dev)
colnames(n_cor) <- selected_path
plot_matrix <- n_cor

pdf("Fig.SelectedDevelopmentCorrelation.TPG.pdf",width=8,height=9)
all_exp <- c(as.matrix(plot_matrix))
zmax <- max(all_exp)
zmin <- min(all_exp)
par(mar=c(8,5,2,2))
layout(matrix(c(1,1,1,1,1,1,1,1,1,2),nrow=10,ncol=1,byrow=F))
ColorRamp <- colorRampPalette(c("white",cccol[2]), bias=1)(100)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=100)   #number sequence
image(1:ncol(plot_matrix), 1:nrow(plot_matrix), t(plot_matrix), axes=F, col=ColorRamp, xlab="",ylab="");box(lwd=2)
axis(side=1,1:ncol(plot_matrix),labels=selected_path,las=2)
axis(side=2,1:nrow(plot_matrix),labels=rev(dev_labels),las=2)
text(matrix(rep(1:ncol(plot_matrix),nrow(plot_matrix)),ncol=ncol(plot_matrix),nrow=nrow(plot_matrix),byrow = T), matrix(rep(1:nrow(plot_matrix),ncol(plot_matrix)),ncol=ncol(plot_matrix),nrow=nrow(plot_matrix)),as.matrix(round(plot_matrix,2)))
par(mar=c(4,10,2,10))
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="log2(FPKM+1)",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,2),round(zmax,1)))
dev.off()


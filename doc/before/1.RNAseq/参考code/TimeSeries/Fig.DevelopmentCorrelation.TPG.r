###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
red_col <- colorRampPalette(c(cccol[1],cccol[3]), bias=1)(5)[c(1,3,4,5)]
blue_col <- colorRampPalette(c(cccol[1],cccol[2]), bias=1)(5)[c(1,3,4,5)]
green_col <- colorRampPalette(c(cccol[4],cccol[6]), bias=1)(7)
grey_col <- colorRampPalette(c(cccol[8],"black"), bias=1)(4)
options(scipen = 200)

library("edgeR")
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
fpkm <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/WQN_20220907_mmESCs_GeneSymbol_fpkm.txt",header=T,row.names=1,sep="\t")
all_genes <- rownames(fpkm)
all_samples <- colnames(fpkm)
logfpkm <- log2(fpkm+1)

timepoints <- c("H0","H18","H30","H43")

mESC_0h <- all_samples[grep("^mESC_0h",all_samples)]
mESC_p18h <- all_samples[grep("^mESC_p18h",all_samples)]
mESC_p30h <- all_samples[grep("^mESC_p30h",all_samples)]
mESC_p42h <- all_samples[grep("^mESC_p42h",all_samples)]
tCLC_p18h <- all_samples[grep("^tCLC_p18h",all_samples)]
tCLC_p30h <- all_samples[grep("^tCLC_p30h",all_samples)]
tCLC_p42h <- all_samples[grep("^tCLC_p42h",all_samples)]

mESC_path <- c(mESC_0h,mESC_p18h,mESC_p30h,mESC_p42h)
tCLC_path <- c(mESC_0h,tCLC_p18h,tCLC_p30h,tCLC_p42h)

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

# 2021- genes enriched in 2-cell embryos and the pluripotent marker genes used for gene enrichment analysis.xlsx
Totipotency_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/Totipotency_genes.txt")[,1])
Pluripotent_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/Pluripotent_genes.txt")[,1])

common_genes <- intersect(all_genes, rownames(fpkm_dev))
TPG_genes <- intersect(c(Totipotency_genes, Pluripotent_genes),common_genes)
all_path <- c(mESC_0h,mESC_p18h,mESC_p30h,mESC_p42h,tCLC_p18h,tCLC_p30h,tCLC_p42h,Zygote,C2,C4,LC4,C8,C16,C32)
data <- cbind(logfpkm[common_genes,c(mESC_0h,mESC_p18h,mESC_p30h,mESC_p42h,tCLC_p18h,tCLC_p30h,tCLC_p42h)],logfpkm_dev[common_genes,c(Zygote,C2,C4,LC4,C8,C16,C32)])

# # batch effect correction
# experiment_batch <- as.factor(c(rep(1,ncol(logfpkm)),rep(2,length(c(Zygote,C2,C4,LC4,C8,C16,C32)))))
# bec_data <- removeBatchEffect(data,batch=experiment_batch)
###################################################################################################
################################                PLOT               ################################
###################################################################################################

COLS <- rep("black",length(all_path))
names(COLS) <- all_path
COLS[mESC_0h] <- red_col[1]
COLS[mESC_p18h] <- red_col[2]
COLS[mESC_p30h] <- red_col[3]
COLS[mESC_p42h] <- red_col[4]
COLS[tCLC_p18h] <- blue_col[2]
COLS[tCLC_p30h] <- blue_col[3]
COLS[tCLC_p42h] <- blue_col[4]
COLS[Zygote] <- green_col[1]
COLS[C2] <- green_col[2]
COLS[C4] <- green_col[3]
COLS[LC4] <- green_col[4]
COLS[C8] <- green_col[5]
COLS[C16] <- green_col[6]
COLS[C32] <- green_col[7]

PCHS <- rep(18,length(all_path))
names(PCHS) <- all_path
PCHS[tCLC_path] <- 15
PCHS[mESC_path] <- 19

pdf("Fig.DevelopmentCorrelationPCA.TPG.pdf",width=6,height=6)
pca <- prcomp(t(data[TPG_genes,]))
plot(pca$x[,"PC1"],pca$x[,"PC2"],pch=PCHS,col=COLS,main="PCA",xlab="PC1",ylab="PC2")
# text(pca$x[,"PC1"],pca$x[,"PC2"],all_path,pch=PCHS,col=COLS,cex=0.3)
plot(1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("center", c("mESC_0h","mESC_p18h","mESC_p30h","mESC_p42h","tCLC_p18h","tCLC_p30h","tCLC_p42h","Zygote","C2","C4","LC4","C8","C16","C32"), col=c(red_col[1],red_col[2],red_col[3],red_col[4],blue_col[2],blue_col[3],blue_col[4],green_col),
	pch=c(rep(19,4),rep(15,3),rep(18,7)),bty="n",ncol = 2)
dev.off()

d <- dist(t(data[TPG_genes,])) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
pdf("Fig.DevelopmentCorrelationMDS.TPG.pdf",width=6,height=6)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS",col=COLS,pch=PCHS,lwd=2)
plot(1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("center", c("mESC_0h","mESC_p18h","mESC_p30h","mESC_p42h","tCLC_p18h","tCLC_p30h","tCLC_p42h","Zygote","C2","C4","LC4","C8","C16","C32"), col=c(red_col[1],red_col[2],red_col[3],red_col[4],blue_col[2],blue_col[3],blue_col[4],green_col),
	pch=c(rep(19,4),rep(15,3),rep(18,7)),bty="n",ncol = 2)
dev.off()










# all_path <- c(mESC_0h,mESC_p18h,mESC_p30h,mESC_p42h,tCLC_p18h,tCLC_p30h,tCLC_p42h,timepoint_dev)
# data <- cbind(logfpkm[common_genes,c(mESC_0h,mESC_p18h,mESC_p30h,mESC_p42h,tCLC_p18h,tCLC_p30h,tCLC_p42h)],logavgfpkm_dev[common_genes,timepoint_dev])

# ###################################################################################################
# ################################                PLOT               ################################
# ###################################################################################################

# COLS <- rep("black",length(all_path))
# names(COLS) <- all_path
# COLS[mESC_0h] <- red_col[1]
# COLS[mESC_p18h] <- red_col[2]
# COLS[mESC_p30h] <- red_col[3]
# COLS[mESC_p42h] <- red_col[4]
# COLS[tCLC_p18h] <- blue_col[2]
# COLS[tCLC_p30h] <- blue_col[3]
# COLS[tCLC_p42h] <- blue_col[4]
# COLS["Zygote"] <- cccol[4]
# COLS["C2"] <- cccol[5]
# COLS["C4"] <- cccol[6]
# COLS["LC4"] <- cccol[7]
# COLS["C8"] <- cccol[8]
# COLS["C16"] <- cccol[9]
# COLS["C32"] <- cccol[10]

# PCHS <- rep(8,length(all_path))
# names(PCHS) <- all_path
# PCHS[tCLC_path] <- 15
# PCHS[mESC_path] <- 19

# pdf("inhousedata.logfpkm.pca.pdf",width=6,height=6)
# pca <- prcomp(t(data[TPG_genes,]))
# plot(pca$x[,"PC1"],pca$x[,"PC2"],pch=PCHS,col=COLS,main="PCA",xlab="PC1",ylab="PC2")
# # text(pca$x[,"PC1"],pca$x[,"PC2"],all_path,pch=PCHS,col=COLS,cex=0.3)
# dev.off()

# d <- dist(t(data[TPG_genes,])) # euclidean distances between the rows
# fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# # plot solution 
# x <- fit$points[,1]
# y <- fit$points[,2]
# pdf("Fig.DevelopmentCorrelationMDS.TPG.pdf",width=6,height=6)
# plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS",col=COLS,pch=PCHS,lwd=2)
# plot(1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
# legend("center", c("mESC_0h","mESC_p18h","mESC_p30h","mESC_p42h","tCLC_p18h","tCLC_p30h","tCLC_p42h","Zygote","C2","C4","LC4","C8","C16","C32"), col=c(red_col[1],red_col[2],red_col[3],red_col[4],blue_col[2],blue_col[3],blue_col[4],cccol[4],cccol[5],cccol[6],cccol[7],cccol[8],cccol[9],cccol[10]),
# 	pch=c(rep(19,4),rep(15,3),rep(8,7)),bty="n",ncol = 2)
# dev.off()





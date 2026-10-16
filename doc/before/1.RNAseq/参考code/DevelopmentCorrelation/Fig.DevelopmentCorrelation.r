###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)
library(quantreg) # https://www.codenong.com/31862766/
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################

############ Expression matrix
fpkm <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/WQN_DUXno_mmESCs_GeneSymbol_fpkm.txt",header=T,row.names=1,sep="\t")
all_samples <- colnames(fpkm)
all_genes <- rownames(fpkm)
logfpkm <- log2(fpkm+1)

dox_C <- all_samples[grep("^dox_C",all_samples)]
dox_E <- all_samples[grep("^dox_E",all_samples)]
pu_C <- all_samples[grep("^pu_C",all_samples)]
pu_E <- all_samples[grep("^pu_E",all_samples)]
pup5_E <- all_samples[grep("^pup5_E",all_samples)]

path <- c(dox_C,dox_E,pu_C,pu_E,pup5_E)
timepoint_ESC <- c("dox_C","dox_E","pu_C","pu_E","pup5_E")


############ development 
fpkm_dev <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/data/GSE66582_mmEmbryo_RNAseq/GSE66582_GeneSymbol_fpkm.txt",header=T,row.names=1)

all_samples_dev <- colnames(fpkm_dev)
all_genes_dev <- rownames(fpkm_dev)
logfpkm_dev <- log2(fpkm_dev+1)
C2 <- all_samples_dev[grep("^GSE66582_2C([^0-9]|$)",all_samples_dev)]

common_genes <- intersect(all_genes, all_genes_dev)

# 2021- genes enriched in 2-cell embryos and the pluripotent marker genes used for gene enrichment analysis.xlsx
Totipotency_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.DiffPattern/RNAseq/Totipotency_genes.txt")[,1])
Pluripotent_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.DiffPattern/RNAseq/Pluripotent_genes.txt")[,1])

TPG <- intersect(c(Totipotency_genes, Pluripotent_genes),common_genes)


###################################################################################################
################################                PLOT               ################################
###################################################################################################

pdf("Fig.DevelopmentCorrelation.EachSample.pdf",width=8,height=8)
par(mar=c(4,4,2,2),mfrow=c(4,4))

x_ALL <- logfpkm[common_genes,"dox_E_r1"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r1"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="dox_E_r1 log2(fpkm+1)", ylab="2C_r1 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"dox_E_r2"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r1"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="dox_E_r2 log2(fpkm+1)", ylab="2C_r1 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"dox_E_r1"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r2"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="dox_E_r1 log2(fpkm+1)", ylab="2C_r2 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"dox_E_r2"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r2"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="dox_E_r2 log2(fpkm+1)", ylab="2C_r2 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])

x_ALL <- logfpkm[common_genes,"dox_C_r1"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r1"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="dox_C_r1 log2(fpkm+1)", ylab="2C_r1 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"dox_C_r2"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r1"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="dox_C_r2 log2(fpkm+1)", ylab="2C_r1 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"dox_C_r1"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r2"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="dox_C_r1 log2(fpkm+1)", ylab="2C_r2 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"dox_C_r2"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r2"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="dox_C_r2 log2(fpkm+1)", ylab="2C_r2 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])

x_ALL <- logfpkm[common_genes,"pu_E_r1"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r1"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="pu_E_r1 log2(fpkm+1)", ylab="2C_r1 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"pu_E_r2"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r1"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="pu_E_r2 log2(fpkm+1)", ylab="2C_r1 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"pu_E_r1"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r2"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="pu_E_r1 log2(fpkm+1)", ylab="2C_r2 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"pu_E_r2"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r2"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="pu_E_r2 log2(fpkm+1)", ylab="2C_r2 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])

x_ALL <- logfpkm[common_genes,"pu_C_r1"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r1"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="pu_C_r1 log2(fpkm+1)", ylab="2C_r1 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"pu_C_r2"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r1"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="pu_C_r2 log2(fpkm+1)", ylab="2C_r1 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"pu_C_r1"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r2"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="pu_C_r1 log2(fpkm+1)", ylab="2C_r2 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])
x_ALL <- logfpkm[common_genes,"pu_C_r2"]; names(x_ALL) <- common_genes
y_ALL <- logfpkm_dev[common_genes,"GSE66582_2C_r2"]; names(y_ALL) <- common_genes
cor_ALL <- round(cor(x_ALL, y_ALL),2)
x_TPG <- x_ALL[TPG]
y_TPG <- y_ALL[TPG]
cor_TPG <- round(cor(x_TPG, y_TPG),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[8], main="", xlab="pu_C_r2 log2(fpkm+1)", ylab="2C_r2 log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_TPG, y_TPG, type="p", pch=16, cex=0.3, col=cccol50[1])
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[8],lwd=2)
lad_reg <- rq(y_TPG ~ x_TPG); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_all=",cor_ALL,sep=""),paste("R_tpg=",cor_TPG,sep="")),bty="n",pch=16, col=cccol[c(8,1)])

dev.off()

pdf("Fig.DevelopmentCorrelation.CombinedSample.pdf",width=8,height=4)
par(mar=c(5,5,5,5),mfrow=c(1,2))

x_ALL <- apply(logfpkm_dev[TPG,c("GSE66582_2C_r1","GSE66582_2C_r2")],1,mean); names(x_ALL) <- TPG
y_ALL <- apply(logfpkm[TPG,c("dox_E_r1","dox_E_r2")],1,mean); names(y_ALL) <- TPG
z_ALL <- apply(logfpkm[TPG,c("dox_C_r1","dox_C_r2")],1,mean); names(z_ALL) <- TPG
cor_y <- round(cor(x_ALL, y_ALL),2)
cor_z <- round(cor(x_ALL, z_ALL),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[2], main="", xlab="2C log2(fpkm+1)", ylab="ESC log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_ALL, z_ALL, type="p", pch=16, cex=0.3, col=cccol50[1])
axis(side=4,at=seq(0,15,5),labels=seq(0,15,5))
mtext("2CLC log2(fpkm+1)", side=4, line=3)
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[2],lwd=2)
lad_reg <- rq(z_ALL ~ x_ALL); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_E=",cor_y,sep=""),paste("R_C=",cor_z,sep="")),bty="n",pch=16, col=cccol[c(2,1)])

x_ALL <- apply(logfpkm_dev[TPG,c("GSE66582_2C_r1","GSE66582_2C_r2")],1,mean); names(x_ALL) <- TPG
y_ALL <- apply(logfpkm[TPG,c("pu_E_r1","pu_E_r2")],1,mean); names(y_ALL) <- TPG
z_ALL <- apply(logfpkm[TPG,c("pu_C_r1","pu_C_r2")],1,mean); names(z_ALL) <- TPG
cor_y <- round(cor(x_ALL, y_ALL),2)
cor_z <- round(cor(x_ALL, z_ALL),2)
plot(x_ALL, y_ALL, type="p", pch=16, cex=0.3, col=cccol50[2], main="", xlab="2C log2(fpkm+1)", ylab="pu ESC log2(fpkm+1)", xlim=c(0,15), ylim=c(0,15))
points(x_ALL, z_ALL, type="p", pch=16, cex=0.3, col=cccol50[1])
axis(side=4,at=seq(0,15,5),labels=seq(0,15,5))
mtext("pu 2CLC log2(fpkm+1)", side=4, line=3)
lad_reg <- rq(y_ALL ~ x_ALL); abline(lad_reg,col=cccol[2],lwd=2)
lad_reg <- rq(z_ALL ~ x_ALL); abline(lad_reg,col=cccol[1],lwd=2)
abline(a=0, b=1, lty=2)
legend("topleft",c(paste("R_E=",cor_y,sep=""),paste("R_C=",cor_z,sep="")),bty="n",pch=16, col=cccol[c(2,1)])
dev.off()
# ############################################################################
# #################                normalise              ####################
# ############################################################################

# # normalise
# library(edgeR)
# experiment_batch <- as.factor(c(rep(1,length(c(dox_C,dox_E,pu_C,pu_E,pup5_E,duxpu_C,duxpu_E))),rep(2,length(c(GSE66582_2C,GSE66582_E2C))),rep(3,length(c(GSE85632_2CLC,GSE85632_ESCs))),rep(4,length(c(GSE94325_ESCs))),rep(5,length(c(Zygote,C2,C4,LC4,C8,C16,C32)))))
# rmbatch_data <- removeBatchEffect(data,batch=experiment_batch)

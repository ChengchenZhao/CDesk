###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
cccol01 <- paste(cccol,"01",sep="")
cccol05 <- paste(cccol,"05",sep="")
cccol30 <- paste(cccol,"30",sep="")
cccol50 <- paste(cccol,"50",sep="")
cccol80 <- paste(cccol,"80",sep="")

options(scipen = 200)

library(amap)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################

############ ESC - 2clc 
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

############ development 
fpkm_dev <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/data/GSE66582_mmEmbryo_RNAseq/GSE66582_GeneSymbol_fpkm.txt",header=T,row.names=1)

all_samples_dev <- colnames(fpkm_dev)
all_genes_dev <- rownames(fpkm_dev)
logfpkm_dev <- log2(fpkm_dev+1)
C2 <- all_samples_dev[grep("^GSE66582_2C([^0-9]|$)",all_samples_dev)]
fpkm_C2 <- apply(fpkm_dev[,C2],1,mean)
names(fpkm_C2) <- all_genes_dev

common_genes <- intersect(all_genes, all_genes_dev)

############ merge data 
fpkm_merge <- cbind(fpkm[common_genes,c(mESC_0h,mESC_p18h,mESC_p30h,mESC_p42h,tCLC_p18h,tCLC_p30h,tCLC_p42h)],fpkm_C2[common_genes])
colnames(fpkm_merge) <- c("mESC","mESC_p18h","mESC_p30h","mESC_p42h","tCLC_p18h","tCLC_p30h","tCLC_p42h","C2")
rownames(fpkm_merge) <- common_genes
logfpkm_merge <- log2(fpkm_merge+1)

timepoints <- c("ESC","18h","30h","42h","2C")

ESC_path <- c("mESC","mESC_p18h","mESC_p30h","mESC_p42h","C2")
C2L_path <- c("mESC","tCLC_p18h","tCLC_p30h","tCLC_p42h","C2")
all_path <- c("mESC","mESC_p18h","mESC_p30h","mESC_p42h","tCLC_p18h","tCLC_p30h","tCLC_p42h","C2")

# DEGs
mESC_0h_vs_mESC_p18h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/mESC_0h_vs_mESC_p18h.deg")[,1])
mESC_p18h_vs_mESC_p30h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/mESC_p18h_vs_mESC_p30h.deg")[,1])
mESC_p30h_vs_mESC_p42h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/mESC_p30h_vs_mESC_p42h.deg")[,1])
mESC_0h_vs_tCLC_p18h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/mESC_0h_vs_tCLC_p18h.deg")[,1])
tCLC_p18h_vs_tCLC_p30h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/tCLC_p18h_vs_tCLC_p30h.deg")[,1])
tCLC_p30h_vs_tCLC_p42h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/tCLC_p30h_vs_tCLC_p42h.deg")[,1])
mESCDEG <- unique(c(mESC_0h_vs_mESC_p18h.deg,mESC_p18h_vs_mESC_p30h.deg,mESC_p30h_vs_mESC_p42h.deg))
tCLCDEG <- unique(c(mESC_0h_vs_tCLC_p18h.deg,tCLC_p18h_vs_tCLC_p30h.deg,tCLC_p30h_vs_tCLC_p42h.deg))

###################################################################################################
################################                PLOT               ################################
###################################################################################################
plot_gene <- intersect(c(mESCDEG, tCLCDEG), common_genes)

clusterData <- logfpkm_merge[plot_gene,all_path]
plotData <- logfpkm_merge[plot_gene,all_path]

k <- 10
set.seed(999)
km <- Kmeans(clusterData,k,method = "correlation",iter.max = 100, nstart = 10)
for (id in seq(k)){
    modGenes = names(which(km$cluster==id))
    fileName = paste("cluster_",id,"_genes.txt",sep="")
    write.table(modGenes,file=fileName,row.names=F,col.names=F,quote=F,sep="\t")
}

selected_cluster <- seq(k)
clusterCol <- cccol[1:k]
pdf(file = "Fig.ExpClustersLinePlot.OnePath.pdf", width = 5.55, height = 3);
for (each_i in seq(length(selected_cluster))){
    each <- selected_cluster[each_i]
    modGenes = names(which(km$cluster==each))
    v1 = apply(plotData[modGenes,],2,mean)
    n <- length(modGenes)
    sd <- apply(plotData[modGenes,],2,sd)
    alpha <- 0.05
    v2 = v1 - sd/sqrt(n)*qt(1-alpha/2,n-1)
    v3 = v1 + sd/sqrt(n)*qt(1-alpha/2,n-1)
    par(mar=c(6.4,4,2,2))
    plot(v1,lwd=3,type="l",col=clusterCol[each_i],ylim=c(min(v2),max(v3)),xlab=NA,ylab="log2(fpkm+1)",xaxt="n",main=paste("C",each,"(n=",length(modGenes),")",sep=""));box(lwd=2)
    axis(side=1,1:length(all_path),all_path,las=2)
    polygon(c(1,1:length(all_path),length(all_path):2),c(v2[1],v3,v2[length(all_path):2]),col=adjustcolor("grey", alpha.f = 0.4),border=NA)
}
dev.off()

pdf(file = "Fig.ExpClustersHeatmap.OnePath.pdf", width = 5.4, height = 2);
plot_data <- plotData
for (each_i in seq(length(selected_cluster))){
    each <- selected_cluster[each_i]
    plotMatrix <- plot_data[names(which(km$cluster==each)),]
    all_exp <- c(as.matrix(plotMatrix))
    zmax <- quantile(all_exp,0.99,na.rm=T)
    zmin <- quantile(all_exp,0.01,na.rm=T)
    ColorRamp <- colorRampPalette(c(cccol[13],"#F2FAB0","red"), bias=1)(10000)   #color list
    ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequencemodGenes <- names(which(km$cluster==each))
    plotMatrix[plotMatrix<zmin] <- zmin
    plotMatrix[plotMatrix>zmax] <- zmax
    par(oma=c(0.5,0.5,0.5,0.5),mar=c(6.4,2,2,2))
    layout(matrix(c(1,1,1,1,1,1,1,1,2),ncol=9,nrow=1,byrow=T))
    image(1:ncol(plotMatrix), 1:nrow(plotMatrix), t(plotMatrix), xaxt="n", yaxt="n", col=ColorRamp, xlab="", ylab="");box(lwd=2)
    axis(side=1,1:ncol(plotMatrix),labels=colnames(plot_data),cex.axis=1.2,las=2)
    image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
    axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
}
dev.off()



# two paths
LastTPESC_path <- c("mESC_p42h","C2")
LastTPC2L_path <- c("tCLC_p42h","C2")

pdf(file = "Fig.ExpClustersLinePlot.TwoPaths.pdf", width = 5.55, height = 3);
for (each_i in seq(length(selected_cluster))){
    each <- selected_cluster[each_i]
    modGenes <- names(which(km$cluster==each))
    alpha <- 0.05
    n <- length(modGenes)
    
    ESC_v1 <- apply(plotData[modGenes,ESC_path[1:4]],2,mean)
    ESC_sd <- apply(plotData[modGenes,ESC_path[1:4]],2,sd)
    ESC_v2 <- ESC_v1 - ESC_sd/sqrt(n)*qt(1-alpha/2,n-1)
    ESC_v3 <- ESC_v1 + ESC_sd/sqrt(n)*qt(1-alpha/2,n-1)
    ESC_x <- c(1:4)
    LastTPESC_v1 <- apply(plotData[modGenes,LastTPESC_path],2,mean)
    LastTPESC_sd <- apply(plotData[modGenes,LastTPESC_path],2,sd)
    LastTPESC_v2 <- LastTPESC_v1 - LastTPESC_sd/sqrt(n)*qt(1-alpha/2,n-1)
    LastTPESC_v3 <- LastTPESC_v1 + LastTPESC_sd/sqrt(n)*qt(1-alpha/2,n-1)
    LastTPESC_x <- c(4,5)

    C2L_v1 <- apply(plotData[modGenes,C2L_path[1:4]],2,mean)
    C2L_sd <- apply(plotData[modGenes,C2L_path[1:4]],2,sd)
    C2L_v2 <- C2L_v1 - C2L_sd/sqrt(n)*qt(1-alpha/2,n-1)
    C2L_v3 <- C2L_v1 + C2L_sd/sqrt(n)*qt(1-alpha/2,n-1)
    C2L_x <- c(1:4)
    LastTPC2L_v1 <- apply(plotData[modGenes,LastTPC2L_path],2,mean)
    LastTPC2L_sd <- apply(plotData[modGenes,LastTPC2L_path],2,sd)
    LastTPC2L_v2 <- LastTPC2L_v1 - LastTPC2L_sd/sqrt(n)*qt(1-alpha/2,n-1)
    LastTPC2L_v3 <- LastTPC2L_v1 + LastTPC2L_sd/sqrt(n)*qt(1-alpha/2,n-1)
    LastTPC2L_x <- c(4,5)

    par(mar=c(6.4,4,2,2))
    plot(ESC_x,ESC_v1,lwd=3,type="l",col=cccol[1],xlim=c(1,length(timepoints)),ylim=c(min(ESC_v2,LastTPESC_v2,C2L_v2,LastTPC2L_v2),max(ESC_v3,LastTPESC_v3,C2L_v3,LastTPC2L_v3)),xlab=NA,ylab="log2(fpkm+1)",xaxt="n",main=paste("C",each,"(n=",length(modGenes),")",sep=""));box(lwd=2)
    axis(side=1,1:length(timepoints),timepoints,las=2)
    polygon(c(ESC_x[1],ESC_x,ESC_x[length(ESC_x):2]),c(ESC_v2[1],ESC_v3,ESC_v2[length(ESC_path[1:4]):2]),col=adjustcolor(cccol[1], alpha.f = 0.4),border=NA)
    points(LastTPESC_x,LastTPESC_v1,lwd=3,type="l",col=cccol[13])
    polygon(c(LastTPESC_x[1],LastTPESC_x,LastTPESC_x[length(LastTPESC_x):2]),c(LastTPESC_v2[1],LastTPESC_v3,LastTPESC_v2[length(LastTPESC_path):2]),col=adjustcolor(cccol[13], alpha.f = 0.4),border=NA)

    points(C2L_x,C2L_v1,lwd=3,type="l",col=cccol[2])
    polygon(c(C2L_x[1],C2L_x,C2L_x[length(C2L_x):2]),c(C2L_v2[1],C2L_v3,C2L_v2[length(C2L_path[1:4]):2]),col=adjustcolor(cccol[2], alpha.f = 0.4),border=NA)
    points(LastTPC2L_x,LastTPC2L_v1,lwd=3,type="l",col=cccol[13])
    polygon(c(LastTPC2L_x[1],LastTPC2L_x,LastTPC2L_x[length(LastTPC2L_x):2]),c(LastTPC2L_v2[1],LastTPC2L_v3,LastTPC2L_v2[length(LastTPC2L_path):2]),col=adjustcolor(cccol[13], alpha.f = 0.4),border=NA)
}
dev.off()

pdf(file = "Fig.ExpClustersHeatmap.TwoPaths.pdf", width = 8, height = 1.8);
plot_data_ESC_all <- c()
plot_data_C2L_all <- c()
tmp_regioncounts <- 0
boundary <- c()

for (each_i in seq(length(selected_cluster))){
    each <- selected_cluster[each_i]
    plot_region <- names(which(km$cluster==each))
    plot_data_ESC <- plotData[plot_region,ESC_path]
    plot_data_C2L <- plotData[plot_region,C2L_path]
    
    plot_data_ESC_all <- rbind(plot_data_ESC_all, plotData[plot_region,ESC_path])
    plot_data_C2L_all <- rbind(plot_data_C2L_all, plotData[plot_region,C2L_path])
    tmp_regioncounts <- tmp_regioncounts+length(plot_region)
    boundary <- c(boundary, tmp_regioncounts)

    layout(matrix(c(rep(1,11),rep(2,11),rep(3,3)),nrow=1,ncol=25))
    all_exp <- c(as.matrix(plot_data_ESC),as.matrix(plot_data_C2L)) # using same scale bar
    zmax <- quantile(all_exp,0.99,na.rm=T)
    zmin <- quantile(all_exp,0.01,na.rm=T)
    
    ColorRamp <- colorRampPalette(c(cccol[13],"#F2FAB0","red"), bias=1)(10000)   #color list
    ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence
    
    tmp_matrix <- plot_data_ESC
    tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
    tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
    tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
    tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
    tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
    tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
    par(mar=c(6,2,2,2))
    image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="ESC",ylab="",useRaster=F);box(lwd=2)
    axis(side=1,1:ncol(tmp_matrix),labels=ESC_path,las=2)
    
    tmp_matrix <- plot_data_C2L
    tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
    tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
    tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
    tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
    tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
    tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
    par(mar=c(6,2,2,2))
    image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="C2L",ylab="",useRaster=F);box(lwd=2)
    axis(side=1,1:ncol(tmp_matrix),labels=C2L_path,las=2)
    
    par(mar=c(6,4,2,2))
    image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="Accessibilty",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
    axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))

}
dev.off()

pdf(file = "Fig.ExpClustersHeatmap.TwoPaths.ALL.pdf", width = 8, height = 6);
layout(matrix(c(rep(1,11),rep(2,11),rep(3,3)),nrow=1,ncol=25))
all_exp <- c(as.matrix(plot_data_ESC_all),as.matrix(plot_data_C2L_all)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c(cccol[13],"#F2FAB0","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_ESC_all
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(6,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="ESC",ylab="",useRaster=F);box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=ESC_path,las=2)
abline(h=boundary-0.5,lwd=2)

tmp_matrix <- plot_data_C2L_all
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(6,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="C2L",ylab="",useRaster=F);box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=C2L_path,las=2)
abline(h=boundary-0.5,lwd=2)

par(mar=c(10,4,20,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="Accessibilty",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))

dev.off()
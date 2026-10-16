###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
red_col <- colorRampPalette(c(cccol[1],cccol[3]), bias=1)(5)[c(1,3,4,5)]
blue_col <- colorRampPalette(c(cccol[1],cccol[2]), bias=1)(5)[c(1,3,4,5)]
grey_col <- colorRampPalette(c(cccol[8],"black"), bias=1)(4)
options(scipen = 200)
library(amap)

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
all_path <- c(mESC_0h,mESC_p18h,mESC_p30h,mESC_p42h,tCLC_p18h,tCLC_p30h,tCLC_p42h)

mESC_0h_vs_mESC_p18h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/mESC_0h_vs_mESC_p18h.deg")[,1])
mESC_p18h_vs_mESC_p30h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/mESC_p18h_vs_mESC_p30h.deg")[,1])
mESC_p30h_vs_mESC_p42h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/mESC_p30h_vs_mESC_p42h.deg")[,1])
mESC_0h_vs_tCLC_p18h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/mESC_0h_vs_tCLC_p18h.deg")[,1])
tCLC_p18h_vs_tCLC_p30h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/tCLC_p18h_vs_tCLC_p30h.deg")[,1])
tCLC_p30h_vs_tCLC_p42h.deg <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/0.mapping/RNAseq/gfold/cutoff_1/tCLC_p30h_vs_tCLC_p42h.deg")[,1])
mESCDEG <- unique(c(mESC_0h_vs_mESC_p18h.deg,mESC_p18h_vs_mESC_p30h.deg,mESC_p30h_vs_mESC_p42h.deg))
tCLCDEG <- unique(c(mESC_0h_vs_tCLC_p18h.deg,tCLC_p18h_vs_tCLC_p30h.deg,tCLC_p30h_vs_tCLC_p42h.deg))

print(length(mESCDEG))
print(length(tCLCDEG))

mESCData <- logfpkm[mESCDEG,mESC_path]
tCLCData <- logfpkm[tCLCDEG,tCLC_path]

###################################################################################################
################################                PLOT               ################################
###################################################################################################

k <- 9
set.seed(999)

km <- Kmeans(mESCData,k,method = "correlation",iter.max =100,nstart=50)
for (id in seq(k)){
	modGenes = names(which(km$cluster==id))
	fileName = paste("cluster_",id,"_gene.mESCData.txt",sep="")
	write.table(cbind(modGenes),file=fileName,row.names=F,col.names=F,quote=F,sep="\t")
}
selected_cluster <- seq(k)
clusterCol <- cccol[1:k]

pdf(file = "Fig.ExpClusters_1.mESCData.pdf", width = 5.55, height = 3);
for (each_i in seq(length(selected_cluster))){
	each <- selected_cluster[each_i]
	modGenes = names(which(km$cluster==each))
	v1 = apply(mESCData[modGenes,],2,mean)
	n <- length(modGenes)
	sd <- apply(mESCData[modGenes,],2,sd)
	alpha <- 0.05
	v2 = v1 - sd/sqrt(n)*qt(1-alpha/2,n-1)
	v3 = v1 + sd/sqrt(n)*qt(1-alpha/2,n-1)
	par(mar=c(6.4,4,2,2))
	plot(v1,lwd=3,type="l",col=clusterCol[each_i],ylim=c(min(v2),max(v3)),xlab=NA,ylab="Log2(FPKM+1)",xaxt="n",main=paste("C",each,"(n=",length(modGenes),")",sep=""));box(lwd=2)
	axis(side=1,1:length(mESC_path),mESC_path,las=2)
	polygon(c(1,1:length(mESC_path),length(mESC_path):2),c(v2[1],v3,v2[length(mESC_path):2]),col=adjustcolor("grey", alpha.f = 0.4),border=NA)
}
dev.off()

pdf(file = "Fig.ExpClusters_2.mESCData.pdf", width = 5.4, height = 2);
plot_data <- mESCData
for (each_i in seq(length(selected_cluster))){
	each <- selected_cluster[each_i]
	plotMatrix <- plot_data[names(which(km$cluster==each)),]
	all_exp <- c(as.matrix(plotMatrix))
	zmax <- quantile(all_exp,0.99,na.rm=T)
	zmin <- quantile(all_exp,0.01,na.rm=T)
	ColorRamp <- colorRampPalette(c("lightblue","red"), bias=1)(10000)   #color list
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

km <- Kmeans(tCLCData,k,method = "correlation",iter.max =100,nstart=50)
for (id in seq(k)){
	modGenes = names(which(km$cluster==id))
	fileName = paste("cluster_",id,"_gene.tCLCData.txt",sep="")
	write.table(cbind(modGenes),file=fileName,row.names=F,col.names=F,quote=F,sep="\t")
}
selected_cluster <- seq(k)
clusterCol <- cccol[1:k]

pdf(file = "Fig.ExpClusters_1.tCLCData.pdf", width = 5.55, height = 3);
for (each_i in seq(length(selected_cluster))){
	each <- selected_cluster[each_i]
	modGenes = names(which(km$cluster==each))
	v1 = apply(tCLCData[modGenes,],2,mean)
	n <- length(modGenes)
	sd <- apply(tCLCData[modGenes,],2,sd)
	alpha <- 0.05
	v2 = v1 - sd/sqrt(n)*qt(1-alpha/2,n-1)
	v3 = v1 + sd/sqrt(n)*qt(1-alpha/2,n-1)
	par(mar=c(6.4,4,2,2))
	plot(v1,lwd=3,type="l",col=clusterCol[each_i],ylim=c(min(v2),max(v3)),xlab=NA,ylab="Log2(FPKM+1)",xaxt="n",main=paste("C",each,"(n=",length(modGenes),")",sep=""));box(lwd=2)
	axis(side=1,1:length(tCLC_path),tCLC_path,las=2)
	polygon(c(1,1:length(tCLC_path),length(tCLC_path):2),c(v2[1],v3,v2[length(tCLC_path):2]),col=adjustcolor("grey", alpha.f = 0.4),border=NA)
}
dev.off()

pdf(file = "Fig.ExpClusters_2.tCLCData.pdf", width = 5.4, height = 2);
plot_data <- tCLCData
for (each_i in seq(length(selected_cluster))){
	each <- selected_cluster[each_i]
	plotMatrix <- plot_data[names(which(km$cluster==each)),]
	all_exp <- c(as.matrix(plotMatrix))
	zmax <- quantile(all_exp,0.99,na.rm=T)
	zmin <- quantile(all_exp,0.01,na.rm=T)
	ColorRamp <- colorRampPalette(c("lightblue","red"), bias=1)(10000)   #color list
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

# selected_cluster <- c(20,17,1,9,2,8,5)
# clusterCol <- c("#CE0013","#FA8072","#32CD32","#7FFFD4","#3A5FCD","#004138","#00CED1")

# pdf(file = "Fig.ExpClusters.pdf", width = 5.55, height = 3);
# for (each_i in seq(length(selected_cluster))){
# 	each <- selected_cluster[each_i]
# 	modGenes = names(which(km$cluster==each))
# 	v1 = apply(mESCData[modGenes,],2,mean)
# 	n <- length(modGenes)
# 	sd <- apply(mESCData[modGenes,],2,sd)
# 	alpha <- 0.05
# 	v2 = v1 - sd/sqrt(n)*qt(1-alpha/2,n-1)
# 	v3 = v1 + sd/sqrt(n)*qt(1-alpha/2,n-1)
# 	# plot(v1,lwd=3,type="l",col=clusterCol[each],xlab=NA,ylab="Expression level(Log FPKM)",ylim=c(0,4),xaxt="n")
# 	plot(v1,lwd=3,type="l",col=clusterCol[each_i],xlab=NA,ylab="Log2(FPKM+1)",xaxt="n")
# 	axis(side=1,1:length(n_time_point),rep("",11),las=2);axis(side=2);box()
# 	polygon(c(1,1:length(n_time_point),length(n_time_point):2),c(v2[1],v3,v2[length(n_time_point):2]),col=adjustcolor("grey", alpha.f = 0.4),border=NA)
# 	# fileName = paste(k,"cluster_",each,"_gene.txt",sep="")
# 	# write.table(cbind(modGenes),file=fileName,row.names=F,col.names=F,quote=F,sep="\t")
# }
# dev.off()

# pdf(file = "SFig2D_2.pdf", width = 5.4, height = 2);
# plot_data <- mESCData
# for (each_i in seq(length(selected_cluster))){
# 	each <- selected_cluster[each_i]
# 	plotMatrix <- plot_data[names(which(km$cluster==each)),]
# 	all_exp <- c(as.matrix(plotMatrix))
# 	zmax <- quantile(all_exp,0.99)
# 	zmin <- quantile(all_exp,0.01)
# 	ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
# 	ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequencemodGenes <- names(which(km$cluster==each))
# 	plotMatrix[plotMatrix<zmin] <- zmin
# 	plotMatrix[plotMatrix>zmax] <- zmax
# 	par(oma=c(0.5,0.5,0.5,0.5),mar=c(8,2,2,2))
# 	layout(matrix(c(1,1,1,1,1,1,1,1,2),ncol=9,nrow=1,byrow=T))
# 	image(1:ncol(plotMatrix), 1:nrow(plotMatrix), t(plotMatrix), xaxt="n", yaxt="n", col=ColorRamp, xlab="", ylab="")
# 	axis(side=1,1:ncol(plotMatrix),labels=colnames(plot_data),cex.axis=1.2,las=2);box()
# 	image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T)
# 	axis(side=2,c(zmin,round((zmax-zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax-zmin)/2,1),round(zmax,1)))
# }
# dev.off()


# pdf(file = "SFig2E.pdf", width = 10, height = 4.5)
# for (each_i in seq(length(selected_cluster))){
# 	each <- selected_cluster[each_i]
# 	file <- read.table(paste("../data/",each,"_GO.txt",sep=''),sep='\t')
# 	file[,3] <- -log10(file[,3])
# 	file <- file[order(file[,3]),]
# 	file <- file[!duplicated(file),]
# 	par(mar=c(8,30,8,5))
# 	# barplot(file[,3],horiz=T,names=file[,2],las=2,border=F,col=clusterCol[each_i],space=0.5,axes=F)
# 	barplot(file[,3],horiz=T,names=file[,2],las=2,border=F,col="grey",space=0.5,axes=F)
# 	axis(side=1);mtext("-log10(p-value)",side=1,line=2)
# }
# dev.off()

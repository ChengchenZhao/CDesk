############################################################################
########################	   function define       #######################
############################################################################
# library(amap)
# cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol <- c("#C7A609","#A40126","#235E7C")
ErrorBar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
HalfErrorBar <- function(x, y, upper, lower=upper, length=0.05,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	for (i in seq(length(x))){
		if (y[i]>0){arrows(x[i],y[i]+upper[i], x[i], y[i], angle=90, code=1, length=length, ...)}
		else{arrows(x[i],y[i],x[i],y[i]-lower[i], angle=90, code=2, length=length, ...)}
	}
}

ExpressionPatternBarplot <- function(genes,cols,data_path,data_mean,data_sd,data_n,data_label,data_mark){
	for(each in genes){
		plot_matrix <- data_mean[each,data_path]
		plot_sd <- data_sd[each,data_path]
		n <- data_n
		ymax <- max(plot_matrix+qt(0.975, n-1) * plot_sd / sqrt(n))
		bp <- barplot(plot_matrix,col=cols,main=each,ylab="log2(fpkm+1)",names.arg=data_label,las=2,bty="l",border = "white",ylim=c(0,ymax*1.3))
		HalfErrorBar(bp,plot_matrix, qt(0.975, n-1) * plot_sd / sqrt(n),col="black",length=0.05,lwd=1.5)
		dist <- ymax/20
		by <- plot_matrix+qt(0.975, n-1) * plot_sd / sqrt(n)
		if(!is.na(data_mark[each]) & data_mark[each]!="?"){
			segments(bp[1], by[1]+dist/2, bp[1], ymax*1.2, lwd=1.5)
			segments(bp[6], by[6]+dist/2, bp[6], ymax*1.2, lwd=1.5)
			segments(bp[1], ymax*1.2, ((bp[1]+bp[6])/2)-0.5, ymax*1.2, lwd=1.5)
			segments(bp[6], ymax*1.2, ((bp[1]+bp[6])/2)+0.5, ymax*1.2, lwd=1.5)
			text(x=(bp[1]+bp[6])/2, y=ymax*1.2, label=data_mark[each], col="red")
		}
	}
}

UnExpressedGene <- function(x,cutoff){
	# the unexpressed gene = (exp of all samples < cutoff)
	return(sum(x >= cutoff) == 0)
}
HighExpressedGene <- function(x,cutoff){
	# the high expressed gene = (exp of all samples >= cutoff)
	return(sum(x >= cutoff) == length(x))
}

ExpressionPatternLines <- function(x,plot_matrix,COL){
	v1 <- apply(plot_matrix,2,mean,na.rm=T)
	v2 <- apply(plot_matrix,2,quantile,probs=c(0.25,0.75))[1,]
	v3 <- apply(plot_matrix,2,quantile,probs=c(0.25,0.75))[2,]
	polygon(c(1,x,rev(x)[1:length(x)-1]),c(v2[1],v3,v2[ncol(plot_matrix):2]),col=adjustcolor(COL, alpha.f = 0.2),border=NA)
	points(x,v1,lwd=2,type="l",col=COL,xaxt="n")
}

Boxplot_List <- function(plot_list,MAIN,NAMES,YLAB,COL){
	n <- length(plot_list)
	xpos <- 0:(n-1)+1.5
	p_value <- c()
	for (i in seq(n-1)){
		tryCatch({p_value <- c(p_value,t.test(plot_list[[i]],plot_list[[i+1]])$p.value)},
			error = function(err){p_value <- c(p_value,1)})
	}
	mark <- symnum(p_value, cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","-"))
	bp <- boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),outline=F,plot=F)
	ylim <- range(bp$stats,na.rm=T)
	dist <- (ylim[2]-ylim[1])/20
	ylim[2] <- ylim[2]+5*dist
	boxplot(plot_list,at=0:(n-1)+1,boxwex=0.5, xlim=c(0.5,n+0.5),ylim=ylim,col=paste(COL,"50",sep=""),border=COL,outline=F,names=NAMES,main=MAIN,lwd=2,ylab=YLAB,las=2)
	box(lwd=2)
	ypos_1 <- bp$stats[5,][1:n-1]
	ypos_2 <- bp$stats[5,][2:n]
	# legend("topleft",c("***:p<0.001","**:p<0.01","*:p<0.05"),text.col="red",bty="n")
	if (length(mark)!=0){
		for(i in 1:length(mark)){
			if(!is.na(mark[i]) & mark[i]!="?"){
				segments(xpos[i]-.45, ypos_1[i]+dist/2, xpos[i]-.45, max(ypos_1[i], ypos_2[i])+dist)
				segments(xpos[i]+.45, ypos_2[i]+dist/2, xpos[i]+.45, max(ypos_1[i], ypos_2[i])+dist)
				segments(xpos[i]-.45, max(ypos_1[i], ypos_2[i])+dist, xpos[i]-0.3, max(ypos_1[i], ypos_2[i])+dist)
				segments(xpos[i]+.45, max(ypos_1[i], ypos_2[i])+dist, xpos[i]+0.3, max(ypos_1[i], ypos_2[i])+dist)
				text(x=xpos[i], y=max(ypos_1[i], ypos_2[i])+dist, label=mark[i], col="black")
			}
		}
	}
	for (i in seq(n)){
		text(x=i, y=bp$stats[3,i], label=length(na.omit(plot_list[[i]])), col="white",cex=0.3)
	}
}
############################################################################
#########################	   read in data       ##########################
############################################################################

fpkm <- read.table("../../Data/inhouse_scRNAseq_FPKM.txt",header=T,row.names=1)
all_samples <- c("Ctrl_1C14hpa_r1","Ctrl_1C14hpa_r2","Ctrl_1C14hpa_r3","Ctrl_1C6hpa_r1","Ctrl_1C6hpa_r2","Ctrl_1C6hpa_r3","Ctrl_4C_r1","Ctrl_4C_r2","Ctrl_4C_r3","Ctrl_8C_r1","Ctrl_8C_r2","Ctrl_8C_r3","Ctrl_CC_r1","Ctrl_CC_r2","Ctrl_CC_r3","Ctrl_E2C_r1","Ctrl_E2C_r2","Ctrl_E2C_r3","Ctrl_ICM_r1","Ctrl_ICM_r2","Ctrl_ICM_r3","Ctrl_L2C_r1","Ctrl_L2C_r2","Ctrl_L2C_r3","Ctrl_Morula_r1","Ctrl_Morula_r2","Ctrl_Morula_r3","Ctrl_TE_r1","Ctrl_TE_r2","Ctrl_TE_r3","M14h_1C14hpa_r1","M14h_1C14hpa_r2","M14h_4C_r1","M14h_4C_r2","M14h_CC_r1","M14h_CC_r2","M14h_E2C_r1","M14h_E2C_r2","M14h_L2C_r1","M14h_L2C_r2")
colnames(fpkm) <- all_samples
all_genes <- row.names(fpkm)
logfpkm <- log2(fpkm+1)

Ct_times <- c("CC","1C6hpa","1C14hpa","E2C","L2C","4C","8C","Morula","ICM","TE")
Tt_times <- c("CC","1C14hpa","E2C","L2C","4C")

Ctrl_CC <- all_samples[grep("Ctrl_CC",all_samples)]
Ctrl_1C6hpa <- all_samples[grep("Ctrl_1C6hpa",all_samples)]
Ctrl_1C14hpa <- all_samples[grep("Ctrl_1C14hpa",all_samples)]
Ctrl_E2C <- all_samples[grep("Ctrl_E2C",all_samples)]
Ctrl_L2C <- all_samples[grep("Ctrl_L2C",all_samples)]
Ctrl_4C <- all_samples[grep("Ctrl_4C",all_samples)]
Ctrl_8C <- all_samples[grep("Ctrl_8C",all_samples)]
Ctrl_Morula <- all_samples[grep("Ctrl_Morula",all_samples)]
Ctrl_ICM <- all_samples[grep("Ctrl_ICM",all_samples)]
Ctrl_TE <- all_samples[grep("Ctrl_TE",all_samples)]
M14h_CC <- all_samples[grep("M14h_CC",all_samples)]
M14h_1C14hpa <- all_samples[grep("M14h_1C14hpa",all_samples)]
M14h_E2C <- all_samples[grep("M14h_E2C",all_samples)]
M14h_L2C <- all_samples[grep("M14h_L2C",all_samples)]
M14h_4C <- all_samples[grep("M14h_4C",all_samples)]

# Ctrl_avglogfpkm <- log2(1+cbind(apply(fpkm[,Ctrl_CC],1,mean),apply(fpkm[,Ctrl_1C6hpa],1,mean),apply(fpkm[,Ctrl_1C14hpa],1,mean),apply(fpkm[,Ctrl_E2C],1,mean),apply(fpkm[,Ctrl_L2C],1,mean),apply(fpkm[,Ctrl_4C],1,mean),apply(fpkm[,Ctrl_8C],1,mean),apply(fpkm[,Ctrl_Morula],1,mean),apply(fpkm[,Ctrl_ICM],1,mean),apply(fpkm[,Ctrl_TE],1,mean)))
Ctrl_avglogfpkm <- cbind(apply(logfpkm[,Ctrl_CC],1,mean),apply(logfpkm[,Ctrl_1C6hpa],1,mean),apply(logfpkm[,Ctrl_1C14hpa],1,mean),apply(logfpkm[,Ctrl_E2C],1,mean),apply(logfpkm[,Ctrl_L2C],1,mean),apply(logfpkm[,Ctrl_4C],1,mean),apply(logfpkm[,Ctrl_8C],1,mean),apply(logfpkm[,Ctrl_Morula],1,mean),apply(logfpkm[,Ctrl_ICM],1,mean),apply(logfpkm[,Ctrl_TE],1,mean))
colnames(Ctrl_avglogfpkm) <- c("Ctrl_CC","Ctrl_1C6hpa","Ctrl_1C14hpa","Ctrl_E2C","Ctrl_L2C","Ctrl_4C","Ctrl_8C","Ctrl_Morula","Ctrl_ICM","Ctrl_TE")
# M14h_avglogfpkm <- log2(1+cbind(apply(fpkm[,M14h_CC],1,mean),apply(fpkm[,M14h_1C14hpa],1,mean),apply(fpkm[,M14h_E2C],1,mean),apply(fpkm[,M14h_L2C],1,mean),apply(fpkm[,M14h_4C],1,mean)))
M14h_avglogfpkm <- cbind(apply(logfpkm[,M14h_CC],1,mean),apply(logfpkm[,M14h_1C14hpa],1,mean),apply(logfpkm[,M14h_E2C],1,mean),apply(logfpkm[,M14h_L2C],1,mean),apply(logfpkm[,M14h_4C],1,mean))
colnames(M14h_avglogfpkm) <- c("M14h_CC","M14h_1C14hpa","M14h_E2C","M14h_L2C","M14h_4C")

Dev_fpkm <- read.table("../../Data/GSE71434_RNAseq_mmEmbryo_FPKM.txt",header=T,row.names=1)
Dev_logfpkm <- log2(Dev_fpkm+1)
# Dev_avglogfpkm <- log2(1+cbind(apply(Dev_fpkm[,c("MII_r1","MII_r2")],1,mean),apply(Dev_fpkm[,c("zygote_r1","zygote_r2")],1,mean),apply(Dev_fpkm[,c("EC2_r1","EC2_r2")],1,mean),apply(Dev_fpkm[,c("LC2_r1","LC2_r2")],1,mean),apply(Dev_fpkm[,c("C4_r1","C4_r2")],1,mean),apply(Dev_fpkm[,c("C8_r1","C8_r2")],1,mean),apply(Dev_fpkm[,c("ICM_r1","ICM_r2")],1,mean)))
Dev_avglogfpkm <- cbind(apply(Dev_logfpkm[,c("MII_r1","MII_r2")],1,mean),apply(Dev_logfpkm[,c("zygote_r1","zygote_r2")],1,mean),apply(Dev_logfpkm[,c("EC2_r1","EC2_r2")],1,mean),apply(Dev_logfpkm[,c("LC2_r1","LC2_r2")],1,mean),apply(Dev_logfpkm[,c("C4_r1","C4_r2")],1,mean),apply(Dev_logfpkm[,c("C8_r1","C8_r2")],1,mean),apply(Dev_logfpkm[,c("ICM_r1","ICM_r2")],1,mean))
colnames(Dev_avglogfpkm) <- c("Dev_MII","Dev_zygote","Dev_E2C","Dev_L2C","Dev_4C","Dev_8C","Dev_ICM")

# Major ZGA_genes 
# ("zygote-E2C DE-up" OR "E2C-L2C DE-up") AND "MII unexpressed"

# MII_unexpressed <- rownames(Dev_fpkm)[which(apply(Dev_fpkm[,c("MII_r1","MII_r2")],1,UnExpressedGene,1))]
# zygote_EC2_up <- read.table("/mnt/Storage/home/zhaochengchen/Work/3.MCNT/data/GSE71434_RNAseq_mmEmbryo/GFOLD/cutoff_0.58/zygote_EC2_up.geneSymbol")[,1]
# EC2_LC2_up <- read.table("/mnt/Storage/home/zhaochengchen/Work/3.MCNT/data/GSE71434_RNAseq_mmEmbryo/GFOLD/cutoff_0.58/EC2_LC2_up.geneSymbol")[,1]
# Major_1 <- intersect(zygote_EC2_up,MII_unexpressed)
# Major_2 <- intersect(EC2_LC2_up,MII_unexpressed)
# Major <- intersect(union(zygote_EC2_up,EC2_LC2_up),MII_unexpressed)

# # Minor ZGA genes
# # "MII-zygote DE-up" AND "MII unexpressed"
# MII_zygote_up <- read.table("/mnt/Storage/home/zhaochengchen/Work/3.MCNT/data/GSE71434_RNAseq_mmEmbryo/GFOLD/cutoff_0.58/MII_zygote_up.geneSymbol")[,1]
# Minor <- intersect(MII_zygote_up,MII_unexpressed)

# # Oocyte expressed genes
# # "MII highly-expressed" AND "MII-zygote DE-down"
# MII_highlyexpressed <- rownames(Dev_fpkm)[which(apply(Dev_fpkm[,c("MII_r1","MII_r2")],1,HighExpressedGene,3))]
# MII_zygote_dn <- read.table("/mnt/Storage/home/zhaochengchen/Work/3.MCNT/data/GSE71434_RNAseq_mmEmbryo/GFOLD/cutoff_0.58/MII_zygote_dn.geneSymbol")[,1]
# MII_specexpressed <- intersect(MII_highlyexpressed, MII_zygote_dn)

# write.table(Major_1,file="Major_1.ZGAGenes.txt",quote=F,sep="\t",col.names=F,row.names=F)
# write.table(Major_2,file="Major_2.ZGAGenes.txt",quote=F,sep="\t",col.names=F,row.names=F)
# write.table(Minor,file="Minor.ZGAGenes.txt",quote=F,sep="\t",col.names=F,row.names=F)

# write.table(MII_highlyexpressed,file="MII.HighlyExpressedGenes.txt",quote=F,sep="\t",col.names=F,row.names=F)
# write.table(MII_specexpressed,file="MII.SpecExpressedGenes.txt",quote=F,sep="\t",col.names=F,row.names=F)

Major_1 <- as.vector(read.table("Major_1.ZGAGenes.txt")[,1])
Major_2 <- as.vector(read.table("Major_2.ZGAGenes.txt")[,1])
Major <- union(Major_1,Major_2)
write.table(file="Major.ZGAGenes.txt",Major,col.names=F,row.names=F,quote=F)

M14h_genes <- as.vector(read.table("MajorinL2C.M14h_genes.txt")[,1])

###################################################################################################
##################################            plot            #####################################
###################################################################################################

# Selected path
Ct_path <- c("Ctrl_1C14hpa","Ctrl_E2C","Ctrl_L2C","Ctrl_4C")
Tt_path <- c("M14h_1C14hpa","M14h_E2C","M14h_L2C","M14h_4C")
De_path <- c("Dev_zygote","Dev_E2C","Dev_L2C","Dev_4C")

pdf("Fig4d.ExpressionChangePatternSelected.Lineplot.Major.pdf", width = 5, height = 4);
par(mar=c(4,4,4,2))
plot(1,type="n",xaxt="n",ylim=c(0,5.5),xlim=c(1,4),xlab="",ylab="Log2(FPKM+1)",main="Major",las=2);box(lwd=2)
axis(side=1,at=seq(4),labels=c("1C","E2C","L2C","4C"),las=2)
plot_matrix <- Ctrl_avglogfpkm[Major,Ct_path]; ExpressionPatternLines(seq(4),plot_matrix,cccol[1])
plot_matrix <- M14h_avglogfpkm[Major,Tt_path]; ExpressionPatternLines(seq(4),plot_matrix,cccol[2])
plot_matrix <- Dev_avglogfpkm[Major,De_path]; ExpressionPatternLines(seq(4),plot_matrix,cccol[3])
legend("topleft",c("Ctrl","M14h","Dev"),col=cccol,bty="n",pch=15,ncol=3)
dev.off()


plot_matrix_1 <- Ctrl_avglogfpkm[Major,Ct_path]
plot_matrix_2 <- M14h_avglogfpkm[Major,Tt_path]
plot_matrix_3 <- Dev_avglogfpkm[Major,De_path]
pdf("Fig4e.ExpressionChangePatternSelected.Heatmap.Major.pdf", width = 6, height = 8);
par(oma=c(1,1,1,1),mar=c(6,8,2,1),mfrow=c(1,3))
tmp_vector <- c(plot_matrix_1,plot_matrix_2,plot_matrix_3)
zmax <- as.numeric(quantile(tmp_vector,0.99,na.rm=T))
zmin <- as.numeric(quantile(tmp_vector,0.01,na.rm=T))
ColorRamp <- colorRampPalette(c("lightblue","white","#F2FAB0","red","#CE0013"), bias=1)(10000)
ColorLevels <- seq(to=zmax,from=zmin, length=10000)

tmp_matrix <- plot_matrix_1
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), xaxt="n", yaxt="n", col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", ylab="",main="Ctrl",cex.main=1);box(lwd=2)
axis(side=1,at=seq(length(Ct_path)),labels=Ct_path,cex.axis=1,las=2)
axis(side=2,at=seq(nrow(tmp_matrix)),labels=rownames(tmp_matrix),cex.axis=1,las=2)

tmp_matrix <- plot_matrix_2
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), xaxt="n", yaxt="n", col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", ylab="",main="M14h",cex.main=1);box(lwd=2)
axis(side=1,at=seq(length(Tt_path)),labels=Tt_path,cex.axis=1,las=2)
# axis(side=2,at=seq(nrow(tmp_matrix)),labels=rownames(tmp_matrix),cex.axis=1,las=2)

tmp_matrix <- plot_matrix_3
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), xaxt="n", yaxt="n", col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", ylab="",main="Dev",cex.main=1);box(lwd=2)
axis(side=1,at=seq(length(De_path)),labels=De_path,cex.axis=1,las=2)
# axis(side=2,at=seq(nrow(tmp_matrix)),labels=rownames(tmp_matrix),cex.axis=1,las=2)
dev.off()

pdf("Fig4e.ExpressionChangePatternSelected.Heatmap.Major.Legend.pdf", width = 3, height = 3);
par(oma=c(1,1,1,1),mar=c(6,2,6,2))
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, ylab="",xlab="Log2(FPKM+1)",xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=1,c(zmin,(zmin+zmax)/2,zmax),labels=c(round(zmin,2),round((zmin+zmax)/2,2),round(zmax,2)),cex.axis=1)
dev.off()




# pdf("ExpressionChangePattern.Lineplot.M14h_genes.pdf", width = 3.8, height = 3.8);
# par(mar=c(4,4,8,2))
# plot(1,type="n",xaxt="n",ylim=c(0,7),xlim=c(0,10),xlab="",ylab="Log2(FPKM+1)",main="M14h_genes");box(lwd=2)
# axis(side=1,at=seq(10),labels=c("CC","1C6hpa","1C14hpa","E2C","L2C","4C","8C","Morula","ICM","TE"),las=2)
# axis(side=3,at=c(0,3,4,5,6,7,9),labels=c("MII","zygote","E2C","L2C","4C","8C","ICM"),las=2)
# plot_matrix <- Ctrl_avglogfpkm[M14h_genes,]; ExpressionPatternLines(seq(10),plot_matrix,cccol[1])
# plot_matrix <- M14h_avglogfpkm[M14h_genes,]; ExpressionPatternLines(c(1,3,4,5,6),plot_matrix,cccol[2])
# plot_matrix <- Dev_avglogfpkm[M14h_genes,]; ExpressionPatternLines(c(0,3,4,5,6,7,9),plot_matrix,cccol[3])
# legend("topleft",c("Ctrl","M14h","Dev"),col=cccol,bty="n",pch=15)
# dev.off()

# Selected path
Ct_path <- c("Ctrl_CC","Ctrl_1C14hpa","Ctrl_E2C","Ctrl_L2C","Ctrl_4C")
Tt_path <- c("M14h_CC","M14h_1C14hpa","M14h_E2C","M14h_L2C","M14h_4C")
De_path <- c("Dev_MII","Dev_zygote","Dev_E2C","Dev_L2C","Dev_4C")

pdf("ExpressionChangePatternSelected.Lineplot.M14h_genes.pdf", width = 3.8, height = 3.8);
par(mar=c(4,4,8,2))
plot(1,type="n",xaxt="n",ylim=c(0,5),xlim=c(0,5),xlab="",ylab="Log2(FPKM+1)",main="M14h_genes");box(lwd=2)
axis(side=1,at=seq(5),labels=c("CC","1C","E2C","L2C","4C"),las=2)
axis(side=3,at=c(0,2,3,4,5),labels=c("MII","zygote","E2C","L2C","4C"),las=2)
plot_matrix <- Ctrl_avglogfpkm[M14h_genes,Ct_path]; ExpressionPatternLines(seq(5),plot_matrix,cccol[1])
plot_matrix <- M14h_avglogfpkm[M14h_genes,Tt_path]; ExpressionPatternLines(seq(5),plot_matrix,cccol[2])
plot_matrix <- Dev_avglogfpkm[M14h_genes,De_path]; ExpressionPatternLines(c(0,2,3,4,5),plot_matrix,cccol[3])
legend("topleft",c("Ctrl","M14h","Dev"),col=cccol,bty="n",pch=15)
dev.off()

Ct_path <- c("Ctrl_1C14hpa","Ctrl_E2C","Ctrl_L2C","Ctrl_4C")
Tt_path <- c("M14h_1C14hpa","M14h_E2C","M14h_L2C","M14h_4C")
De_path <- c("Dev_zygote","Dev_E2C","Dev_L2C","Dev_4C")
plot_matrix_1 <- Ctrl_avglogfpkm[M14h_genes,Ct_path]
plot_matrix_2 <- M14h_avglogfpkm[M14h_genes,Tt_path]
plot_matrix_3 <- Dev_avglogfpkm[M14h_genes,De_path]
pdf("ExpressionChangePatternSelected.Heatmap.M14h_genes.pdf", width = 5.2, height = 4);
par(oma=c(1,1,1,1),mar=c(6,8,2,1),mfrow=c(1,3))
tmp_vector <- c(plot_matrix_1,plot_matrix_2,plot_matrix_3)
zmax <- as.numeric(quantile(tmp_vector,0.99,na.rm=T))
zmin <- as.numeric(quantile(tmp_vector,0.01,na.rm=T))
ColorRamp <- colorRampPalette(c("lightblue","white","#F2FAB0","red","#CE0013"), bias=1)(10000)
ColorLevels <- seq(to=zmax,from=zmin, length=10000)

tmp_matrix <- plot_matrix_1
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), xaxt="n", yaxt="n", col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", ylab="",main="Ctrl",cex.main=1);box(lwd=2)
axis(side=1,at=seq(length(Ct_path)),labels=Ct_path,cex.axis=1,las=2)
axis(side=2,at=seq(nrow(tmp_matrix)),labels=rownames(tmp_matrix),cex.axis=1,las=2)

tmp_matrix <- plot_matrix_2
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), xaxt="n", yaxt="n", col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", ylab="",main="M14h",cex.main=1);box(lwd=2)
axis(side=1,at=seq(length(Tt_path)),labels=Tt_path,cex.axis=1,las=2)
# axis(side=2,at=seq(nrow(tmp_matrix)),labels=rownames(tmp_matrix),cex.axis=1,las=2)

tmp_matrix <- plot_matrix_3
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), xaxt="n", yaxt="n", col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", ylab="",main="Dev",cex.main=1);box(lwd=2)
axis(side=1,at=seq(length(De_path)),labels=De_path,cex.axis=1,las=2)
# axis(side=2,at=seq(nrow(tmp_matrix)),labels=rownames(tmp_matrix),cex.axis=1,las=2)
dev.off()

pdf("ExpressionChangePatternSelected.Heatmap.M14h_genes.Legend.pdf", width = 3, height = 3);
par(oma=c(1,1,1,1),mar=c(6,2,6,2))
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, ylab="",xlab="Log2(FPKM+1)",xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=1,c(zmin,(zmin+zmax)/2,zmax),labels=c(round(zmin,2),round((zmin+zmax)/2,2),round(zmax,2)),cex.axis=1)
dev.off()



###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
library(scales)

sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")

NormalizeCounts <- function(indata){
  return(indata*1e6/sum(indata) + 1)
}

MarkerGenesExpressionPlot_1D <- function(genes,data,x,y,cols,name){
	pdf(paste(name,"_","MarkerGenesExpressionPlot_1D.pdf",sep=""),width=5,height=4)
	plot(x,y,xlab="Dimension 1",ylab="Dimension 2",pch=20,cex=0,type="n",yaxt="n")
	axis(side=2,at=seq(7),labels=rev(c("d0","d1","d3","d6","d9","d12","F10")),las=2)
	box(lwd=2)
	for (i in seq(length(genes))){
		l <- data[genes[i],]
		initCol <- c("white",cols[i])
		mi <- min(l,na.rm=TRUE)
		ma <- max(l,na.rm=TRUE)
		v <- round((l - mi)/(ma - mi),2)
		for ( k in 1:length(v) ){
		    points(x[k],y[k],col=alpha(cols[i],as.numeric(v[k])),pch=16)
		}
	}
	legend("topright",bty="n",c(genes),col=cols,pch=16)
	dev.off()
	
	pdf(paste(name,"_","CellTypeDecision_legend.pdf",sep=""),width=3,height=4)
	par(mfrow=c(3,1))
	for (i in seq(length(genes))){
		l <- data[genes[i],]
		initCol <- c("white",cols[i])
		mi <- min(l,na.rm=TRUE)
		ma <- max(l,na.rm=TRUE)
		ColorRamp <- colorRampPalette(initCol)(100)
		ColorLevels <- seq(mi, ma, length=length(ColorRamp))
		image(ColorLevels,1,
			matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),
			col=ColorRamp,
			xlab=paste(genes[i]," Expression Log2(fpkm+1)",sep=""),ylab="",
			yaxt="n")
		box(lwd=2)
	}
	dev.off()
}

MarkerGenesExpressionPlot_2D <- function(genes,data,x,y,cols,name){
	pdf(paste(name,"_","MarkerGenesExpressionPlot_2D.pdf",sep=""),width=5,height=5)
	plot(x,y,xlab="Dimension 1",ylab="Dimension 2",pch=20,cex=0,type="n")
	box(lwd=2)
	for (i in seq(length(genes))){
		l <- data[genes[i],]
		initCol <- c("white",cols[i])
		mi <- min(l,na.rm=TRUE)
		ma <- max(l,na.rm=TRUE)
		v <- round((l - mi)/(ma - mi),2)
		for ( k in 1:length(v) ){
		    points(x[k],y[k],col=alpha(cols[i],as.numeric(v[k])),pch=16)
		}
	}
	legend("topright",bty="n",c(genes),col=cols,pch=16)
	dev.off()
	
	pdf(paste(name,"_","CellTypeDecision_legend.pdf",sep=""),width=3,height=4)
	par(mfrow=c(3,1))
	for (i in seq(length(genes))){
		l <- data[genes[i],]
		initCol <- c("white",cols[i])
		mi <- min(l,na.rm=TRUE)
		ma <- max(l,na.rm=TRUE)
		ColorRamp <- colorRampPalette(initCol)(100)
		ColorLevels <- seq(mi, ma, length=length(ColorRamp))
		image(ColorLevels,1,
			matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),
			col=ColorRamp,
			xlab=paste(genes[i]," Expression Log2(fpkm+1)",sep=""),ylab="",
			yaxt="n")
		box(lwd=2)
	}
	dev.off()
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################

##### expression matrix
raw_data <- read.table(sys_argv[1],row.names=1,header=T)
normalized_counts <- apply(raw_data,2,NormalizeCounts)
logfpkm <- log2(normalized_counts)
cluster_result <- read.table(sys_argv[2],row.names=1,header=T)
marker_gene <- sys_argv[3]
x_1 <- cluster_result[,1]
y_1 <- cluster_result[,2]
###################################################################################################
################################                PLOT               ################################
###################################################################################################
for (each in seq(length(marker_gene))){
	MarkerGenesExpressionPlot_2D(marker_gene[each],logfpkm,x_1,y_1,cccol[each],marker_gene[each])
}

# MarkerGenesExpressionPlot_1D(c("CBX5","DNMT3L"),logfpkm,x_1,y_1,cccol[1:2],"CBX5_DNMT3L")
# MarkerGenesExpressionPlot_2D(c("CBX5","DNMT3L"),logfpkm,x_2,y_2,cccol[1:2],"CBX5_DNMT3L")

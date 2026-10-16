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

Boxplot_List <- function(plot_list,MAIN,NAMES,YLAB,COLS,COLS50){
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
	ylim[2] <- ylim[2]+1*dist
	boxplot(plot_list,at=0:(n-1)+1,boxwex=0.6, xlim=c(0.5,n+0.5),ylim=ylim,col=COLS50,border=COLS,outline=F,names=NAMES,main=MAIN,lwd=2,ylab=YLAB,las=2)
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


# 2021- genes enriched in 2-cell embryos and the pluripotent marker genes used for gene enrichment analysis.xlsx
Totipotency_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/Totipotency_genes.txt")[,1])
Pluripotent_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/Pluripotent_genes.txt")[,1])

plot_genes <- intersect(c(Totipotency_genes, Pluripotent_genes),all_genes)

Apoptosis_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/ApoptosisGenes.txt")[,1])
Translation_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/TranslationGenes.txt")[,1])
###################################################################################################
################################                PLOT               ################################
###################################################################################################

pdf("Fig.MarkerGeneExpressionPattern.pdf",width=6,height=5)
par(mar=c(6,4,4,2),mfrow=c(1,3))
plot_list <- list(logfpkm[Totipotency_genes,tCLC_path[1]],logfpkm[Totipotency_genes,tCLC_path[2]],logfpkm[Totipotency_genes,tCLC_path[3]],logfpkm[Totipotency_genes,tCLC_path[4]])
Boxplot_List(plot_list,"Totipotency",tCLC_path,"log2(fpkm+1)",blue_col,paste(blue_col,"50",sep=""))
plot_list <- list(logfpkm[Totipotency_genes,mESC_path[1]],logfpkm[Totipotency_genes,mESC_path[2]],logfpkm[Totipotency_genes,mESC_path[3]],logfpkm[Totipotency_genes,mESC_path[4]])
Boxplot_List(plot_list,"Totipotency",mESC_path,"log2(fpkm+1)",blue_col,paste(blue_col,"50",sep=""))
Totipotency_fc <- log2(fpkm[Totipotency_genes,tCLC_path[2:4]]/fpkm[Totipotency_genes,mESC_path[2:4]])
boxplot(Totipotency_fc,col=cccol50[c(4,5,6)],border=cccol[c(4,5,6)],main="Totipotency",ylab="log2(2CLC/mESC)",xaxt="n",ylim=c(-10,10),outline=F,las=2);box(lwd=2)
axis(side=1,1:3,c("H18","H30","H43"),las=2)
abline(h=0,lty=2)

par(mar=c(6,4,4,2),mfrow=c(1,3))
plot_list <- list(logfpkm[Pluripotent_genes,tCLC_path[1]],logfpkm[Pluripotent_genes,tCLC_path[2]],logfpkm[Pluripotent_genes,tCLC_path[3]],logfpkm[Pluripotent_genes,tCLC_path[4]])
Boxplot_List(plot_list,"Pluripotent",tCLC_path,"log2(fpkm+1)",blue_col,paste(blue_col,"50",sep=""))
plot_list <- list(logfpkm[Pluripotent_genes,mESC_path[1]],logfpkm[Pluripotent_genes,mESC_path[2]],logfpkm[Pluripotent_genes,mESC_path[3]],logfpkm[Pluripotent_genes,mESC_path[4]])
Boxplot_List(plot_list,"Pluripotent",mESC_path,"log2(fpkm+1)",blue_col,paste(blue_col,"50",sep=""))
Pluripotent_fc <- log2(fpkm[Pluripotent_genes,tCLC_path[2:4]]/fpkm[Pluripotent_genes,mESC_path[2:4]])
boxplot(Pluripotent_fc,col=cccol50[c(4,5,6)],border=cccol[c(4,5,6)],main="Pluripotent",ylab="log2(2CLC/mESC)",xaxt="n",ylim=c(-10,10),outline=F,las=2);box(lwd=2)
axis(side=1,1:3,c("H18","H30","H43"),las=2)
abline(h=0,lty=2)

par(mar=c(6,4,4,2),mfrow=c(1,3))
plot_list <- list(logfpkm[Apoptosis_genes,tCLC_path[1]],logfpkm[Apoptosis_genes,tCLC_path[2]],logfpkm[Apoptosis_genes,tCLC_path[3]],logfpkm[Apoptosis_genes,tCLC_path[4]])
Boxplot_List(plot_list,"Apoptosis",tCLC_path,"log2(fpkm+1)",blue_col,paste(blue_col,"50",sep=""))
plot_list <- list(logfpkm[Apoptosis_genes,mESC_path[1]],logfpkm[Apoptosis_genes,mESC_path[2]],logfpkm[Apoptosis_genes,mESC_path[3]],logfpkm[Apoptosis_genes,mESC_path[4]])
Boxplot_List(plot_list,"Apoptosis",mESC_path,"log2(fpkm+1)",blue_col,paste(blue_col,"50",sep=""))
Apoptosis_fc <- log2(fpkm[Apoptosis_genes,tCLC_path[2:4]]/fpkm[Apoptosis_genes,mESC_path[2:4]])
boxplot(Apoptosis_fc,col=cccol50[c(4,5,6)],border=cccol[c(4,5,6)],main="Apoptosis",ylab="log2(2CLC/mESC)",xaxt="n",ylim=c(-10,10),outline=F,las=2);box(lwd=2)
axis(side=1,1:3,c("H18","H30","H43"),las=2)
abline(h=0,lty=2)

par(mar=c(6,4,4,2),mfrow=c(1,3))
plot_list <- list(logfpkm[Translation_genes,tCLC_path[1]],logfpkm[Translation_genes,tCLC_path[2]],logfpkm[Translation_genes,tCLC_path[3]],logfpkm[Translation_genes,tCLC_path[4]])
Boxplot_List(plot_list,"Translation",tCLC_path,"log2(fpkm+1)",blue_col,paste(blue_col,"50",sep=""))
plot_list <- list(logfpkm[Translation_genes,mESC_path[1]],logfpkm[Translation_genes,mESC_path[2]],logfpkm[Translation_genes,mESC_path[3]],logfpkm[Translation_genes,mESC_path[4]])
Boxplot_List(plot_list,"Translation",mESC_path,"log2(fpkm+1)",blue_col,paste(blue_col,"50",sep=""))
Translation_fc <- log2(fpkm[Translation_genes,tCLC_path[2:4]]/fpkm[Translation_genes,mESC_path[2:4]])
boxplot(Translation_fc,col=cccol50[c(4,5,6)],border=cccol[c(4,5,6)],main="Translation",ylab="log2(2CLC/mESC)",xaxt="n",ylim=c(-10,10),outline=F,las=2);box(lwd=2)
axis(side=1,1:3,c("H18","H30","H43"),las=2)
abline(h=0,lty=2)
dev.off()
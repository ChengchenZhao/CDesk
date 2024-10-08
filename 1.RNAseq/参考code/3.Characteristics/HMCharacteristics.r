###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

cccol1 <- colorRampPalette(c(cccol[1],cccol[5]), bias=1)(7); cccol150 <- paste(cccol1, "50", sep="")
cccol2 <- colorRampPalette(c(cccol[2],cccol[6]), bias=1)(7); cccol250 <- paste(cccol2, "50", sep="")
cccol3 <- colorRampPalette(c(cccol[3],cccol[7]), bias=1)(7); cccol350 <- paste(cccol3, "50", sep="")
cccol4 <- colorRampPalette(c(cccol[4],cccol[8]), bias=1)(7); cccol450 <- paste(cccol4, "50", sep="")

Boxplot_List <- function(PLOTLIST,MAIN,COL,BORDER,NAMES,YLAB){
	n <- length(PLOTLIST)
	xpos <- 0:(n-1)+1.5
	p_value <- c()
	for (i in seq(n-1)){
		tryCatch({p_value <- c(p_value,t.test(PLOTLIST[[i]],PLOTLIST[[i+1]])$p.value)},
			error = function(err){p_value <- c(p_value,1)})
	}
	mark <- symnum(p_value, cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","-"))
	bp <- boxplot(PLOTLIST,at=0:(n-1)+1, xlim=c(0.5,n+0.5),outline=F,plot=F)
	ylim <- range(bp$stats,na.rm=T)
	dist <- (ylim[2]-ylim[1])/20
	ylim[2] <- ylim[2]+1*dist
	boxplot(PLOTLIST,at=0:(n-1)+1,boxwex=0.6,xlim=c(0.5,n+0.5),ylim=ylim,col=COL,border=BORDER,outline=F,names=NAMES,main=MAIN,lwd=2,ylab=YLAB,las=2)
	box(lwd=2)
	ypos_1 <- bp$stats[5,][1:n-1]
	ypos_2 <- bp$stats[5,][2:n]
	# legend("topleft",c("***:p<0.001","**:p<0.01","*:p<0.05"),text.col=cccol[1],bty="n")
	for(i in 1:length(mark)){
		if(!is.na(mark[i])){
			segments(xpos[i]-.45, ypos_1[i]+dist/2, xpos[i]-.45, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]+.45, ypos_2[i]+dist/2, xpos[i]+.45, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]-.45, max(ypos_1[i], ypos_2[i])+dist, xpos[i]-0.25, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]+.45, max(ypos_1[i], ypos_2[i])+dist, xpos[i]+0.25, max(ypos_1[i], ypos_2[i])+dist)
			text(x=xpos[i], y=max(ypos_1[i], ypos_2[i])+dist, label=mark[i], col="black")
		}
	}
	# average line
	tmp_mean <- c()
	for (i in seq(n)){
		tmp_mean <- c(tmp_mean, mean(PLOTLIST[[i]], na.rm=T))
	}
	points(seq(n),tmp_mean,type="b",bg=COL,pch=21,lwd=2)
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
## Chromatin accessibility
MII_acc <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE116857_hgmmEmbryo_ATACseq/GSE116857.mm.MII.ATAC.1.PromoterSignal", row.names=1); MII_acc_norm <- as.data.frame(scale(MII_acc))
Sperm_acc <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE136403_mmEmbryo_ATACseq/GSE136403.mm.sperm.ATAC.PromoterSignal", row.names=1); Sperm_acc_norm <- as.data.frame(scale(Sperm_acc))
Zygote_acc <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE136403_mmEmbryo_ATACseq/GSE136403.mm.zygote.ATAC.PromoterSignal", row.names=1); Zygote_acc_norm <- as.data.frame(scale(Zygote_acc))
C2E_acc <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE66390_mmEmbryo_ATACseq/GSE66390_mm2cellEarly.ATAC.PromoterSignal", row.names=1); C2E_acc_norm <- as.data.frame(scale(C2E_acc))
C2L_acc <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE66390_mmEmbryo_ATACseq/GSE66390_mm2cellLate.ATAC.PromoterSignal", row.names=1); C2L_acc_norm <- as.data.frame(scale(C2L_acc))
C4_acc <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE66390_mmEmbryo_ATACseq/GSE66390_mm4cell.ATAC.PromoterSignal", row.names=1); C4_acc_norm <- as.data.frame(scale(C4_acc))
C8_acc <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE66390_mmEmbryo_ATACseq/GSE66390_mm8cell.ATAC.PromoterSignal", row.names=1); C8_acc_norm <- as.data.frame(scale(C8_acc))

## H3K4me3
MII_k4 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE73952_mmEmbryo_K4K27ChIPseq/GSE73952.mm.MII.K4me3.PromoterSignal", row.names=1); MII_k4_norm <- as.data.frame(scale(MII_k4))
Zygote_k4 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE73952_mmEmbryo_K4K27ChIPseq/GSE73952.mm.zygote.K4me3.PromoterSignal", row.names=1); Zygote_k4_norm <- as.data.frame(scale(Zygote_k4))
C2_k4 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE73952_mmEmbryo_K4K27ChIPseq/GSE73952.mm.2cell.K4me3.PromoterSignal", row.names=1); C2_k4_norm <- as.data.frame(scale(C2_k4))
C4_k4 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE73952_mmEmbryo_K4K27ChIPseq/GSE73952.mm.4cell.K4me3.PromoterSignal", row.names=1); C4_k4_norm <- as.data.frame(scale(C4_k4))
C8_k4 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE73952_mmEmbryo_K4K27ChIPseq/GSE73952.mm.8cell.K4me3.PromoterSignal", row.names=1); C8_k4_norm <- as.data.frame(scale(C8_k4))

## H3K27me3
MII_k27 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE73952_mmEmbryo_K4K27ChIPseq/GSE73952.mm.MII.K27me3.PromoterSignal", row.names=1); MII_k27_norm <- as.data.frame(scale(MII_k27))
Zygote_k27 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE73952_mmEmbryo_K4K27ChIPseq/GSE73952.mm.zygote.K27me3.PromoterSignal", row.names=1); Zygote_k27_norm <- as.data.frame(scale(Zygote_k27))
C2_k27 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE73952_mmEmbryo_K4K27ChIPseq/GSE73952.mm.2cell.K27me3.PromoterSignal", row.names=1); C2_k27_norm <- as.data.frame(scale(C2_k27))
C4_k27 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE73952_mmEmbryo_K4K27ChIPseq/GSE73952.mm.4cell.K27me3.PromoterSignal", row.names=1); C4_k27_norm <- as.data.frame(scale(C4_k27))
C8_k27 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE73952_mmEmbryo_K4K27ChIPseq/GSE73952.mm.8cell.K27me3.PromoterSignal", row.names=1); C8_k27_norm <- as.data.frame(scale(C8_k27))

## H3K9me3
MII_k9 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE97778_mmEmbryo_K9ChIPseq/GSE73952.mm.MII.K9me3.PromoterSignal", row.names=1); MII_k9_norm <- as.data.frame(scale(MII_k9))
Zygote_k9 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE97778_mmEmbryo_K9ChIPseq/GSE73952.mm.zygote.K9me3.PromoterSignal", row.names=1); Zygote_k9_norm <- as.data.frame(scale(Zygote_k9))
C2_k9 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE97778_mmEmbryo_K9ChIPseq/GSE73952.mm.2cell.K9me3.PromoterSignal", row.names=1); C2_k9_norm <- as.data.frame(scale(C2_k9))
C4_k9 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE97778_mmEmbryo_K9ChIPseq/GSE73952.mm.4cell.K9me3.PromoterSignal", row.names=1); C4_k9_norm <- as.data.frame(scale(C4_k9))
C8_k9 <- read.table("/data/zhaochengchen/Work/1.CSwitch/data/GSE97778_mmEmbryo_K9ChIPseq/GSE73952.mm.8cell.K9me3.PromoterSignal", row.names=1); C8_k9_norm <- as.data.frame(scale(C8_k9))

## Genes 
pu_C_vs_dox_C_Up_Totipotency <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/pu_C_vs_dox_C_Up_Totipotency.txt")[,1])
pu_C_vs_dox_C_Dn_Totipotency <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/pu_C_vs_dox_C_Dn_Totipotency.txt")[,1])
pu_C_vs_dox_C_Up_Pluripotent <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/pu_C_vs_dox_C_Up_Pluripotent.txt")[,1])
pu_C_vs_dox_C_Dn_Pluripotent <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/pu_C_vs_dox_C_Dn_Pluripotent.txt")[,1])

Totipotency_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/Totipotency_genes.txt")[,1])
set.seed(1); Totipotency_others <- setdiff(Totipotency_genes,c(pu_C_vs_dox_C_Up_Totipotency,pu_C_vs_dox_C_Dn_Totipotency))
Pluripotent_genes <- as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.ExpressionPattern/Pluripotent_genes.txt")[,1])
set.seed(1); Pluripotent_others <- setdiff(Pluripotent_genes,c(pu_C_vs_dox_C_Up_Pluripotent,pu_C_vs_dox_C_Dn_Pluripotent))

###################################################################################################
################################                PLOT               ################################
###################################################################################################

pdf("HMCharacteristics_Totipotency.pdf",width=7,height=4)
par(mfrow=c(1,3))
plot_list <- list(MII_acc[pu_C_vs_dox_C_Up_Totipotency,],Sperm_acc[pu_C_vs_dox_C_Up_Totipotency,],Zygote_acc[pu_C_vs_dox_C_Up_Totipotency,],C2E_acc[pu_C_vs_dox_C_Up_Totipotency,],C2L_acc[pu_C_vs_dox_C_Up_Totipotency,],C4_acc[pu_C_vs_dox_C_Up_Totipotency,],C8_acc[pu_C_vs_dox_C_Up_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Totipotency",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")
plot_list <- list(MII_acc[pu_C_vs_dox_C_Dn_Totipotency,],Sperm_acc[pu_C_vs_dox_C_Dn_Totipotency,],Zygote_acc[pu_C_vs_dox_C_Dn_Totipotency,],C2E_acc[pu_C_vs_dox_C_Dn_Totipotency,],C2L_acc[pu_C_vs_dox_C_Dn_Totipotency,],C4_acc[pu_C_vs_dox_C_Dn_Totipotency,],C8_acc[pu_C_vs_dox_C_Dn_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Totipotency",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")
plot_list <- list(MII_acc[Totipotency_others,],Sperm_acc[Totipotency_others,],Zygote_acc[Totipotency_others,],C2E_acc[Totipotency_others,],C2L_acc[Totipotency_others,],C4_acc[Totipotency_others,],C8_acc[Totipotency_others,])
Boxplot_List(plot_list,"Totipotency_others",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")

par(mfrow=c(1,3))
plot_list <- list(MII_k4[pu_C_vs_dox_C_Up_Totipotency,],Zygote_k4[pu_C_vs_dox_C_Up_Totipotency,],C2_k4[pu_C_vs_dox_C_Up_Totipotency,],C4_k4[pu_C_vs_dox_C_Up_Totipotency,],C8_k4[pu_C_vs_dox_C_Up_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Totipotency",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")
plot_list <- list(MII_k4[pu_C_vs_dox_C_Dn_Totipotency,],Zygote_k4[pu_C_vs_dox_C_Dn_Totipotency,],C2_k4[pu_C_vs_dox_C_Dn_Totipotency,],C4_k4[pu_C_vs_dox_C_Dn_Totipotency,],C8_k4[pu_C_vs_dox_C_Dn_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Totipotency",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")
plot_list <- list(MII_k4[Totipotency_others,],Zygote_k4[Totipotency_others,],C2_k4[Totipotency_others,],C4_k4[Totipotency_others,],C8_k4[Totipotency_others,])
Boxplot_List(plot_list,"Totipotency_others",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")

par(mfrow=c(1,3))
plot_list <- list(MII_k27[pu_C_vs_dox_C_Up_Totipotency,],Zygote_k27[pu_C_vs_dox_C_Up_Totipotency,],C2_k27[pu_C_vs_dox_C_Up_Totipotency,],C4_k27[pu_C_vs_dox_C_Up_Totipotency,],C8_k27[pu_C_vs_dox_C_Up_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Totipotency",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")
plot_list <- list(MII_k27[pu_C_vs_dox_C_Dn_Totipotency,],Zygote_k27[pu_C_vs_dox_C_Dn_Totipotency,],C2_k27[pu_C_vs_dox_C_Dn_Totipotency,],C4_k27[pu_C_vs_dox_C_Dn_Totipotency,],C8_k27[pu_C_vs_dox_C_Dn_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Totipotency",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")
plot_list <- list(MII_k27[Totipotency_others,],Zygote_k27[Totipotency_others,],C2_k27[Totipotency_others,],C4_k27[Totipotency_others,],C8_k27[Totipotency_others,])
Boxplot_List(plot_list,"Totipotency_others",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")

par(mfrow=c(1,3))
plot_list <- list(MII_k9[pu_C_vs_dox_C_Up_Totipotency,],Zygote_k9[pu_C_vs_dox_C_Up_Totipotency,],C2_k9[pu_C_vs_dox_C_Up_Totipotency,],C4_k9[pu_C_vs_dox_C_Up_Totipotency,],C8_k9[pu_C_vs_dox_C_Up_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Totipotency",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
plot_list <- list(MII_k9[pu_C_vs_dox_C_Dn_Totipotency,],Zygote_k9[pu_C_vs_dox_C_Dn_Totipotency,],C2_k9[pu_C_vs_dox_C_Dn_Totipotency,],C4_k9[pu_C_vs_dox_C_Dn_Totipotency,],C8_k9[pu_C_vs_dox_C_Dn_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Totipotency",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
plot_list <- list(MII_k9[Totipotency_others,],Zygote_k9[Totipotency_others,],C2_k9[Totipotency_others,],C4_k9[Totipotency_others,],C8_k9[Totipotency_others,])
Boxplot_List(plot_list,"Totipotency_others",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
dev.off()


###################################################################################################
pdf("HMCharacteristics.Norm_Totipotency.pdf",width=7,height=4)
par(mfrow=c(1,3))
plot_list <- list(MII_acc_norm[pu_C_vs_dox_C_Up_Totipotency,],Sperm_acc_norm[pu_C_vs_dox_C_Up_Totipotency,],Zygote_acc_norm[pu_C_vs_dox_C_Up_Totipotency,],C2E_acc_norm[pu_C_vs_dox_C_Up_Totipotency,],C2L_acc_norm[pu_C_vs_dox_C_Up_Totipotency,],C4_acc_norm[pu_C_vs_dox_C_Up_Totipotency,],C8_acc_norm[pu_C_vs_dox_C_Up_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Totipotency",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")
plot_list <- list(MII_acc_norm[pu_C_vs_dox_C_Dn_Totipotency,],Sperm_acc_norm[pu_C_vs_dox_C_Dn_Totipotency,],Zygote_acc_norm[pu_C_vs_dox_C_Dn_Totipotency,],C2E_acc_norm[pu_C_vs_dox_C_Dn_Totipotency,],C2L_acc_norm[pu_C_vs_dox_C_Dn_Totipotency,],C4_acc_norm[pu_C_vs_dox_C_Dn_Totipotency,],C8_acc_norm[pu_C_vs_dox_C_Dn_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Totipotency",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")
plot_list <- list(MII_acc_norm[Totipotency_others,],Sperm_acc_norm[Totipotency_others,],Zygote_acc_norm[Totipotency_others,],C2E_acc_norm[Totipotency_others,],C2L_acc_norm[Totipotency_others,],C4_acc_norm[Totipotency_others,],C8_acc_norm[Totipotency_others,])
Boxplot_List(plot_list,"Totipotency_others",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")

par(mfrow=c(1,3))
plot_list <- list(MII_k4_norm[pu_C_vs_dox_C_Up_Totipotency,],Zygote_k4_norm[pu_C_vs_dox_C_Up_Totipotency,],C2_k4_norm[pu_C_vs_dox_C_Up_Totipotency,],C4_k4_norm[pu_C_vs_dox_C_Up_Totipotency,],C8_k4_norm[pu_C_vs_dox_C_Up_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Totipotency",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")
plot_list <- list(MII_k4_norm[pu_C_vs_dox_C_Dn_Totipotency,],Zygote_k4_norm[pu_C_vs_dox_C_Dn_Totipotency,],C2_k4_norm[pu_C_vs_dox_C_Dn_Totipotency,],C4_k4_norm[pu_C_vs_dox_C_Dn_Totipotency,],C8_k4_norm[pu_C_vs_dox_C_Dn_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Totipotency",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")
plot_list <- list(MII_k4_norm[Totipotency_others,],Zygote_k4_norm[Totipotency_others,],C2_k4_norm[Totipotency_others,],C4_k4_norm[Totipotency_others,],C8_k4_norm[Totipotency_others,])
Boxplot_List(plot_list,"Totipotency_others",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")

par(mfrow=c(1,3))
plot_list <- list(MII_k27_norm[pu_C_vs_dox_C_Up_Totipotency,],Zygote_k27_norm[pu_C_vs_dox_C_Up_Totipotency,],C2_k27_norm[pu_C_vs_dox_C_Up_Totipotency,],C4_k27_norm[pu_C_vs_dox_C_Up_Totipotency,],C8_k27_norm[pu_C_vs_dox_C_Up_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Totipotency",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")
plot_list <- list(MII_k27_norm[pu_C_vs_dox_C_Dn_Totipotency,],Zygote_k27_norm[pu_C_vs_dox_C_Dn_Totipotency,],C2_k27_norm[pu_C_vs_dox_C_Dn_Totipotency,],C4_k27_norm[pu_C_vs_dox_C_Dn_Totipotency,],C8_k27_norm[pu_C_vs_dox_C_Dn_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Totipotency",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")
plot_list <- list(MII_k27_norm[Totipotency_others,],Zygote_k27_norm[Totipotency_others,],C2_k27_norm[Totipotency_others,],C4_k27_norm[Totipotency_others,],C8_k27_norm[Totipotency_others,])
Boxplot_List(plot_list,"Totipotency_others",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")

par(mfrow=c(1,3))
plot_list <- list(MII_k9_norm[pu_C_vs_dox_C_Up_Totipotency,],Zygote_k9_norm[pu_C_vs_dox_C_Up_Totipotency,],C2_k9_norm[pu_C_vs_dox_C_Up_Totipotency,],C4_k9_norm[pu_C_vs_dox_C_Up_Totipotency,],C8_k9_norm[pu_C_vs_dox_C_Up_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Totipotency",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
plot_list <- list(MII_k9_norm[pu_C_vs_dox_C_Dn_Totipotency,],Zygote_k9_norm[pu_C_vs_dox_C_Dn_Totipotency,],C2_k9_norm[pu_C_vs_dox_C_Dn_Totipotency,],C4_k9_norm[pu_C_vs_dox_C_Dn_Totipotency,],C8_k9_norm[pu_C_vs_dox_C_Dn_Totipotency,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Totipotency",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
plot_list <- list(MII_k9_norm[Totipotency_others,],Zygote_k9_norm[Totipotency_others,],C2_k9_norm[Totipotency_others,],C4_k9_norm[Totipotency_others,],C8_k9_norm[Totipotency_others,])
Boxplot_List(plot_list,"Totipotency_others",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
dev.off()



###################################################################################################
pdf("HMCharacteristics_Pluripotent.pdf",width=7,height=4)
par(mfrow=c(1,3))
plot_list <- list(MII_acc[pu_C_vs_dox_C_Up_Pluripotent,],Sperm_acc[pu_C_vs_dox_C_Up_Pluripotent,],Zygote_acc[pu_C_vs_dox_C_Up_Pluripotent,],C2E_acc[pu_C_vs_dox_C_Up_Pluripotent,],C2L_acc[pu_C_vs_dox_C_Up_Pluripotent,],C4_acc[pu_C_vs_dox_C_Up_Pluripotent,],C8_acc[pu_C_vs_dox_C_Up_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Pluripotent",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")
plot_list <- list(MII_acc[pu_C_vs_dox_C_Dn_Pluripotent,],Sperm_acc[pu_C_vs_dox_C_Dn_Pluripotent,],Zygote_acc[pu_C_vs_dox_C_Dn_Pluripotent,],C2E_acc[pu_C_vs_dox_C_Dn_Pluripotent,],C2L_acc[pu_C_vs_dox_C_Dn_Pluripotent,],C4_acc[pu_C_vs_dox_C_Dn_Pluripotent,],C8_acc[pu_C_vs_dox_C_Dn_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Pluripotent",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")
plot_list <- list(MII_acc[Pluripotent_others,],Sperm_acc[Pluripotent_others,],Zygote_acc[Pluripotent_others,],C2E_acc[Pluripotent_others,],C2L_acc[Pluripotent_others,],C4_acc[Pluripotent_others,],C8_acc[Pluripotent_others,])
Boxplot_List(plot_list,"Pluripotent_others",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")

par(mfrow=c(1,3))
plot_list <- list(MII_k4[pu_C_vs_dox_C_Up_Pluripotent,],Zygote_k4[pu_C_vs_dox_C_Up_Pluripotent,],C2_k4[pu_C_vs_dox_C_Up_Pluripotent,],C4_k4[pu_C_vs_dox_C_Up_Pluripotent,],C8_k4[pu_C_vs_dox_C_Up_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Pluripotent",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")
plot_list <- list(MII_k4[pu_C_vs_dox_C_Dn_Pluripotent,],Zygote_k4[pu_C_vs_dox_C_Dn_Pluripotent,],C2_k4[pu_C_vs_dox_C_Dn_Pluripotent,],C4_k4[pu_C_vs_dox_C_Dn_Pluripotent,],C8_k4[pu_C_vs_dox_C_Dn_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Pluripotent",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")
plot_list <- list(MII_k4[Pluripotent_others,],Zygote_k4[Pluripotent_others,],C2_k4[Pluripotent_others,],C4_k4[Pluripotent_others,],C8_k4[Pluripotent_others,])
Boxplot_List(plot_list,"Pluripotent_others",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")

par(mfrow=c(1,3))
plot_list <- list(MII_k27[pu_C_vs_dox_C_Up_Pluripotent,],Zygote_k27[pu_C_vs_dox_C_Up_Pluripotent,],C2_k27[pu_C_vs_dox_C_Up_Pluripotent,],C4_k27[pu_C_vs_dox_C_Up_Pluripotent,],C8_k27[pu_C_vs_dox_C_Up_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Pluripotent",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")
plot_list <- list(MII_k27[pu_C_vs_dox_C_Dn_Pluripotent,],Zygote_k27[pu_C_vs_dox_C_Dn_Pluripotent,],C2_k27[pu_C_vs_dox_C_Dn_Pluripotent,],C4_k27[pu_C_vs_dox_C_Dn_Pluripotent,],C8_k27[pu_C_vs_dox_C_Dn_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Pluripotent",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")
plot_list <- list(MII_k27[Pluripotent_others,],Zygote_k27[Pluripotent_others,],C2_k27[Pluripotent_others,],C4_k27[Pluripotent_others,],C8_k27[Pluripotent_others,])
Boxplot_List(plot_list,"Pluripotent_others",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")

par(mfrow=c(1,3))
plot_list <- list(MII_k9[pu_C_vs_dox_C_Up_Pluripotent,],Zygote_k9[pu_C_vs_dox_C_Up_Pluripotent,],C2_k9[pu_C_vs_dox_C_Up_Pluripotent,],C4_k9[pu_C_vs_dox_C_Up_Pluripotent,],C8_k9[pu_C_vs_dox_C_Up_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Pluripotent",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
plot_list <- list(MII_k9[pu_C_vs_dox_C_Dn_Pluripotent,],Zygote_k9[pu_C_vs_dox_C_Dn_Pluripotent,],C2_k9[pu_C_vs_dox_C_Dn_Pluripotent,],C4_k9[pu_C_vs_dox_C_Dn_Pluripotent,],C8_k9[pu_C_vs_dox_C_Dn_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Pluripotent",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
plot_list <- list(MII_k9[Pluripotent_others,],Zygote_k9[Pluripotent_others,],C2_k9[Pluripotent_others,],C4_k9[Pluripotent_others,],C8_k9[Pluripotent_others,])
Boxplot_List(plot_list,"Pluripotent_others",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
dev.off()


###################################################################################################
pdf("HMCharacteristics.Norm_Pluripotent.pdf",width=7,height=4)
par(mfrow=c(1,3))
plot_list <- list(MII_acc_norm[pu_C_vs_dox_C_Up_Pluripotent,],Sperm_acc_norm[pu_C_vs_dox_C_Up_Pluripotent,],Zygote_acc_norm[pu_C_vs_dox_C_Up_Pluripotent,],C2E_acc_norm[pu_C_vs_dox_C_Up_Pluripotent,],C2L_acc_norm[pu_C_vs_dox_C_Up_Pluripotent,],C4_acc_norm[pu_C_vs_dox_C_Up_Pluripotent,],C8_acc_norm[pu_C_vs_dox_C_Up_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Pluripotent",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")
plot_list <- list(MII_acc_norm[pu_C_vs_dox_C_Dn_Pluripotent,],Sperm_acc_norm[pu_C_vs_dox_C_Dn_Pluripotent,],Zygote_acc_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C2E_acc_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C2L_acc_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C4_acc_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C8_acc_norm[pu_C_vs_dox_C_Dn_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Pluripotent",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")
plot_list <- list(MII_acc_norm[Pluripotent_others,],Sperm_acc_norm[Pluripotent_others,],Zygote_acc_norm[Pluripotent_others,],C2E_acc_norm[Pluripotent_others,],C2L_acc_norm[Pluripotent_others,],C4_acc_norm[Pluripotent_others,],C8_acc_norm[Pluripotent_others,])
Boxplot_List(plot_list,"Pluripotent_others",cccol150,cccol1,c("MII","Sperm","Zygote","C2E","C2L","C4","C8"),"Chromatin accessibility")

par(mfrow=c(1,3))
plot_list <- list(MII_k4_norm[pu_C_vs_dox_C_Up_Pluripotent,],Zygote_k4_norm[pu_C_vs_dox_C_Up_Pluripotent,],C2_k4_norm[pu_C_vs_dox_C_Up_Pluripotent,],C4_k4_norm[pu_C_vs_dox_C_Up_Pluripotent,],C8_k4_norm[pu_C_vs_dox_C_Up_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Pluripotent",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")
plot_list <- list(MII_k4_norm[pu_C_vs_dox_C_Dn_Pluripotent,],Zygote_k4_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C2_k4_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C4_k4_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C8_k4_norm[pu_C_vs_dox_C_Dn_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Pluripotent",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")
plot_list <- list(MII_k4_norm[Pluripotent_others,],Zygote_k4_norm[Pluripotent_others,],C2_k4_norm[Pluripotent_others,],C4_k4_norm[Pluripotent_others,],C8_k4_norm[Pluripotent_others,])
Boxplot_List(plot_list,"Pluripotent_others",cccol250,cccol2,c("MII","Zygote","C2","C4","C8"),"H3K4me3 signal")

par(mfrow=c(1,3))
plot_list <- list(MII_k27_norm[pu_C_vs_dox_C_Up_Pluripotent,],Zygote_k27_norm[pu_C_vs_dox_C_Up_Pluripotent,],C2_k27_norm[pu_C_vs_dox_C_Up_Pluripotent,],C4_k27_norm[pu_C_vs_dox_C_Up_Pluripotent,],C8_k27_norm[pu_C_vs_dox_C_Up_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Pluripotent",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")
plot_list <- list(MII_k27_norm[pu_C_vs_dox_C_Dn_Pluripotent,],Zygote_k27_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C2_k27_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C4_k27_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C8_k27_norm[pu_C_vs_dox_C_Dn_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Pluripotent",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")
plot_list <- list(MII_k27_norm[Pluripotent_others,],Zygote_k27_norm[Pluripotent_others,],C2_k27_norm[Pluripotent_others,],C4_k27_norm[Pluripotent_others,],C8_k27_norm[Pluripotent_others,])
Boxplot_List(plot_list,"Pluripotent_others",cccol350,cccol3,c("MII","Zygote","C2","C4","C8"),"H3K27me3 signal")

par(mfrow=c(1,3))
plot_list <- list(MII_k9_norm[pu_C_vs_dox_C_Up_Pluripotent,],Zygote_k9_norm[pu_C_vs_dox_C_Up_Pluripotent,],C2_k9_norm[pu_C_vs_dox_C_Up_Pluripotent,],C4_k9_norm[pu_C_vs_dox_C_Up_Pluripotent,],C8_k9_norm[pu_C_vs_dox_C_Up_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Up_Pluripotent",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
plot_list <- list(MII_k9_norm[pu_C_vs_dox_C_Dn_Pluripotent,],Zygote_k9_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C2_k9_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C4_k9_norm[pu_C_vs_dox_C_Dn_Pluripotent,],C8_k9_norm[pu_C_vs_dox_C_Dn_Pluripotent,])
Boxplot_List(plot_list,"pu_C_vs_dox_C_Dn_Pluripotent",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
plot_list <- list(MII_k9_norm[Pluripotent_others,],Zygote_k9_norm[Pluripotent_others,],C2_k9_norm[Pluripotent_others,],C4_k9_norm[Pluripotent_others,],C8_k9_norm[Pluripotent_others,])
Boxplot_List(plot_list,"Pluripotent_others",cccol450,cccol4,c("MII","Zygote","C2","C4","C8"),"H3K9me3 signal")
dev.off()

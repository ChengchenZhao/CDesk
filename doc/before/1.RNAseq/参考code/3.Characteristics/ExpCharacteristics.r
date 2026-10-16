###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)
library(beeswarm)

GeneExpressionPlot <-function(gene_list,out_name){
    pdf(paste(out_name,".pdf",sep=""),width=10,height=3.4)
    for (each in gene_list){
        par(mar=c(6,4,4,2))
        layout(matrix(c(rep(1,7),rep(2,3)),nrow=1,ncol=10))
        plot_list <- list(apply(logfpkm[each,Zygote],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e1],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e2],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e3],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e4],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e5],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e6],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e7],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e8],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e9],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e10],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e11],2,mean,na.rm=T),
                  apply(logfpkm[each,C2_e12],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e1],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e2],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e3],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e4],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e5],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e6],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e7],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e8],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e9],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e10],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e11],2,mean,na.rm=T),
                  apply(logfpkm[each,C4_e12],2,mean,na.rm=T),
                  apply(logfpkm[each,LC4_e1],2,mean,na.rm=T),
                  apply(logfpkm[each,LC4_e2],2,mean,na.rm=T),
                  apply(logfpkm[each,LC4_e3],2,mean,na.rm=T),
                  apply(logfpkm[each,LC4_e4],2,mean,na.rm=T),
                  apply(logfpkm[each,LC4_e5],2,mean,na.rm=T),
                  apply(logfpkm[each,LC4_e6],2,mean,na.rm=T),
                  apply(logfpkm[each,C8_e1],2,mean,na.rm=T),
                  apply(logfpkm[each,C8_e2],2,mean,na.rm=T),
                  apply(logfpkm[each,C8_e3],2,mean,na.rm=T),
                  apply(logfpkm[each,C8_e4],2,mean,na.rm=T),
                  apply(logfpkm[each,C16_e1],2,mean,na.rm=T),
                  apply(logfpkm[each,C16_e2],2,mean,na.rm=T),
                  apply(logfpkm[each,C16_e3],2,mean,na.rm=T),
                  apply(logfpkm[each,C32_e1],2,mean,na.rm=T),
                  apply(logfpkm[each,C32_e2],2,mean,na.rm=T),
                  apply(logfpkm[each,C32_e4],2,mean,na.rm=T),
                  apply(logfpkm[each,C32_e6],2,mean,na.rm=T),
                  apply(logfpkm[each,C32_eMix],2,mean,na.rm=T))
        beeswarm(plot_list,xaxt="n",col = cccol[c(9,rep(1,12),rep(2,12),rep(3,6),rep(4,4),rep(5,3),rep(6,5))],pch = 16,method = "swarm",corral = "wrap",main=each,xlab="",ylab="log(fpkm+1)",ylim=c(0,12),las=2)
        box(lwd=2)
        axis(side=1,at=1:length(plot_list),labels=embryo_labels,las=2)
        plot_mean <- mean_matrix[each,]
        points(seq(length(plot_mean)),plot_mean,type="b",col="black",pch="*",lwd=2)
        plot_vector <- CV2[each,]
        plot(seq(length(plot_vector)),plot_vector,xaxt="n",col=cccol[c(9,rep(1,12),rep(2,12),rep(3,6),rep(4,4),rep(5,3),rep(6,5))],type="h",lwd=3,xlab="",ylab="CV2",main="",ylim=c(0,12))
        axis(side=1,at=1:length(plot_vector),labels=embryo_labels,las=2)
        box(lwd=2)
    }
    dev.off()
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
## expression
fpkm <- read.table("/data/zhaochengchen/Work/1.CSwitch/results/0.AverageSignal/GSE136714_RNAseqFPKM_Combined.txt",header=T,row.names=1)
all_cells <- colnames(fpkm)
logfpkm <- log2(fpkm+1)

Zygote <- all_cells[grep("^Zygote",all_cells)]
C2_e1 <- all_cells[grep("^C2_e1([^0-9]|$)",all_cells)]
C2_e2 <- all_cells[grep("^C2_e2([^0-9]|$)",all_cells)]
C2_e3 <- all_cells[grep("^C2_e3([^0-9]|$)",all_cells)]
C2_e4 <- all_cells[grep("^C2_e4([^0-9]|$)",all_cells)]
C2_e5 <- all_cells[grep("^C2_e5([^0-9]|$)",all_cells)]
C2_e6 <- all_cells[grep("^C2_e6([^0-9]|$)",all_cells)]
C2_e7 <- all_cells[grep("^C2_e7([^0-9]|$)",all_cells)]
C2_e8 <- all_cells[grep("^C2_e8([^0-9]|$)",all_cells)]
C2_e9 <- all_cells[grep("^C2_e9([^0-9]|$)",all_cells)]
C2_e10 <- all_cells[grep("^C2_e10([^0-9]|$)",all_cells)]
C2_e11 <- all_cells[grep("^C2_e11([^0-9]|$)",all_cells)]
C2_e12 <- all_cells[grep("^C2_e12([^0-9]|$)",all_cells)]
C4_e1 <- all_cells[grep("^C4_e1([^0-9]|$)",all_cells)]
C4_e2 <- all_cells[grep("^C4_e2([^0-9]|$)",all_cells)]
C4_e3 <- all_cells[grep("^C4_e3([^0-9]|$)",all_cells)]
C4_e4 <- all_cells[grep("^C4_e4([^0-9]|$)",all_cells)]
C4_e5 <- all_cells[grep("^C4_e5([^0-9]|$)",all_cells)]
C4_e6 <- all_cells[grep("^C4_e6([^0-9]|$)",all_cells)]
C4_e7 <- all_cells[grep("^C4_e7([^0-9]|$)",all_cells)]
C4_e8 <- all_cells[grep("^C4_e8([^0-9]|$)",all_cells)]
C4_e9 <- all_cells[grep("^C4_e9([^0-9]|$)",all_cells)]
C4_e10 <- all_cells[grep("^C4_e10([^0-9]|$)",all_cells)]
C4_e11 <- all_cells[grep("^C4_e11([^0-9]|$)",all_cells)]
C4_e12 <- all_cells[grep("^C4_e12([^0-9]|$)",all_cells)]
LC4_e1 <- all_cells[grep("^LC4_e1([^0-9]|$)",all_cells)]
LC4_e2 <- all_cells[grep("^LC4_e2([^0-9]|$)",all_cells)]
LC4_e3 <- all_cells[grep("^LC4_e3([^0-9]|$)",all_cells)]
LC4_e4 <- all_cells[grep("^LC4_e4([^0-9]|$)",all_cells)]
LC4_e5 <- all_cells[grep("^LC4_e5([^0-9]|$)",all_cells)]
LC4_e6 <- all_cells[grep("^LC4_e6([^0-9]|$)",all_cells)]
C8_e1 <- all_cells[grep("^C8_e1([^0-9]|$)",all_cells)]
C8_e2 <- all_cells[grep("^C8_e2([^0-9]|$)",all_cells)]
C8_e3 <- all_cells[grep("^C8_e3([^0-9]|$)",all_cells)]
C8_e4 <- all_cells[grep("^C8_e4([^0-9]|$)",all_cells)]
C16_e1 <- all_cells[grep("^C16_e1([^0-9]|$)",all_cells)]
C16_e2 <- all_cells[grep("^C16_e2([^0-9]|$)",all_cells)]
C16_e3 <- all_cells[grep("^C16_e3([^0-9]|$)",all_cells)]
C32_e1 <- all_cells[grep("^C32_e1([^0-9]|$)",all_cells)]
C32_e2 <- all_cells[grep("^C32_e2([^0-9]|$)",all_cells)]
C32_e4 <- all_cells[grep("^C32_e4([^0-9]|$)",all_cells)]
C32_e6 <- all_cells[grep("^C32_e6([^0-9]|$)",all_cells)]
C32_eMix <- all_cells[grep("^C32_eMix([^0-9]|$)",all_cells)]

mean_matrix <- cbind(apply(logfpkm[,Zygote],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e1],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e2],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e3],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e4],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e5],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e6],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e7],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e8],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e9],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e10],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e11],1,mean,na.rm=T),
                     apply(logfpkm[,C2_e12],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e1],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e2],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e3],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e4],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e5],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e6],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e7],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e8],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e9],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e10],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e11],1,mean,na.rm=T),
                     apply(logfpkm[,C4_e12],1,mean,na.rm=T),
                     apply(logfpkm[,LC4_e1],1,mean,na.rm=T),
                     apply(logfpkm[,LC4_e2],1,mean,na.rm=T),
                     apply(logfpkm[,LC4_e3],1,mean,na.rm=T),
                     apply(logfpkm[,LC4_e4],1,mean,na.rm=T),
                     apply(logfpkm[,LC4_e5],1,mean,na.rm=T),
                     apply(logfpkm[,LC4_e6],1,mean,na.rm=T),
                     apply(logfpkm[,C8_e1],1,mean,na.rm=T),
                     apply(logfpkm[,C8_e2],1,mean,na.rm=T),
                     apply(logfpkm[,C8_e3],1,mean,na.rm=T),
                     apply(logfpkm[,C8_e4],1,mean,na.rm=T),
                     apply(logfpkm[,C16_e1],1,mean,na.rm=T),
                     apply(logfpkm[,C16_e2],1,mean,na.rm=T),
                     apply(logfpkm[,C16_e3],1,mean,na.rm=T),
                     apply(logfpkm[,C32_e1],1,mean,na.rm=T),
                     apply(logfpkm[,C32_e2],1,mean,na.rm=T),
                     apply(logfpkm[,C32_e4],1,mean,na.rm=T),
                     apply(logfpkm[,C32_e6],1,mean,na.rm=T),
                     apply(logfpkm[,C32_eMix],1,mean,na.rm=T))
var_matrix <-  cbind(apply(logfpkm[,Zygote],1,var,na.rm=T),
                     apply(logfpkm[,C2_e1],1,var,na.rm=T),
                     apply(logfpkm[,C2_e2],1,var,na.rm=T),
                     apply(logfpkm[,C2_e3],1,var,na.rm=T),
                     apply(logfpkm[,C2_e4],1,var,na.rm=T),
                     apply(logfpkm[,C2_e5],1,var,na.rm=T),
                     apply(logfpkm[,C2_e6],1,var,na.rm=T),
                     apply(logfpkm[,C2_e7],1,var,na.rm=T),
                     apply(logfpkm[,C2_e8],1,var,na.rm=T),
                     apply(logfpkm[,C2_e9],1,var,na.rm=T),
                     apply(logfpkm[,C2_e10],1,var,na.rm=T),
                     apply(logfpkm[,C2_e11],1,var,na.rm=T),
                     apply(logfpkm[,C2_e12],1,var,na.rm=T),
                     apply(logfpkm[,C4_e1],1,var,na.rm=T),
                     apply(logfpkm[,C4_e2],1,var,na.rm=T),
                     apply(logfpkm[,C4_e3],1,var,na.rm=T),
                     apply(logfpkm[,C4_e4],1,var,na.rm=T),
                     apply(logfpkm[,C4_e5],1,var,na.rm=T),
                     apply(logfpkm[,C4_e6],1,var,na.rm=T),
                     apply(logfpkm[,C4_e7],1,var,na.rm=T),
                     apply(logfpkm[,C4_e8],1,var,na.rm=T),
                     apply(logfpkm[,C4_e9],1,var,na.rm=T),
                     apply(logfpkm[,C4_e10],1,var,na.rm=T),
                     apply(logfpkm[,C4_e11],1,var,na.rm=T),
                     apply(logfpkm[,C4_e12],1,var,na.rm=T),
                     apply(logfpkm[,LC4_e1],1,var,na.rm=T),
                     apply(logfpkm[,LC4_e2],1,var,na.rm=T),
                     apply(logfpkm[,LC4_e3],1,var,na.rm=T),
                     apply(logfpkm[,LC4_e4],1,var,na.rm=T),
                     apply(logfpkm[,LC4_e5],1,var,na.rm=T),
                     apply(logfpkm[,LC4_e6],1,var,na.rm=T),
                     apply(logfpkm[,C8_e1],1,var,na.rm=T),
                     apply(logfpkm[,C8_e2],1,var,na.rm=T),
                     apply(logfpkm[,C8_e3],1,var,na.rm=T),
                     apply(logfpkm[,C8_e4],1,var,na.rm=T),
                     apply(logfpkm[,C16_e1],1,var,na.rm=T),
                     apply(logfpkm[,C16_e2],1,var,na.rm=T),
                     apply(logfpkm[,C16_e3],1,var,na.rm=T),
                     apply(logfpkm[,C32_e1],1,var,na.rm=T),
                     apply(logfpkm[,C32_e2],1,var,na.rm=T),
                     apply(logfpkm[,C32_e4],1,var,na.rm=T),
                     apply(logfpkm[,C32_e6],1,var,na.rm=T),
                     apply(logfpkm[,C32_eMix],1,var,na.rm=T))
CV2 <- var_matrix/(mean_matrix)^2
embryo_labels <- c("Zygote","C2_e1","C2_e2","C2_e3","C2_e4","C2_e5","C2_e6","C2_e7","C2_e8","C2_e9","C2_e10","C2_e11","C2_e12","C4_e1","C4_e2","C4_e3","C4_e4","C4_e5","C4_e6","C4_e7","C4_e8","C4_e9","C4_e10","C4_e11","C4_e12","LC4_e1","LC4_e2","LC4_e3","LC4_e4","LC4_e5","LC4_e6","C8_e1","C8_e2","C8_e3","C8_e4","C16_e1","C16_e2","C16_e3","C32_e1","C32_e2","C32_e4","C32_e6","C32_eMix")
colnames(CV2) <- embryo_labels
colnames(var_matrix) <- embryo_labels
colnames(mean_matrix) <- embryo_labels

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
GeneExpressionPlot(pu_C_vs_dox_C_Up_Totipotency,"ExpCharacteristics.pu_C_vs_dox_C_Up_Totipotency")
GeneExpressionPlot(pu_C_vs_dox_C_Dn_Totipotency,"ExpCharacteristics.pu_C_vs_dox_C_Dn_Totipotency")
GeneExpressionPlot(pu_C_vs_dox_C_Up_Pluripotent,"ExpCharacteristics.pu_C_vs_dox_C_Up_Pluripotent")
GeneExpressionPlot(pu_C_vs_dox_C_Dn_Pluripotent,"ExpCharacteristics.pu_C_vs_dox_C_Dn_Pluripotent")


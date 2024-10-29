###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
sys_argv <- commandArgs(T)
# sys_argv <- c("/data/zhaochengchen/Work/2.CEmbryo/results/0.mapping/scRNAseq/matrix/EGKX_0h_WE_c1_qcmat.txt","/data/zhaochengchen/Work/2.CEmbryo/results/0.mapping/scRNAseq/matrix/EGKX_0h_WE_c1_expmat.txt","EGKX_0h_WE_c1")
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
data <- read.table(sys_argv[1],header=T,row.names=1)
exp <- read.table(sys_argv[2],header=T,row.names=1)
out_name <- sys_argv[3]
###################################################################################################
################################                PLOT               ################################
###################################################################################################

reads_num <- sort(data[,2],decreasing=T)

reads_num_cut <- reads_num[1:min(10000,length(reads_num))] # excluding the > 10000 cells

# y_norm <- (log10(reads_num_cut)-min(log10(reads_num_cut)))/(max(log10(reads_num_cut))-min(log10(reads_num_cut)))
# x_norm <- seq(length(reads_num_cut))/length(reads_num_cut)
# lenscore <- 10^9
# for(each in seq(reads_num_cut)){
# 	tmp_score <- (1-x_norm[each])^2+(1-y_norm[each])^2
# 	if (lenscore > tmp_score){
# 		lenscore <- tmp_score
# 		lenid <- each
# 	}
# }


tmp_score_list <- c()
y_norm <- (log10(reads_num_cut)-min(log10(reads_num_cut)))/(max(log10(reads_num_cut))-min(log10(reads_num_cut)))
if (length(reads_num) < 10000){
	y_norm <- c(y_norm, rep(0,10000-length(reads_num)))
}
x_norm <- seq(10000)/10000

lenscore <- 10^9
for(each in seq(10000)){
	tmp_score <- (1-x_norm[each])^2+(1-y_norm[each])^2
	if (lenscore > tmp_score){
		lenscore <- tmp_score
		lenid <- each
	}
	tmp_score_list <- c(tmp_score_list,tmp_score)
}

# pdf(paste(out_name,"_SelecetedCellsThresholdScore.pdf",sep=""),width=4.3,height=4.5)
# plot(x_norm,y_norm,type="l")
# plot(x_norm,tmp_score_list,type="l")
# dev.off()

# pdf(paste(out_name,"_CellReadsDistribution.pdf",sep=""),width=4.3,height=4.5)
# plot(seq(length(reads_num)),log10(reads_num),type="l",col=cccol[2],lwd=2,xlab="Cell barcode id",ylab="Reads number(log10)",main=out_name);box(lwd=2)
# points(seq(4),log10(reads_num[1:4]),pch=16,col=cccol[1])
# text(seq(4)+(length(reads_num)/12),log10(reads_num[1:4]),round(log10(reads_num[1:4]),2))
# abline(v=lenid,lty=2,lwd=2)
# text(lenid,max(log10(reads_num))*0.9,lenid,col=cccol[1])
# print(out_name)
# print(lenid)

# plot(seq(length(reads_num_cut)),log10(reads_num_cut),type="l",col=cccol[2],lwd=2,xlab="Cell barcode id",ylab="Reads number(log10)",main=out_name);box(lwd=2)
# points(seq(4),log10(reads_num_cut[1:4]),pch=16,col=cccol[1])
# text(seq(4)+(length(reads_num_cut)/12),log10(reads_num_cut[1:4]),round(log10(reads_num_cut[1:4]),2))
# abline(v=lenid,lty=2,lwd=2)
# text(lenid,max(log10(reads_num_cut))*0.9,lenid,col=cccol[1])
# dev.off()

selected_cells <- row.names(data)[which(data[,2]>=reads_num[lenid])]
write.table(selected_cells,file=paste(out_name,"_SelectedCellBarcodes.txt",sep=""),quote=F,sep="\t",col.names=F,row.names=F)

# out_exp <- exp[,selected_cells]
# colnames(out_exp) <- paste(out_name,"_",selected_cells,sep="")
# write.table(cbind(colnames(out_exp)),file=paste(out_name,"_SelectedCells.txt",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
# write.table(out_exp,file=paste(out_name,"_SelectedCellsExpression.txt",sep=""),quote=F,sep="\t",col.names=T,row.names=T)

# GN_cut <- 1000
# GNselected_cells <- row.names(data)[which(data[,9]>GN_cut)]
# out_exp <- exp[,GNselected_cells]
# colnames(out_exp) <- paste(out_name,"_",GNselected_cells,sep="")
# write.table(cbind(colnames(out_exp)),file=paste(out_name,"_GNSelectedCells.txt",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
# write.table(out_exp,file=paste(out_name,"_GNSelectedCellsExpression.txt",sep=""),quote=F,sep="\t",col.names=T,row.names=T)

# pdf(paste(out_name,"_CoveredGeneNumberDistribution.pdf",sep=""),width=2.6,height=4.5)
# boxplot(as.numeric(data[which(data[,2]>=reads_num[lenid]),9]),col=cccol50[1],border=cccol[1],outline=F)
# boxplot(as.numeric(data[which(data[,9]>GN_cut),9]),col=cccol50[1],border=cccol[1],outline=F)
# dev.off()

############################################################################
########################	   function define       #######################
############################################################################
# library(amap)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")

############################################################################
#########################	   read in data       ##########################
############################################################################

############ Expression matrix
Counts <- read.table("inhouse_RNAseq_Counts_20211029.txt",header=T,row.names=1)
all_samples <- c("CC_10h_r1","CC_10h_r2","CC_10h_r3","CC_12h_r1","CC_12h_r2","CC_14h_r1","CC_14h_r2","CC_14h_r3","CC_16h_r1","CC_16h_r2","CC_16h_r3","CC_2h_r1","CC_2h_r2","CC_2h_r3","CC_4h_r1","CC_4h_r2","CC_4h_r3","CC_6h_r1","CC_6h_r2","CC_8h_r1","CC_8h_r2","CC_8h_r3","CC_ctrl_r1","CC_ctrl_r2","CC_ctrl_r3")
colnames(Counts) <- all_samples
all_genes <- row.names(Counts)
write.table(all_genes,file="AllGenes.txt",quote=F,col.names=F,row.names=F)

CC_ctrl <- c("CC_ctrl_r1","CC_ctrl_r2")
CC_14h <- all_samples[grep("_14h_",all_samples)]

avgCounts <- cbind(apply(Counts[all_genes,CC_ctrl],1,mean),apply(Counts[all_genes,CC_14h],1,mean))
colnames(avgCounts) <- c("CC_ctrl","CC_14h")
rownames(avgCounts) <- all_genes

# genes
GAP <- as.vector(read.table("GAP.OverlappedGenes.txt")[,1])
GAP_e <- as.vector(read.table("GAP.OverlappedGenes.Enhancer.txt")[,1])

############################################################################
####################                plot               #####################
############################################################################

pdf("Fig2d.CorrespondingGenesMAPlot.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
x <- avgCounts[all_genes,"CC_14h"]
y <- avgCounts[all_genes,"CC_ctrl"]
M <- log2(x/y)
A <- (log2(x)+log2(y))/2
plot(A,M,pch=".",col="grey",xlab="A (log2 mean expression)",ylab="M (log2 fold change)",ylim=c(-10,10));box(lwd=2)
abline(h=c(0),lty=2)
points(A[GAP],M[GAP],pch=".",col=cccol[1])
legend("topleft",legend=c(sum(M[GAP]>0,na.rm=T)),col=cccol[1],bty="n")
legend("bottomleft",legend=c(sum(M[GAP]<0,na.rm=T)),col=cccol[1],bty="n")
dev.off()

pdf("Fig2d.CorrespondingGenesMAPlot.Enhancer.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
x <- avgCounts[all_genes,"CC_14h"]
y <- avgCounts[all_genes,"CC_ctrl"]
M <- log2(x/y)
A <- (log2(x)+log2(y))/2
plot(A,M,pch=".",col="grey",xlab="A (log2 mean expression)",ylab="M (log2 fold change)",ylim=c(-10,10));box(lwd=2)
abline(h=c(0),lty=2)
points(A[GAP_e],M[GAP_e],pch=".",col=cccol[1])
legend("topleft",legend=c(sum(M[GAP_e]>0,na.rm=T)),col=cccol[1],bty="n")
legend("bottomleft",legend=c(sum(M[GAP_e]<0,na.rm=T)),col=cccol[1],bty="n")
dev.off()
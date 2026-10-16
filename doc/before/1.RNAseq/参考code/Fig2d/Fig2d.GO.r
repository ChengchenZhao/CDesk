###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
# genes
GAP <- as.vector(read.table("GAP.OverlappedGenes.txt")[,1])
GAP_e <- as.vector(read.table("GAP.OverlappedGenes.Enhancer.txt")[,1])

set.seed(1)
write.table(sample(GAP,3000),file="GAP.OverlappedGenes.Promoter.sample1.txt",quote=F,col.names=F,row.names=F)
set.seed(2)
write.table(sample(GAP,3000),file="GAP.OverlappedGenes.Promoter.sample2.txt",quote=F,col.names=F,row.names=F)
set.seed(3)
write.table(sample(GAP,3000),file="GAP.OverlappedGenes.Promoter.sample3.txt",quote=F,col.names=F,row.names=F)

set.seed(1)
write.table(sample(GAP_e,3000),file="GAP.OverlappedGenes.Enhancer.sample1.txt",quote=F,col.names=F,row.names=F)
set.seed(2)
write.table(sample(GAP_e,3000),file="GAP.OverlappedGenes.Enhancer.sample2.txt",quote=F,col.names=F,row.names=F)
set.seed(3)
write.table(sample(GAP_e,3000),file="GAP.OverlappedGenes.Enhancer.sample3.txt",quote=F,col.names=F,row.names=F)

GO_Promoter <- read.table("Fig2d.GO.Promoter.txt",sep='\t',header=T)
GO_Enhancer <- read.table("Fig2d.GO.Enhancer.txt",sep='\t',header=T)

###################################################################################################
################################                PLOT               ################################
###################################################################################################
pdf("Fig2d.GO.Promoter.pdf",width=10,height=4)
GO_Promoter[,3] <- -log10(GO_Promoter[,3])
GO_Promoter <- GO_Promoter[order(GO_Promoter[,3]),]
GO_Promoter <- GO_Promoter[!duplicated(GO_Promoter),]
par(mar=c(5,30,4,5))
barplot(GO_Promoter[,3],horiz=T,names=GO_Promoter[,2],las=2,border=F,col="grey",space=0.5,axes=F)
axis(side=1);mtext("-log10(p-value)",side=1,line=2)
dev.off()

pdf("Fig2d.GO.Enhancer.pdf",width=10,height=4)
GO_Enhancer[,3] <- -log10(GO_Enhancer[,3])
GO_Enhancer <- GO_Enhancer[order(GO_Enhancer[,3]),]
GO_Enhancer <- GO_Enhancer[!duplicated(GO_Enhancer),]
par(mar=c(5,30,4,5))
barplot(GO_Enhancer[,3],horiz=T,names=GO_Enhancer[,2],las=2,border=F,col="grey",space=0.5,axes=F)
axis(side=1);mtext("-log10(p-value)",side=1,line=2)
dev.off()
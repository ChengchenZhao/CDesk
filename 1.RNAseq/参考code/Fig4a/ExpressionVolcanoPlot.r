###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)
library(qvalue)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
fpkm <- read.table("../../Data/inhouse_scRNAseq_FPKM.txt",header=T,row.names=1)
all_samples <- c("Ctrl_1C14hpa_r1","Ctrl_1C14hpa_r2","Ctrl_1C14hpa_r3","Ctrl_1C6hpa_r1","Ctrl_1C6hpa_r2","Ctrl_1C6hpa_r3","Ctrl_4C_r1","Ctrl_4C_r2","Ctrl_4C_r3","Ctrl_8C_r1","Ctrl_8C_r2","Ctrl_8C_r3","Ctrl_CC_r1","Ctrl_CC_r2","Ctrl_CC_r3","Ctrl_E2C_r1","Ctrl_E2C_r2","Ctrl_E2C_r3","Ctrl_ICM_r1","Ctrl_ICM_r2","Ctrl_ICM_r3","Ctrl_L2C_r1","Ctrl_L2C_r2","Ctrl_L2C_r3","Ctrl_Morula_r1","Ctrl_Morula_r2","Ctrl_Morula_r3","Ctrl_TE_r1","Ctrl_TE_r2","Ctrl_TE_r3","M14h_1C14hpa_r1","M14h_1C14hpa_r2","M14h_4C_r1","M14h_4C_r2","M14h_CC_r1","M14h_CC_r2","M14h_E2C_r1","M14h_E2C_r2","M14h_L2C_r1","M14h_L2C_r2")
colnames(fpkm) <- all_samples

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

Dev_fpkm <- read.table("../../Data/GSE71434_RNAseq_mmEmbryo_FPKM.txt",header=T,row.names=1)
Dev_MII <- colnames(Dev_fpkm)[grep("MII",colnames(Dev_fpkm))]
Dev_1C <- colnames(Dev_fpkm)[grep("zygote",colnames(Dev_fpkm))]
Dev_E2C <- colnames(Dev_fpkm)[grep("EC2",colnames(Dev_fpkm))]
Dev_L2C <- colnames(Dev_fpkm)[grep("LC2",colnames(Dev_fpkm))]
Dev_4C <- colnames(Dev_fpkm)[grep("C4",colnames(Dev_fpkm))]
Dev_8C <- colnames(Dev_fpkm)[grep("C8",colnames(Dev_fpkm))]
Dev_ICM <- colnames(Dev_fpkm)[grep("ICM",colnames(Dev_fpkm))]

all_genes <- intersect(rownames(fpkm),rownames(Dev_fpkm))

Ctrl_avgfpkm <- cbind(apply(fpkm[all_genes,Ctrl_CC],1,mean),apply(fpkm[all_genes,Ctrl_1C14hpa],1,mean),apply(fpkm[all_genes,Ctrl_E2C],1,mean),apply(fpkm[all_genes,Ctrl_L2C],1,mean),apply(fpkm[all_genes,Ctrl_4C],1,mean))
colnames(Ctrl_avgfpkm) <- c("Ctrl_CC","Ctrl_1C","Ctrl_E2C","Ctrl_L2C","Ctrl_4C");rownames(Ctrl_avgfpkm) <- all_genes
M14h_avgfpkm <- cbind(apply(fpkm[all_genes,M14h_CC],1,mean),apply(fpkm[all_genes,M14h_1C14hpa],1,mean),apply(fpkm[all_genes,M14h_E2C],1,mean),apply(fpkm[all_genes,M14h_L2C],1,mean),apply(fpkm[all_genes,M14h_4C],1,mean))
colnames(M14h_avgfpkm) <- c("M14h_CC","M14h_1C","M14h_E2C","M14h_L2C","M14h_4C");rownames(Ctrl_avgfpkm) <- all_genes
Dev_avgfpkm <- cbind(apply(Dev_fpkm[all_genes,Dev_MII],1,mean),apply(Dev_fpkm[all_genes,Dev_1C],1,mean),apply(Dev_fpkm[all_genes,Dev_E2C],1,mean),apply(Dev_fpkm[all_genes,Dev_L2C],1,mean),apply(Dev_fpkm[all_genes,Dev_4C],1,mean))
colnames(Dev_avgfpkm) <- c("Dev_MII","Dev_1C","Dev_E2C","Dev_L2C","Dev_4C");rownames(Ctrl_avgfpkm) <- all_genes

###################################################################################################
################################                PLOT               ################################
###################################################################################################

#### Falsely Activated/Suppressed during Ctrl SCNT vs Development 
### 1C 
log2fc <- log2(Ctrl_avgfpkm[all_genes,"Ctrl_1C"]/Dev_avgfpkm[all_genes,"Dev_1C"])
mlog10pvalue <- c()
for (each in all_genes){
	tmp_x <- table(as.numeric(fpkm[each,Ctrl_1C14hpa]))
	tmp_y <- table(as.numeric(Dev_fpkm[each,Dev_1C]))
	if (length(tmp_x)==1 & length(tmp_y)==1 ){
		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
	}else{
	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,Ctrl_1C14hpa],Dev_fpkm[each,Dev_1C])$p.value))
	}
}
mlog10pvalue[mlog10pvalue > 6] <- 6
mlog10pvalue[mlog10pvalue == -1] <- 6
log2fc[log2fc > 10] <- 10
log2fc[log2fc < -10] <- -10
names(mlog10pvalue) <- all_genes
names(log2fc) <- all_genes
Ctrl_1C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
Ctrl_1C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]

# markers <- c()
pdf("FalselyActSup.Ctrl_1C.VolcanoPlot.pdf",width=4,height=4.5)
par(mar=c(6,4,4,2))
plot(log2fc,mlog10pvalue,pch=".",main="",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="",ylab="",bty="n",xaxt="n",yaxt="n");
# plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_1C/Dev_1C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
abline(h=-log10(0.01), lty=2)
abline(v=c(-0.5,0.5), lty=2)
points(log2fc[Ctrl_1C_dn_genes],mlog10pvalue[Ctrl_1C_dn_genes],pch=16,col="#235E7C",cex=0.3)
points(log2fc[Ctrl_1C_up_genes],mlog10pvalue[Ctrl_1C_up_genes],pch=16,col="#A40126",cex=0.3)
# points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(Ctrl_1C_dn_genes),length(Ctrl_1C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
dev.off()
write.table(Ctrl_1C_dn_genes,file="Ctrl_1C_dn_genes.txt",col.names=F,row.names=F,quote=F)
write.table(Ctrl_1C_up_genes,file="Ctrl_1C_up_genes.txt",col.names=F,row.names=F,quote=F)

### E2C 
log2fc <- log2(Ctrl_avgfpkm[all_genes,"Ctrl_E2C"]/Dev_avgfpkm[all_genes,"Dev_E2C"])
mlog10pvalue <- c()
for (each in all_genes){
	tmp_x <- table(as.numeric(fpkm[each,Ctrl_E2C]))
	tmp_y <- table(as.numeric(Dev_fpkm[each,Dev_E2C]))
	if (length(tmp_x)==1 & length(tmp_y)==1 ){
		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
	}else{
	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,Ctrl_E2C],Dev_fpkm[each,Dev_E2C])$p.value))
	}
}
mlog10pvalue[mlog10pvalue > 6] <- 6
mlog10pvalue[mlog10pvalue == -1] <- 6
log2fc[log2fc > 10] <- 10
log2fc[log2fc < -10] <- -10
names(mlog10pvalue) <- all_genes
names(log2fc) <- all_genes
Ctrl_E2C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
Ctrl_E2C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]

# markers <- c()
pdf("FalselyActSup.Ctrl_E2C.VolcanoPlot.pdf",width=4,height=4.5)
par(mar=c(6,4,4,2))
plot(log2fc,mlog10pvalue,pch=".",main="",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="",ylab="",bty="n",xaxt="n",yaxt="n");
# plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
abline(h=-log10(0.01), lty=2)
abline(v=c(-0.5,0.5), lty=2)
points(log2fc[Ctrl_E2C_dn_genes],mlog10pvalue[Ctrl_E2C_dn_genes],pch=16,col="#235E7C",cex=0.3)
points(log2fc[Ctrl_E2C_up_genes],mlog10pvalue[Ctrl_E2C_up_genes],pch=16,col="#A40126",cex=0.3)
# points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(Ctrl_E2C_dn_genes),length(Ctrl_E2C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
dev.off()
write.table(Ctrl_E2C_dn_genes,file="Ctrl_E2C_dn_genes.txt",col.names=F,row.names=F,quote=F)
write.table(Ctrl_E2C_up_genes,file="Ctrl_E2C_up_genes.txt",col.names=F,row.names=F,quote=F)

### L2C 
log2fc <- log2(Ctrl_avgfpkm[all_genes,"Ctrl_L2C"]/Dev_avgfpkm[all_genes,"Dev_L2C"])
mlog10pvalue <- c()
for (each in all_genes){
	tmp_x <- table(as.numeric(fpkm[each,Ctrl_L2C]))
	tmp_y <- table(as.numeric(Dev_fpkm[each,Dev_L2C]))
	if (length(tmp_x)==1 & length(tmp_y)==1 ){
		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
	}else{
	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,Ctrl_L2C],Dev_fpkm[each,Dev_L2C])$p.value))
	}
}
mlog10pvalue[mlog10pvalue > 6] <- 6
mlog10pvalue[mlog10pvalue == -1] <- 6
log2fc[log2fc > 10] <- 10
log2fc[log2fc < -10] <- -10
names(mlog10pvalue) <- all_genes
names(log2fc) <- all_genes
Ctrl_L2C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
Ctrl_L2C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]

# markers <- c()
pdf("FalselyActSup.Ctrl_L2C.VolcanoPlot.pdf",width=4,height=4.5)
par(mar=c(6,4,4,2))
plot(log2fc,mlog10pvalue,pch=".",main="",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="",ylab="",bty="n",xaxt="n",yaxt="n");
# plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_L2C/Dev_L2C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
abline(h=-log10(0.01), lty=2)
abline(v=c(-0.5,0.5), lty=2)
points(log2fc[Ctrl_L2C_dn_genes],mlog10pvalue[Ctrl_L2C_dn_genes],pch=16,col="#235E7C",cex=0.3)
points(log2fc[Ctrl_L2C_up_genes],mlog10pvalue[Ctrl_L2C_up_genes],pch=16,col="#A40126",cex=0.3)
# points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(Ctrl_L2C_dn_genes),length(Ctrl_L2C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
dev.off()
write.table(Ctrl_L2C_dn_genes,file="Ctrl_L2C_dn_genes.txt",col.names=F,row.names=F,quote=F)
write.table(Ctrl_L2C_up_genes,file="Ctrl_L2C_up_genes.txt",col.names=F,row.names=F,quote=F)


### 4C 
log2fc <- log2(Ctrl_avgfpkm[all_genes,"Ctrl_4C"]/Dev_avgfpkm[all_genes,"Dev_4C"])
mlog10pvalue <- c()
for (each in all_genes){
	tmp_x <- table(as.numeric(fpkm[each,Ctrl_4C]))
	tmp_y <- table(as.numeric(Dev_fpkm[each,Dev_4C]))
	if (length(tmp_x)==1 & length(tmp_y)==1 ){
		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
	}else{
	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,Ctrl_4C],Dev_fpkm[each,Dev_4C])$p.value))
	}
}
mlog10pvalue[mlog10pvalue > 6] <- 6
mlog10pvalue[mlog10pvalue == -1] <- 6
log2fc[log2fc > 10] <- 10
log2fc[log2fc < -10] <- -10
names(mlog10pvalue) <- all_genes
names(log2fc) <- all_genes
Ctrl_4C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
Ctrl_4C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]

# markers <- c()
pdf("FalselyActSup.Ctrl_4C.VolcanoPlot.pdf",width=4,height=4.5)
par(mar=c(6,4,4,2))
plot(log2fc,mlog10pvalue,pch=".",main="",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="",ylab="",bty="n",xaxt="n",yaxt="n");
# plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_4C/Dev_4C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
abline(h=-log10(0.01), lty=2)
abline(v=c(-0.5,0.5), lty=2)
points(log2fc[Ctrl_4C_dn_genes],mlog10pvalue[Ctrl_4C_dn_genes],pch=16,col="#235E7C",cex=0.3)
points(log2fc[Ctrl_4C_up_genes],mlog10pvalue[Ctrl_4C_up_genes],pch=16,col="#A40126",cex=0.3)
# points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(Ctrl_4C_dn_genes),length(Ctrl_4C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
dev.off()
write.table(Ctrl_4C_dn_genes,file="Ctrl_4C_dn_genes.txt",col.names=F,row.names=F,quote=F)
write.table(Ctrl_4C_up_genes,file="Ctrl_4C_up_genes.txt",col.names=F,row.names=F,quote=F)



#### Falsely Activated/Suppressed during M14h SCNT vs Development 
### 1C 
log2fc <- log2(M14h_avgfpkm[all_genes,"M14h_1C"]/Dev_avgfpkm[all_genes,"Dev_1C"])
mlog10pvalue <- c()
for (each in all_genes){
	tmp_x <- table(as.numeric(fpkm[each,M14h_1C14hpa]))
	tmp_y <- table(as.numeric(Dev_fpkm[each,Dev_1C]))
	if (length(tmp_x)==1 & length(tmp_y)==1 ){
		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
	}else{
	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,M14h_1C14hpa],Dev_fpkm[each,Dev_1C])$p.value))
	}
}
mlog10pvalue[mlog10pvalue > 6] <- 6
mlog10pvalue[mlog10pvalue == -1] <- 6
log2fc[log2fc > 10] <- 10
log2fc[log2fc < -10] <- -10
names(mlog10pvalue) <- all_genes
names(log2fc) <- all_genes
M14h_1C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
M14h_1C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]
M14h_1C_nd_genes <- all_genes[which(log2fc <= 0.5 & log2fc >= -0.5 & mlog10pvalue < -log10(0.01))]

# markers <- c()
pdf("FalselyActSup.M14h_1C.VolcanoPlot.pdf",width=4,height=4.5)
par(mar=c(6,4,4,2))
plot(log2fc,mlog10pvalue,pch=".",main="",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="",ylab="",bty="n",xaxt="n",yaxt="n");
# plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(M14h_1C/Dev_1C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
abline(h=-log10(0.01), lty=2)
abline(v=c(-0.5,0.5), lty=2)
points(log2fc[M14h_1C_dn_genes],mlog10pvalue[M14h_1C_dn_genes],pch=16,col="#235E7C",cex=0.3)
points(log2fc[M14h_1C_up_genes],mlog10pvalue[M14h_1C_up_genes],pch=16,col="#A40126",cex=0.3)
# points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(M14h_1C_dn_genes),length(M14h_1C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
dev.off()
write.table(M14h_1C_dn_genes,file="M14h_1C_dn_genes.txt",col.names=F,row.names=F,quote=F)
write.table(M14h_1C_up_genes,file="M14h_1C_up_genes.txt",col.names=F,row.names=F,quote=F)
write.table(M14h_1C_nd_genes,file="M14h_1C_nd_genes.txt",col.names=F,row.names=F,quote=F)

### E2C 
log2fc <- log2(M14h_avgfpkm[all_genes,"M14h_E2C"]/Dev_avgfpkm[all_genes,"Dev_E2C"])
mlog10pvalue <- c()
for (each in all_genes){
	tmp_x <- table(as.numeric(fpkm[each,M14h_E2C]))
	tmp_y <- table(as.numeric(Dev_fpkm[each,Dev_E2C]))
	if (length(tmp_x)==1 & length(tmp_y)==1 ){
		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
	}else{
	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,M14h_E2C],Dev_fpkm[each,Dev_E2C])$p.value))
	}
}
mlog10pvalue[mlog10pvalue > 6] <- 6
mlog10pvalue[mlog10pvalue == -1] <- 6
log2fc[log2fc > 10] <- 10
log2fc[log2fc < -10] <- -10
names(mlog10pvalue) <- all_genes
names(log2fc) <- all_genes
M14h_E2C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
M14h_E2C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]
M14h_E2C_nd_genes <- all_genes[which(log2fc <= 0.5 & log2fc >= -0.5 & mlog10pvalue < -log10(0.01))]

# markers <- c()
pdf("FalselyActSup.M14h_E2C.VolcanoPlot.pdf",width=4,height=4.5)
par(mar=c(6,4,4,2))
plot(log2fc,mlog10pvalue,pch=".",main="",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="",ylab="",bty="n",xaxt="n",yaxt="n");
# plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(M14h_E2C/Dev_E2C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
abline(h=-log10(0.01), lty=2)
abline(v=c(-0.5,0.5), lty=2)
points(log2fc[M14h_E2C_dn_genes],mlog10pvalue[M14h_E2C_dn_genes],pch=16,col="#235E7C",cex=0.3)
points(log2fc[M14h_E2C_up_genes],mlog10pvalue[M14h_E2C_up_genes],pch=16,col="#A40126",cex=0.3)
# points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(M14h_E2C_dn_genes),length(M14h_E2C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
dev.off()
write.table(M14h_E2C_dn_genes,file="M14h_E2C_dn_genes.txt",col.names=F,row.names=F,quote=F)
write.table(M14h_E2C_up_genes,file="M14h_E2C_up_genes.txt",col.names=F,row.names=F,quote=F)
write.table(M14h_E2C_nd_genes,file="M14h_E2C_nd_genes.txt",col.names=F,row.names=F,quote=F)

### L2C 
log2fc <- log2(M14h_avgfpkm[all_genes,"M14h_L2C"]/Dev_avgfpkm[all_genes,"Dev_L2C"])
mlog10pvalue <- c()
for (each in all_genes){
	tmp_x <- table(as.numeric(fpkm[each,M14h_L2C]))
	tmp_y <- table(as.numeric(Dev_fpkm[each,Dev_L2C]))
	if (length(tmp_x)==1 & length(tmp_y)==1 ){
		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
	}else{
	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,M14h_L2C],Dev_fpkm[each,Dev_L2C])$p.value))
	}
}
mlog10pvalue[mlog10pvalue > 6] <- 6
mlog10pvalue[mlog10pvalue == -1] <- 6
log2fc[log2fc > 10] <- 10
log2fc[log2fc < -10] <- -10
names(mlog10pvalue) <- all_genes
names(log2fc) <- all_genes
M14h_L2C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
M14h_L2C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]
M14h_L2C_nd_genes <- all_genes[which(log2fc <= 0.5 & log2fc >= -0.5 & mlog10pvalue < -log10(0.01))]

# markers <- c()
pdf("FalselyActSup.M14h_L2C.VolcanoPlot.pdf",width=4,height=4.5)
par(mar=c(6,4,4,2))
plot(log2fc,mlog10pvalue,pch=".",main="",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="",ylab="",bty="n",xaxt="n",yaxt="n");
# plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(M14h_L2C/Dev_L2C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
abline(h=-log10(0.01), lty=2)
abline(v=c(-0.5,0.5), lty=2)
points(log2fc[M14h_L2C_dn_genes],mlog10pvalue[M14h_L2C_dn_genes],pch=16,col="#235E7C",cex=0.3)
points(log2fc[M14h_L2C_up_genes],mlog10pvalue[M14h_L2C_up_genes],pch=16,col="#A40126",cex=0.3)
# points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(M14h_L2C_dn_genes),length(M14h_L2C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
dev.off()
write.table(M14h_L2C_dn_genes,file="M14h_L2C_dn_genes.txt",col.names=F,row.names=F,quote=F)
write.table(M14h_L2C_up_genes,file="M14h_L2C_up_genes.txt",col.names=F,row.names=F,quote=F)
write.table(M14h_L2C_nd_genes,file="M14h_L2C_nd_genes.txt",col.names=F,row.names=F,quote=F)


### 4C 
log2fc <- log2(M14h_avgfpkm[all_genes,"M14h_4C"]/Dev_avgfpkm[all_genes,"Dev_4C"])
mlog10pvalue <- c()
for (each in all_genes){
	tmp_x <- table(as.numeric(fpkm[each,M14h_4C]))
	tmp_y <- table(as.numeric(Dev_fpkm[each,Dev_4C]))
	if (length(tmp_x)==1 & length(tmp_y)==1 ){
		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
	}else{
	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,M14h_4C],Dev_fpkm[each,Dev_4C])$p.value))
	}
}
mlog10pvalue[mlog10pvalue > 6] <- 6
mlog10pvalue[mlog10pvalue == -1] <- 6
log2fc[log2fc > 10] <- 10
log2fc[log2fc < -10] <- -10
names(mlog10pvalue) <- all_genes
names(log2fc) <- all_genes
M14h_4C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
M14h_4C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]
M14h_4C_nd_genes <- all_genes[which(log2fc <= 0.5 & log2fc >= -0.5 & mlog10pvalue < -log10(0.01))]

# markers <- c()
pdf("FalselyActSup.M14h_4C.VolcanoPlot.pdf",width=4,height=4.5)
par(mar=c(6,4,4,2))
plot(log2fc,mlog10pvalue,pch=".",main="",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="",ylab="",bty="n",xaxt="n",yaxt="n");
# plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(M14h_4C/Dev_4C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
abline(h=-log10(0.01), lty=2)
abline(v=c(-0.5,0.5), lty=2)
points(log2fc[M14h_4C_dn_genes],mlog10pvalue[M14h_4C_dn_genes],pch=16,col="#235E7C",cex=0.3)
points(log2fc[M14h_4C_up_genes],mlog10pvalue[M14h_4C_up_genes],pch=16,col="#A40126",cex=0.3)
# points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(M14h_4C_dn_genes),length(M14h_4C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
dev.off()
write.table(M14h_4C_dn_genes,file="M14h_4C_dn_genes.txt",col.names=F,row.names=F,quote=F)
write.table(M14h_4C_up_genes,file="M14h_4C_up_genes.txt",col.names=F,row.names=F,quote=F)
write.table(M14h_4C_nd_genes,file="M14h_4C_nd_genes.txt",col.names=F,row.names=F,quote=F)






# #### Falsely Activated/Suppressed during M14h SCNT vs Ctrl SCNT 
# ### 1C 
# log2fc <- log2(M14h_avgfpkm[all_genes,"M14h_1C"]/Ctrl_avgfpkm[all_genes,"Ctrl_1C"])
# mlog10pvalue <- c()
# for (each in all_genes){
# 	tmp_x <- table(as.numeric(fpkm[each,M14h_1C14hpa]))
# 	tmp_y <- table(as.numeric(fpkm[each,Ctrl_1C14hpa]))
# 	if (length(tmp_x)==1 & length(tmp_y)==1 ){
# 		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
# 	}else{
# 	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,M14h_1C14hpa],fpkm[each,Ctrl_1C14hpa])$p.value))
# 	}
# }
# mlog10pvalue[mlog10pvalue > 6] <- 6
# mlog10pvalue[mlog10pvalue == -1] <- 6
# log2fc[log2fc > 10] <- 10
# log2fc[log2fc < -10] <- -10
# names(mlog10pvalue) <- all_genes
# names(log2fc) <- all_genes
# M14h_1C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
# M14h_1C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]

# # markers <- c()
# pdf("FalselyActSup.M14hCtrl_1C.VolcanoPlot.pdf",width=4,height=4.5)
# par(mar=c(6,4,4,2))
# plot(log2fc,mlog10pvalue,pch=".", change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD"ylab,bty="n",xaxt="n",yaxt="n"="");box(lwd=1)
plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(M14h_1C/Ctrl_1C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
# abline(h=-log10(0.01), lty=2)
# abline(v=c(-0.5,0.5), lty=2)
# points(log2fc[M14h_1C_dn_genes],mlog10pvalue[M14h_1C_dn_genes],pch=16,col="#235E7C",cex=0.3)
# points(log2fc[M14h_1C_up_genes],mlog10pvalue[M14h_1C_up_genes],pch=16,col="#A40126",cex=0.3)
# # points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# # text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
# text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(M14h_1C_dn_genes),length(M14h_1C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
# dev.off()
# write.table(M14h_1C_dn_genes,file="M14hCtrl_1C_dn_genes.txt",col.names=F,row.names=F,quote=F)
# write.table(M14h_1C_up_genes,file="M14hCtrl_1C_up_genes.txt",col.names=F,row.names=F,quote=F)

# ### E2C 
# log2fc <- log2(M14h_avgfpkm[all_genes,"M14h_E2C"]/Ctrl_avgfpkm[all_genes,"Ctrl_E2C"])
# mlog10pvalue <- c()
# for (each in all_genes){
# 	tmp_x <- table(as.numeric(fpkm[each,M14h_E2C]))
# 	tmp_y <- table(as.numeric(fpkm[each,Ctrl_E2C]))
# 	if (length(tmp_x)==1 & length(tmp_y)==1 ){
# 		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
# 	}else{
# 	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,M14h_E2C],fpkm[each,Ctrl_E2C])$p.value))
# 	}
# }
# mlog10pvalue[mlog10pvalue > 6] <- 6
# mlog10pvalue[mlog10pvalue == -1] <- 6
# log2fc[log2fc > 10] <- 10
# log2fc[log2fc < -10] <- -10
# names(mlog10pvalue) <- all_genes
# names(log2fc) <- all_genes
# M14h_E2C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
# M14h_E2C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]

# # markers <- c()
# pdf("FalselyActSup.M14hCtrl_E2C.VolcanoPlot.pdf",width=4,height=4.5)
# par(mar=c(6,4,4,2))
# plot(log2fc,mlog10pvalue,pch=".", change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD"ylab,bty="n",xaxt="n",yaxt="n"="");box(lwd=1)
plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(M14h_E2C/Ctrl_E2C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
# abline(h=-log10(0.01), lty=2)
# abline(v=c(-0.5,0.5), lty=2)
# points(log2fc[M14h_E2C_dn_genes],mlog10pvalue[M14h_E2C_dn_genes],pch=16,col="#235E7C",cex=0.3)
# points(log2fc[M14h_E2C_up_genes],mlog10pvalue[M14h_E2C_up_genes],pch=16,col="#A40126",cex=0.3)
# # points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# # text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
# text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(M14h_E2C_dn_genes),length(M14h_E2C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
# dev.off()
# write.table(M14h_E2C_dn_genes,file="M14hCtrl_E2C_dn_genes.txt",col.names=F,row.names=F,quote=F)
# write.table(M14h_E2C_up_genes,file="M14hCtrl_E2C_up_genes.txt",col.names=F,row.names=F,quote=F)

# ### L2C 
# log2fc <- log2(M14h_avgfpkm[all_genes,"M14h_L2C"]/Ctrl_avgfpkm[all_genes,"Ctrl_L2C"])
# mlog10pvalue <- c()
# for (each in all_genes){
# 	tmp_x <- table(as.numeric(fpkm[each,M14h_L2C]))
# 	tmp_y <- table(as.numeric(fpkm[each,Ctrl_L2C]))
# 	if (length(tmp_x)==1 & length(tmp_y)==1 ){
# 		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
# 	}else{
# 	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,M14h_L2C],fpkm[each,Ctrl_L2C])$p.value))
# 	}
# }
# mlog10pvalue[mlog10pvalue > 6] <- 6
# mlog10pvalue[mlog10pvalue == -1] <- 6
# log2fc[log2fc > 10] <- 10
# log2fc[log2fc < -10] <- -10
# names(mlog10pvalue) <- all_genes
# names(log2fc) <- all_genes
# M14h_L2C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
# M14h_L2C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]

# # markers <- c()
# pdf("FalselyActSup.M14hCtrl_L2C.VolcanoPlot.pdf",width=4,height=4.5)
# par(mar=c(6,4,4,2))
# plot(log2fc,mlog10pvalue,pch=".", change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD"ylab,bty="n",xaxt="n",yaxt="n"="");box(lwd=1)
plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(M14h_L2C/Ctrl_L2C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
# abline(h=-log10(0.01), lty=2)
# abline(v=c(-0.5,0.5), lty=2)
# points(log2fc[M14h_L2C_dn_genes],mlog10pvalue[M14h_L2C_dn_genes],pch=16,col="#235E7C",cex=0.3)
# points(log2fc[M14h_L2C_up_genes],mlog10pvalue[M14h_L2C_up_genes],pch=16,col="#A40126",cex=0.3)
# # points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# # text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
# text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(M14h_L2C_dn_genes),length(M14h_L2C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
# dev.off()
# write.table(M14h_L2C_dn_genes,file="M14hCtrl_L2C_dn_genes.txt",col.names=F,row.names=F,quote=F)
# write.table(M14h_L2C_up_genes,file="M14hCtrl_L2C_up_genes.txt",col.names=F,row.names=F,quote=F)


# ### 4C 
# log2fc <- log2(M14h_avgfpkm[all_genes,"M14h_4C"]/Ctrl_avgfpkm[all_genes,"Ctrl_4C"])
# mlog10pvalue <- c()
# for (each in all_genes){
# 	tmp_x <- table(as.numeric(fpkm[each,M14h_4C]))
# 	tmp_y <- table(as.numeric(fpkm[each,Ctrl_4C]))
# 	if (length(tmp_x)==1 & length(tmp_y)==1 ){
# 		if (tmp_x!=tmp_y){mlog10pvalue <- c(mlog10pvalue, -1)} else {mlog10pvalue <- c(mlog10pvalue, 0)}
# 	}else{
# 	mlog10pvalue <- c(mlog10pvalue, -log10(t.test(fpkm[each,M14h_4C],fpkm[each,Ctrl_4C])$p.value))
# 	}
# }
# mlog10pvalue[mlog10pvalue > 6] <- 6
# mlog10pvalue[mlog10pvalue == -1] <- 6
# log2fc[log2fc > 10] <- 10
# log2fc[log2fc < -10] <- -10
# names(mlog10pvalue) <- all_genes
# names(log2fc) <- all_genes
# M14h_4C_up_genes <- all_genes[which(log2fc > 0.5 & mlog10pvalue >= -log10(0.01))]
# M14h_4C_dn_genes <- all_genes[which(log2fc < -0.5 & mlog10pvalue >= -log10(0.01))]

# # markers <- c()
# pdf("FalselyActSup.M14hCtrl_4C.VolcanoPlot.pdf",width=4,height=4.5)
# par(mar=c(6,4,4,2))
# plot(log2fc,mlog10pvalue,pch=".", change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD"ylab,bty="n",xaxt="n",yaxt="n"="");box(lwd=1)
plot(log2fc,mlog10pvalue,pch=".",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(M14h_4C/Ctrl_4C)",ylab="-log10(p-value)");box(lwd=1)
plot(1,type="n",main="Expression change",xlim=c(-10,10),ylim=c(0,6),col="#ADADAD",xlab="log2(Ctrl_E2C/Dev_E2C)",ylab="-log10(p-value)")
# abline(h=-log10(0.01), lty=2)
# abline(v=c(-0.5,0.5), lty=2)
# points(log2fc[M14h_4C_dn_genes],mlog10pvalue[M14h_4C_dn_genes],pch=16,col="#235E7C",cex=0.3)
# points(log2fc[M14h_4C_up_genes],mlog10pvalue[M14h_4C_up_genes],pch=16,col="#A40126",cex=0.3)
# # points(log2fc[markers],mlog10pvalue[markers],pch="*",col=cccol[4],cex=1)
# # text(log2fc[markers],mlog10pvalue[markers]+0.2,labels=markers,col=cccol[4],cex=0.5)
# text(c(-8,8),c(5.8,5.8),labels=paste("n = ",c(length(M14h_4C_dn_genes),length(M14h_4C_up_genes)),sep=""),col=c("#235E7C","#A40126"),cex=0.8)
# dev.off()
# write.table(M14h_4C_dn_genes,file="M14hCtrl_4C_dn_genes.txt",col.names=F,row.names=F,quote=F)
# write.table(M14h_4C_up_genes,file="M14hCtrl_4C_up_genes.txt",col.names=F,row.names=F,quote=F)

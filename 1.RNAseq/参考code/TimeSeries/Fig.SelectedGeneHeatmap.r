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

############ Selected genes
Totipotent <- intersect(all_genes, c("Zscan4c","Lmx1a","Gm21319","Usp17lc","Dcc","Ttc30b","Olfr815","Zfp352","Gm8994","Fgf1"))
Pluripotent <- intersect(all_genes, c("Esrrb","Rest","Dppa4","Tcf15","Tet2","Bcl3","Fgf4","Sfrp1","Pou2f3","Pou5f1"))
Translation <- intersect(all_genes, c("Rxra","Dnajb14","Uchl1os","Bcl3","Patl2","Thbs1","Hspa1a","Eif4e1b","Hspa1b","Gch1","Dazl","Per2","Cpeb4","Pabpc6","Ago4","App","Rbm24","Mrpl53","Tia1","Krt17","Gadd45gip1","Crtap","Eif3i","Elavl4","Eif2s3y","Pde12"))
Apoptosis <- intersect(all_genes, c("Plscr2","Xkr9","Naip1","Gper1","Thap2","Naip6","Xkr5","Cidec","Api5","Naip3l1","cd24a"))
HistoneModification <- intersect(all_genes, c("Dpy30","Prmt7","Bmp4","Smyd3","Tssk6","Uty","Kdm5c","Tbl1xr1","Rest","Phf20","Kansl1","Carm1","Sap30","H3f3a","H2bc1","Kmt2b","Gps2","Jade3","Lin28a"))
NucleosomeAssembly <- intersect(all_genes, c("H1f1","H1f4","H2ac10","H2ac11","H2ac12","H2ac13","H2ac15","H2ac4","H2ac7","H2ac8","H2bc21","H3c1","H3c10","H3c2","H4c14","H4c2","Macroh2a2"))
ChromatinRemodeling <- intersect(all_genes, c("Cenpv","Ino80","Wapl","Smarcal1"))

pu_C_up_specific_TranslationGenes_tmp <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.DiffPattern/RNAseq/GeneFilter_1/gfold/pu_C_up_specific_TranslationGenes.txt",header=F,row.names=1)
pu_C_up_specific_TranslationGenes <- rownames(pu_C_up_specific_TranslationGenes_tmp)[order(pu_C_up_specific_TranslationGenes_tmp[,2])]
pu_C_up_specific_ApoptosisGenes_tmp <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.DiffPattern/RNAseq/GeneFilter_1/gfold/pu_C_up_specific_ApoptosisGenes.txt",header=F,row.names=1)
pu_C_up_specific_ApoptosisGenes <- rownames(pu_C_up_specific_ApoptosisGenes_tmp)[order(pu_C_up_specific_ApoptosisGenes_tmp[,2])]
pu_C_dn_specific_TranslationGenes_tmp <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.DiffPattern/RNAseq/GeneFilter_1/gfold/pu_C_dn_specific_TranslationGenes.txt",header=F,row.names=1)
pu_C_dn_specific_TranslationGenes <- rownames(pu_C_dn_specific_TranslationGenes_tmp)[order(pu_C_dn_specific_TranslationGenes_tmp[,2])]
pu_C_dn_specific_ApoptosisGenes_tmp <- read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.DiffPattern/RNAseq/GeneFilter_1/gfold/pu_C_dn_specific_ApoptosisGenes.txt",header=F,row.names=1)
pu_C_dn_specific_ApoptosisGenes <- rownames(pu_C_dn_specific_ApoptosisGenes_tmp)[order(pu_C_dn_specific_ApoptosisGenes_tmp[,2])]

pu_C_up_specific <- intersect(all_genes, as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.DiffPattern/RNAseq/GeneFilter_2/pu_C_up_specific.txt")[,1]))
pu_C_dn_specific <- intersect(all_genes, as.vector(read.table("/data/zhaochengchen/Work/OTHER/01.WQN/2.WQN_DUXno_mmESCs_RNAseq/result/2.DiffPattern/RNAseq/GeneFilter_2/pu_C_dn_specific.txt")[,1]))

###################################################################################################
################################                PLOT               ################################
###################################################################################################

# Totipotent
pdf("Fig.SelectedGeneHeatmap.Totipotent.pdf",width=8,height=length(Totipotent)/5+1.2)
plot_data_mESC <- logfpkm[Totipotent,mESC_path]
plot_data_tCLC <- logfpkm[Totipotent,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=Totipotent,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=Totipotent,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()


# Pluripotent
pdf("Fig.SelectedGeneHeatmap.Pluripotent.pdf",width=8,height=length(Pluripotent)/5+1.2)
plot_data_mESC <- logfpkm[Pluripotent,mESC_path]
plot_data_tCLC <- logfpkm[Pluripotent,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=Pluripotent,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=Pluripotent,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()


# Translation
pdf("Fig.SelectedGeneHeatmap.Translation.pdf",width=8,height=length(Translation)/5+1.2)
plot_data_mESC <- logfpkm[Translation,mESC_path]
plot_data_tCLC <- logfpkm[Translation,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=Translation,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=Translation,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()

# Apoptosis
pdf("Fig.SelectedGeneHeatmap.Apoptosis.pdf",width=8,height=length(Apoptosis)/5+1.2)
plot_data_mESC <- logfpkm[Apoptosis,mESC_path]
plot_data_tCLC <- logfpkm[Apoptosis,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=Apoptosis,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=Apoptosis,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()



# HistoneModification
pdf("Fig.SelectedGeneHeatmap.HistoneModification.pdf",width=8,height=length(HistoneModification)/5+1.2)
plot_data_mESC <- logfpkm[HistoneModification,mESC_path]
plot_data_tCLC <- logfpkm[HistoneModification,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=HistoneModification,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=HistoneModification,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()


# NucleosomeAssembly
pdf("Fig.SelectedGeneHeatmap.NucleosomeAssembly.pdf",width=8,height=length(NucleosomeAssembly)/5+1.2)
plot_data_mESC <- logfpkm[NucleosomeAssembly,mESC_path]
plot_data_tCLC <- logfpkm[NucleosomeAssembly,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=NucleosomeAssembly,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=NucleosomeAssembly,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()


# ChromatinRemodeling
pdf("Fig.SelectedGeneHeatmap.ChromatinRemodeling.pdf",width=8,height=length(ChromatinRemodeling)/5+1.2)
plot_data_mESC <- logfpkm[ChromatinRemodeling,mESC_path]
plot_data_tCLC <- logfpkm[ChromatinRemodeling,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=ChromatinRemodeling,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=ChromatinRemodeling,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()




###############################################
# pu_C_up_specific_TranslationGenes
pdf("Fig.SelectedGeneHeatmap.pu_C_up_specific_TranslationGenes.pdf",width=8,height=length(pu_C_up_specific_TranslationGenes)/5+1.2)
plot_data_mESC <- logfpkm[pu_C_up_specific_TranslationGenes,mESC_path]
plot_data_tCLC <- logfpkm[pu_C_up_specific_TranslationGenes,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_up_specific_TranslationGenes,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_up_specific_TranslationGenes,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()

# pu_C_up_specific_ApoptosisGenes
pdf("Fig.SelectedGeneHeatmap.pu_C_up_specific_ApoptosisGenes.pdf",width=8,height=length(pu_C_up_specific_ApoptosisGenes)/5+1.2)
plot_data_mESC <- logfpkm[pu_C_up_specific_ApoptosisGenes,mESC_path]
plot_data_tCLC <- logfpkm[pu_C_up_specific_ApoptosisGenes,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_up_specific_ApoptosisGenes,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_up_specific_ApoptosisGenes,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()

# pu_C_dn_specific_TranslationGenes
pdf("Fig.SelectedGeneHeatmap.pu_C_dn_specific_TranslationGenes.pdf",width=8,height=length(pu_C_dn_specific_TranslationGenes)/5+1.2)
plot_data_mESC <- logfpkm[pu_C_dn_specific_TranslationGenes,mESC_path]
plot_data_tCLC <- logfpkm[pu_C_dn_specific_TranslationGenes,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_dn_specific_TranslationGenes,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_dn_specific_TranslationGenes,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()

# pu_C_dn_specific_ApoptosisGenes
pdf("Fig.SelectedGeneHeatmap.pu_C_dn_specific_ApoptosisGenes.pdf",width=8,height=length(pu_C_dn_specific_ApoptosisGenes)/5+1.2)
plot_data_mESC <- logfpkm[pu_C_dn_specific_ApoptosisGenes,mESC_path]
plot_data_tCLC <- logfpkm[pu_C_dn_specific_ApoptosisGenes,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_dn_specific_ApoptosisGenes,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_dn_specific_ApoptosisGenes,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()


# pu_C_up_specific
pdf("Fig.SelectedGeneHeatmap.pu_C_up_specific.pdf",width=8,height=length(pu_C_up_specific)/5+1.2)
plot_data_mESC <- logfpkm[pu_C_up_specific,mESC_path]
plot_data_tCLC <- logfpkm[pu_C_up_specific,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_up_specific,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_up_specific,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()


# pu_C_dn_specific
pdf("Fig.SelectedGeneHeatmap.pu_C_dn_specific.pdf",width=8,height=length(pu_C_dn_specific)/5+1.2)
plot_data_mESC <- logfpkm[pu_C_dn_specific,mESC_path]
plot_data_tCLC <- logfpkm[pu_C_dn_specific,tCLC_path]
layout(matrix(c(rep(1,8),rep(2,8),rep(3,2)),nrow=1,ncol=18))
all_exp <- c(as.matrix(plot_data_mESC),as.matrix(plot_data_tCLC)) # using same scale bar
zmax <- quantile(all_exp,0.99,na.rm=T)
zmin <- quantile(all_exp,0.01,na.rm=T)

ColorRamp <- colorRampPalette(c("lightblue","white","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data_mESC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="mESC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_dn_specific,las=2)

tmp_matrix <- plot_data_tCLC
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(4,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp, xlab="", main="2CLC",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=timepoints,las=2)
axis(side=2,1:nrow(tmp_matrix),labels=pu_C_dn_specific,las=2)

par(mar=c(4,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()
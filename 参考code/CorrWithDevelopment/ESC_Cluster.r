###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
allcol <- c("#f7acbc","#deab8a","#817936","#444693","#ef5b9c","#fedcbd","#7f7522","#2b4490","#f47920","#80752c","#2a5caa","#f05b72","#905a3d","#87843b","#224b8f","#f15b6c","#8f4b2e","#726930","#003a6c","#f8aba6","#87481f","#454926","#102b6a","#f69c9f","#5f3c23","#2e3a1f","#426ab3","#f58f98","#6b473c","#4d4f36","#46485f","#ca8687","#faa755","#b7ba6b","#4e72b8","#f391a9","#fab27b","#b2d235","#181d4b","#bd6758","#f58220","#5c7a29","#1a2933","#d71345","#843900","#bed742","#121a2a","#d64f44","#905d1d","#7fb80e","#0c212b","#d93a49","#8a5d19","#a3cf62","#6a6da9","#b3424a","#8c531b","#769149","#585eaa","#c76968","#826858","#6d8346","#494e8f","#bb505d","#64492b","#78a355","#afb4db","#987165","#ae6642","#abc88b","#9b95c9","#ac6767","#56452d","#74905d","#6950a1","#973c3f","#96582a","#cde6c7","#6f60aa","#b22c46","#705628","#1d953f","#867892","#a7324a","#4a3113","#77ac98","#918597","#aa363d","#412f1f","#007d65","#6f6d85","#ed1941","#845538","#84bf96","#594c6d","#f26522","#8e7437","#45b97c","#694d9f","#d2553d","#69541b","#225a1f","#6f599c","#b4534b","#d5c59f","#367459","#8552a1","#ef4136","#cd9a5b","#007947","#543044","#c63c26","#cd9a5b","#40835e","#63434f","#f3715c","#b36d41","#2b6447","#7d5886","#a7573b","#df9464","#005831","#401c44","#aa2116","#b76f40","#006c54","#472d56","#b64533","#ad8b3d","#375830","#45224a","#b54334","#dea32c","#274d3d","#411445","#853f04","#d1923f","#375830","#4b2f3d","#840228","#c88400","#27342b","#402e4c","#7a1723","#c37e00","#65c294","#c77eb5","#a03939","#c37e00","#73b9a2","#ea66a6","#8a2e3b","#e0861a","#72baa7","#f173ac","#8e453f","#ffce7b","#005344","#8f4b4a","#fcaf17","#122e29","#892f1b","#ba8448","#293047","#f6f5ec","#6b2c25","#896a45","#00ae9d","#d9d6c3","#733a31","#76624c","#508a88","#d1c7b7","#54211d","#6d5826","#70a19f","#f2eada","#78331e","#ffc20e","#50b7c1","#d3d7d4","#53261f","#fdb933","#00a6ac","#999d9c","#f15a22","#d3c6a6","#78cdd1","#a1a3a6","#b4533c","#c7a252","#008792","#9d9087","#84331f","#dec674","#94d6da","#8a8c8e","#f47a55","#b69968","#afdfe4","#74787c","#f15a22","#c1a173","#5e7c85","#7c8577","#f3704b","#dbce8f","#76becc","#72777b","#da765b","#ffd400","#90d7ec","#77787b","#c85d44","#ffd400","#009ad6","#4f5555","#ae5039","#ffe600","#145b7d","#6c4c49","#6a3427","#f0dc70","#11264f","#563624","#8f4b38","#fcf16e","#7bbfea","#3e4145","#8e3e1f","#decb00","#33a3dc","#3c3645","#f36c21","#cbc547","#228fbd","#464547","#b4532a","#6e6b41","#2468a2","#130c0e","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
cccol01 <- paste(cccol,"01",sep="")
cccol05 <- paste(cccol,"05",sep="")
cccol30 <- paste(cccol,"30",sep="")
cccol50 <- paste(cccol,"50",sep="")
cccol80 <- paste(cccol,"80",sep="")
options(scipen = 5)
red_col <- colorRampPalette(c(cccol[1],cccol[8],cccol[3],cccol[6]), bias=1)(4)
blue_col <- c(colorRampPalette(c(cccol[4],cccol[7]), bias=1)(5),colorRampPalette(c(cccol[2],cccol[11],cccol[9],cccol[14]), bias=1)(5),colorRampPalette(c(cccol[16],cccol[15]), bias=1)(5))
grey_col <- colorRampPalette(c(cccol[8],"black"), bias=1)(10)

# blue_col <- paste(blue_col,"50",sep="")


library(Seurat)
library(hdf5r)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
ESC <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/ESC/GeneSymbol_counts.txt", header=T, row.names=1)
ESC <- ESC[,c("PAP25","PA_4P7","IVFP25")]

P1 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271345_expmat.txt",row.names=1,header=T)
P2 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271346_expmat.txt",row.names=1,header=T)
P3 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271347_expmat.txt",row.names=1,header=T)
P4 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271348_expmat.txt",row.names=1,header=T)
P5 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271349_expmat.txt",row.names=1,header=T)
P6 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271350_expmat.txt",row.names=1,header=T)
P7 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271351_expmat.txt",row.names=1,header=T)
P8 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271352_expmat.txt",row.names=1,header=T)
P9 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271353_expmat.txt",row.names=1,header=T)
P10 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271354_expmat.txt",row.names=1,header=T)
P11 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271355_expmat.txt",row.names=1,header=T)
P12 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271356_expmat.txt",row.names=1,header=T)
P13 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271357_expmat.txt",row.names=1,header=T)
P14 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271358_expmat.txt",row.names=1,header=T)
P15 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271359_expmat.txt",row.names=1,header=T)
P16 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271360_expmat.txt",row.names=1,header=T)
P17 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271361_expmat.txt",row.names=1,header=T)
P18 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271362_expmat.txt",row.names=1,header=T)
P19 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271363_expmat.txt",row.names=1,header=T)
P20 <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/CRR271364_expmat.txt",row.names=1,header=T)

meta_data <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/meta.txt",header=T,row.names=2,sep="\t")

obj_ESC <- CreateSeuratObject(counts=ESC)
obj_P1 <- CreateSeuratObject(counts=P1)
obj_P2 <- CreateSeuratObject(counts=P2)
obj_P3 <- CreateSeuratObject(counts=P3)
obj_P4 <- CreateSeuratObject(counts=P4)
obj_P5 <- CreateSeuratObject(counts=P5)
obj_P6 <- CreateSeuratObject(counts=P6)
obj_P7 <- CreateSeuratObject(counts=P7)
obj_P8 <- CreateSeuratObject(counts=P8)
obj_P9 <- CreateSeuratObject(counts=P9)
obj_P10 <- CreateSeuratObject(counts=P10)
obj_P11 <- CreateSeuratObject(counts=P11)
obj_P12 <- CreateSeuratObject(counts=P12)
obj_P13 <- CreateSeuratObject(counts=P13)
obj_P14 <- CreateSeuratObject(counts=P14)
obj_P15 <- CreateSeuratObject(counts=P15)
obj_P16 <- CreateSeuratObject(counts=P16)
obj_P17 <- CreateSeuratObject(counts=P17)
obj_P18 <- CreateSeuratObject(counts=P18)
obj_P19 <- CreateSeuratObject(counts=P19)
obj_P20 <- CreateSeuratObject(counts=P20)

all_genes <- intersect(rownames(obj_ESC), rownames(obj_P1))

rm(ESC)
rm(P1)
rm(P2)
rm(P3)
rm(P4)
rm(P5)
rm(P6)
rm(P7)
rm(P8)
rm(P9)
rm(P10)
rm(P11)
rm(P12)
rm(P13)
rm(P14)
rm(P15)
rm(P16)
rm(P17)
rm(P18)
rm(P19)
rm(P20)

# normalization
obj_ESC <- NormalizeData(object = obj_ESC, normalization.method = "LogNormalize")
obj_P1 <- NormalizeData(object = obj_P1, normalization.method = "LogNormalize")
obj_P2 <- NormalizeData(object = obj_P2, normalization.method = "LogNormalize")
obj_P3 <- NormalizeData(object = obj_P3, normalization.method = "LogNormalize")
obj_P4 <- NormalizeData(object = obj_P4, normalization.method = "LogNormalize")
obj_P5 <- NormalizeData(object = obj_P5, normalization.method = "LogNormalize")
obj_P6 <- NormalizeData(object = obj_P6, normalization.method = "LogNormalize")
obj_P7 <- NormalizeData(object = obj_P7, normalization.method = "LogNormalize")
obj_P8 <- NormalizeData(object = obj_P8, normalization.method = "LogNormalize")
obj_P9 <- NormalizeData(object = obj_P9, normalization.method = "LogNormalize")
obj_P10 <- NormalizeData(object = obj_P10, normalization.method = "LogNormalize")
obj_P11 <- NormalizeData(object = obj_P11, normalization.method = "LogNormalize")
obj_P12 <- NormalizeData(object = obj_P12, normalization.method = "LogNormalize")
obj_P13 <- NormalizeData(object = obj_P13, normalization.method = "LogNormalize")
obj_P14 <- NormalizeData(object = obj_P14, normalization.method = "LogNormalize")
obj_P15 <- NormalizeData(object = obj_P15, normalization.method = "LogNormalize")
obj_P16 <- NormalizeData(object = obj_P16, normalization.method = "LogNormalize")
obj_P17 <- NormalizeData(object = obj_P17, normalization.method = "LogNormalize")
obj_P18 <- NormalizeData(object = obj_P18, normalization.method = "LogNormalize")
obj_P19 <- NormalizeData(object = obj_P19, normalization.method = "LogNormalize")
obj_P20 <- NormalizeData(object = obj_P20, normalization.method = "LogNormalize")

obj_ESC <- FindVariableFeatures(object = obj_ESC, selection.method = "vst", nfeatures = 2000)
obj_P1 <- FindVariableFeatures(object = obj_P1, selection.method = "vst", nfeatures = 2000)
obj_P2 <- FindVariableFeatures(object = obj_P2, selection.method = "vst", nfeatures = 2000)
obj_P3 <- FindVariableFeatures(object = obj_P3, selection.method = "vst", nfeatures = 2000)
obj_P4 <- FindVariableFeatures(object = obj_P4, selection.method = "vst", nfeatures = 2000)
obj_P5 <- FindVariableFeatures(object = obj_P5, selection.method = "vst", nfeatures = 2000)
obj_P6 <- FindVariableFeatures(object = obj_P6, selection.method = "vst", nfeatures = 2000)
obj_P7 <- FindVariableFeatures(object = obj_P7, selection.method = "vst", nfeatures = 2000)
obj_P8 <- FindVariableFeatures(object = obj_P8, selection.method = "vst", nfeatures = 2000)
obj_P9 <- FindVariableFeatures(object = obj_P9, selection.method = "vst", nfeatures = 2000)
obj_P10 <- FindVariableFeatures(object = obj_P10, selection.method = "vst", nfeatures = 2000)
obj_P11 <- FindVariableFeatures(object = obj_P11, selection.method = "vst", nfeatures = 2000)
obj_P12 <- FindVariableFeatures(object = obj_P12, selection.method = "vst", nfeatures = 2000)
obj_P13 <- FindVariableFeatures(object = obj_P13, selection.method = "vst", nfeatures = 2000)
obj_P14 <- FindVariableFeatures(object = obj_P14, selection.method = "vst", nfeatures = 2000)
obj_P15 <- FindVariableFeatures(object = obj_P15, selection.method = "vst", nfeatures = 2000)
obj_P16 <- FindVariableFeatures(object = obj_P16, selection.method = "vst", nfeatures = 2000)
obj_P17 <- FindVariableFeatures(object = obj_P17, selection.method = "vst", nfeatures = 2000)
obj_P18 <- FindVariableFeatures(object = obj_P18, selection.method = "vst", nfeatures = 2000)
obj_P19 <- FindVariableFeatures(object = obj_P19, selection.method = "vst", nfeatures = 2000)
obj_P20 <- FindVariableFeatures(object = obj_P20, selection.method = "vst", nfeatures = 2000)

# marker genes
EPI <- cbind("SOX2", "POU5F1", "OTX2")
TE <- cbind("GATA3", "GATA2")
PE <- cbind("GATA6", "GATA4")
EPI <- intersect(all_genes,EPI)
TE <- intersect(all_genes,TE)
PE <- intersect(all_genes,PE)

###################################################################################################
################################            development            ################################
###################################################################################################

# RNA
p_integr_dev_list <- list(obj_P1,obj_P2,obj_P3,obj_P4,obj_P5,obj_P6,obj_P7,obj_P8,obj_P9,obj_P10,obj_P11,obj_P12,obj_P13,obj_P14,obj_P15,obj_P16,obj_P17,obj_P18,obj_P19,obj_P20)
p_integr_dev_features <- SelectIntegrationFeatures(object.list = p_integr_dev_list, nfeatures = 2000)
p_integr_dev_anchors <- FindIntegrationAnchors(object.list = p_integr_dev_list, reference = c(1, 3), anchor.features = p_integr_dev_features)
p_integr_dev <- IntegrateData(anchorset = p_integr_dev_anchors, normalization.method = "LogNormalize", features.to.integrate = all_genes)

samples <- c(rep("P1",nrow(obj_P1@meta.data)),rep("P2",nrow(obj_P2@meta.data)),rep("P3",nrow(obj_P3@meta.data)),rep("P4",nrow(obj_P4@meta.data)),rep("P5",nrow(obj_P5@meta.data)),rep("P6",nrow(obj_P6@meta.data)),rep("P7",nrow(obj_P7@meta.data)),rep("P8",nrow(obj_P8@meta.data)),rep("P9",nrow(obj_P9@meta.data)),rep("P10",nrow(obj_P10@meta.data)),rep("P11",nrow(obj_P11@meta.data)),rep("P12",nrow(obj_P12@meta.data)),rep("P13",nrow(obj_P13@meta.data)),rep("P14",nrow(obj_P14@meta.data)),rep("P15",nrow(obj_P15@meta.data)),rep("P16",nrow(obj_P16@meta.data)),rep("P17",nrow(obj_P17@meta.data)),rep("P18",nrow(obj_P18@meta.data)),rep("P19",nrow(obj_P19@meta.data)),rep("P20",nrow(obj_P20@meta.data)))
p_integr_dev@meta.data$samples <- samples

dev_cells <- c(intersect(colnames(p_integr_dev),rownames(meta_data)))
p_integr_dev <- subset(p_integr_dev, cells = dev_cells)

celltypes <- rep("NA", ncol(p_integr_dev)); names(celltypes) <- colnames(p_integr_dev)
celltypes[intersect(colnames(p_integr_dev),rownames(meta_data))] <- meta_data[intersect(colnames(p_integr_dev),rownames(meta_data)),2]
p_integr_dev@meta.data$celltypes <- celltypes

p_integr_dev <- FindVariableFeatures(object = p_integr_dev, assay = "RNA")
p_integr_dev <- ScaleData(object = p_integr_dev, assay = "RNA")
p_integr_dev <- RunPCA(object = p_integr_dev, assay = "RNA")
p_integr_dev <- RunUMAP(object = p_integr_dev, assay = "RNA", reduction = "pca", dims = 1:30, n.neighbors = 30)
p_integr_dev <- RunTSNE(object = p_integr_dev, assay = "RNA", reduction = "pca", dims = 1:30)
write.table(p_integr_dev@reductions$umap@cell.embeddings,file="p_integr_dev_DimensinalResult_UMAP.Raw.txt",quote=F,sep="\t",col.names=T,row.names=T)
write.table(p_integr_dev@reductions$tsne@cell.embeddings,file="p_integr_dev_DimensinalResult_tSNE.Raw.txt",quote=F,sep="\t",col.names=T,row.names=T)

pdf("p_integr_dev_DimensinalPlot_ColoredBySamples.Raw.pdf",width=6.5,height=5)
DimPlot(p_integr_dev, reduction = "pca", group.by = "samples", cols = allcol)
DimPlot(p_integr_dev, reduction = "umap", group.by = "samples", cols = allcol)
DimPlot(p_integr_dev, reduction = "tsne", group.by = "samples", cols = allcol)
DimPlot(p_integr_dev, reduction = "pca", group.by = "celltypes", cols = c(blue_col[c(3,4,5,11,12,13,14,15,7,8,9,10)],blue_col[c(6,1)],blue_col[2]))
DimPlot(p_integr_dev, reduction = "umap", group.by = "celltypes", cols = c(blue_col[c(3,4,5,11,12,13,14,15,7,8,9,10)],blue_col[c(6,1)],blue_col[2]))
DimPlot(p_integr_dev, reduction = "tsne", group.by = "celltypes", cols = c(blue_col[c(3,4,5,11,12,13,14,15,7,8,9,10)],blue_col[c(6,1)],blue_col[2]))
dev.off()

DefaultAssay(p_integr_dev) <- "RNA"
p_integr_dev_umap <- FindNeighbors(object = p_integr_dev, reduction = "umap", dims = 1:2, k.param = 30)
p_integr_dev_umap <- FindClusters(object = p_integr_dev_umap, resolution = 0.15)
clusters <- paste("C",as.numeric(p_integr_dev_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_integr_dev_umap@assays$RNA[,])
write.table(clusters,file="p_integr_dev_umap_ClusteringResult.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_integr_dev_umap_DimensinalPlot_ColoredByClusters.pdf",width=5.6,height=5)
DimPlot(object = p_integr_dev_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_integr_dev_umap) <- "RNA"
pdf("p_integr_dev_umap_DimensinalPlot_ColoredByEPI.pdf",width=5,height=5)
FeaturePlot(p_integr_dev_umap, features = EPI[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_dev_umap, features = EPI[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_integr_dev_umap, features = EPI, cols = c("lightgrey",cccol[2]))
dev.off()
pdf("p_integr_dev_umap_DimensinalPlot_ColoredByTE.pdf",width=5,height=5)
FeaturePlot(p_integr_dev_umap, features = TE[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_dev_umap, features = TE[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_integr_dev_umap, features = TE, cols = c("lightgrey",cccol[2]))
dev.off()
pdf("p_integr_dev_umap_DimensinalPlot_ColoredByPE.pdf",width=5,height=5)
FeaturePlot(p_integr_dev_umap, features = PE[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_dev_umap, features = PE[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_integr_dev_umap, features = PE, cols = c("lightgrey",cccol[2]))
dev.off()

EPI_score <- apply(p_integr_dev_umap@assays$RNA[EPI,], 2, mean, na.rm=T)
TE_score <- apply(p_integr_dev_umap@assays$RNA[TE,], 2, mean, na.rm=T)
PE_score <- apply(p_integr_dev_umap@assays$RNA[PE,], 2, mean, na.rm=T)
p_integr_dev_umap@meta.data$EPI_score <- EPI_score
p_integr_dev_umap@meta.data$TE_score <- TE_score
p_integr_dev_umap@meta.data$PE_score <- PE_score
pdf("p_integr_dev_umap_DimensinalPlot_ColoredByCellTypeFeatures.pdf",width=5,height=5)
plot(FeaturePlot(p_integr_dev_umap, features = "EPI_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
plot(FeaturePlot(p_integr_dev_umap, features = "TE_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
plot(FeaturePlot(p_integr_dev_umap, features = "PE_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
DotPlot(object = p_integr_dev_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

p_integr_dev_tsne <- FindNeighbors(object = p_integr_dev, reduction = "tsne", dims = 1:2, k.param = 30)
p_integr_dev_tsne <- FindClusters(object = p_integr_dev_tsne, resolution = 0.15)
clusters <- paste("C",as.numeric(p_integr_dev_tsne@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_integr_dev_tsne@assays$RNA[,])
write.table(clusters,file="p_integr_dev_tsne_ClusteringResult.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_integr_dev_tsne_DimensinalPlot_ColoredByClusters.pdf",width=5.6,height=5)
DimPlot(object = p_integr_dev_tsne, reduction = "tsne", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_integr_dev_tsne) <- "RNA"
pdf("p_integr_dev_tsne_DimensinalPlot_ColoredByEPI.pdf",width=5,height=5)
FeaturePlot(p_integr_dev_tsne, features = EPI[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_dev_tsne, features = EPI[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_integr_dev_tsne, features = EPI, cols = c("lightgrey",cccol[2]))
dev.off()
pdf("p_integr_dev_tsne_DimensinalPlot_ColoredByTE.pdf",width=5,height=5)
FeaturePlot(p_integr_dev_tsne, features = TE[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_dev_tsne, features = TE[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_integr_dev_tsne, features = TE, cols = c("lightgrey",cccol[2]))
dev.off()
pdf("p_integr_dev_tsne_DimensinalPlot_ColoredByPE.pdf",width=5,height=5)
FeaturePlot(p_integr_dev_tsne, features = PE[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_dev_tsne, features = PE[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_integr_dev_tsne, features = PE, cols = c("lightgrey",cccol[2]))
dev.off()

EPI_score <- apply(p_integr_dev_tsne@assays$RNA[EPI,], 2, mean, na.rm=T)
TE_score <- apply(p_integr_dev_tsne@assays$RNA[TE,], 2, mean, na.rm=T)
PE_score <- apply(p_integr_dev_tsne@assays$RNA[PE,], 2, mean, na.rm=T)
p_integr_dev_tsne@meta.data$EPI_score <- EPI_score
p_integr_dev_tsne@meta.data$TE_score <- TE_score
p_integr_dev_tsne@meta.data$PE_score <- PE_score
pdf("p_integr_dev_tsne_DimensinalPlot_ColoredByCellTypeFeatures.pdf",width=5,height=5)
plot(FeaturePlot(p_integr_dev_tsne, features = "EPI_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
plot(FeaturePlot(p_integr_dev_tsne, features = "TE_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
plot(FeaturePlot(p_integr_dev_tsne, features = "PE_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
DotPlot(object = p_integr_dev_tsne, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

# p_integr_ESC_umap_anno <- RenameIdents(p_integr_ESC_umap,
#                           `0` = "TE",
#                           `1` = "TE",
#                           `2` = "TE",
#                           `3` = "TE",
#                           `4` = "TE",
#                           `5` = "EPI",
#                           `6` = "TE",
#                           `7` = "PE",
#                           `8` = "Undefined1",
#                           `9` = "Undefined2",
#                           `10` = "PE",
#                           `11` = "Undefined3",
#                           `12` = "TE")
# pdf("p_integr_ESC_umap_DimensinalPlot_ColoredByCellTypes.pdf",width=6.2,height=5)
# DimPlot(object = p_integr_ESC_umap_anno, reduction = "umap", cols = cccol, label = TRUE, label.color = "white", label.box = TRUE, repel = TRUE)
# DimPlot(object = p_integr_ESC_umap_anno, reduction = "umap", cols = cccol, repel = TRUE)
# dev.off()
# saveRDS(p_integr_ESC_umap_anno, "p_integr_ESC_umap_IntergratedMatrixAnnotated.rds")

###################################################################################################
################################          including ESC            ################################
###################################################################################################
p_integr_ESC <- CreateSeuratObject(counts=cbind(obj_ESC@assays$RNA[all_genes,],p_integr_dev@assays$RNA[all_genes,]))


# # identify anchors
# p_integr_ESC_list <- list(obj_ESC,obj_P1,obj_P2,obj_P3,obj_P4,obj_P5,obj_P6,obj_P7,obj_P8,obj_P9,obj_P10,obj_P11,obj_P12,obj_P13,obj_P14,obj_P15,obj_P16,obj_P17,obj_P18,obj_P19,obj_P20)
# p_integr_ESC_features <- SelectIntegrationFeatures(object.list = p_integr_ESC_list, nfeatures = 2000)
#     # p_integr_ESC_anchors <- FindIntegrationAnchors(object.list = p_integr_ESC_list, anchor.features = p_integr_ESC_features)
#     # If you have many replicates, setting a reference object will effectively save memory usage.
# p_integr_ESC_anchors <- FindIntegrationAnchors(object.list = p_integr_ESC_list, reference = c(1, 2), anchor.features = p_integr_ESC_features, dims=1:2, k.anchor = 2, k.filter = 2, k.score = 2)

# # integrated (or batch-corrected)
# p_integr_ESC <- IntegrateData(anchorset = p_integr_ESC_anchors, normalization.method = "LogNormalize", features.to.integrate = all_genes)

# # annotated infomation
# samples <- c(rep("ESC",nrow(obj_ESC@meta.data)),rep("P1",nrow(obj_P1@meta.data)),rep("P2",nrow(obj_P2@meta.data)),rep("P3",nrow(obj_P3@meta.data)),rep("P4",nrow(obj_P4@meta.data)),rep("P5",nrow(obj_P5@meta.data)),rep("P6",nrow(obj_P6@meta.data)),rep("P7",nrow(obj_P7@meta.data)),rep("P8",nrow(obj_P8@meta.data)),rep("P9",nrow(obj_P9@meta.data)),rep("P10",nrow(obj_P10@meta.data)),rep("P11",nrow(obj_P11@meta.data)),rep("P12",nrow(obj_P12@meta.data)),rep("P13",nrow(obj_P13@meta.data)),rep("P14",nrow(obj_P14@meta.data)),rep("P15",nrow(obj_P15@meta.data)),rep("P16",nrow(obj_P16@meta.data)),rep("P17",nrow(obj_P17@meta.data)),rep("P18",nrow(obj_P18@meta.data)),rep("P19",nrow(obj_P19@meta.data)),rep("P20",nrow(obj_P20@meta.data)))
# p_integr_ESC@meta.data$samples <- samples

# dev_cells <- c(intersect(colnames(p_integr_ESC),rownames(meta_data)))
# p_integr_ESC <- subset(p_integr_ESC, cells = dev_cells)

celltypes <- rep("0.ESC", ncol(p_integr_ESC)); names(celltypes) <- colnames(p_integr_ESC)
celltypes[intersect(colnames(p_integr_ESC),rownames(meta_data))] <- meta_data[intersect(colnames(p_integr_ESC),rownames(meta_data)),2]
p_integr_ESC@meta.data$celltypes <- celltypes

ptsharp <- rep("1.ESC", ncol(p_integr_ESC)); names(ptsharp) <- colnames(p_integr_ESC)
ptsharp[intersect(colnames(p_integr_ESC),rownames(meta_data))]<- "0.Embryo"
p_integr_ESC@meta.data$ptsharp <- ptsharp


# RNA
p_integr_ESC <- FindVariableFeatures(object = p_integr_ESC, assay = "RNA")
p_integr_ESC <- ScaleData(object = p_integr_ESC, assay = "RNA")
p_integr_ESC <- RunPCA(object = p_integr_ESC, assay = "RNA")
p_integr_ESC <- RunUMAP(object = p_integr_ESC, assay = "RNA", reduction = "pca", dims = 1:30, n.neighbors = 30)
p_integr_ESC <- RunTSNE(object = p_integr_ESC, assay = "RNA", reduction = "pca", dims = 1:30)
write.table(p_integr_ESC@reductions$umap@cell.embeddings,file="p_integr_ESC_DimensinalResult_UMAP.Raw.txt",quote=F,sep="\t",col.names=T,row.names=T)
write.table(p_integr_ESC@reductions$tsne@cell.embeddings,file="p_integr_ESC_DimensinalResult_tSNE.Raw.txt",quote=F,sep="\t",col.names=T,row.names=T)

pdf("p_integr_ESC_DimensinalPlot_ColoredBySamples.Raw.pdf",width=6.5,height=5)
DimPlot(p_integr_ESC, reduction = "pca", group.by = "celltypes",shape.by = "ptsharp", cols = c(red_col[1],blue_col[c(3,4,5,11,12,13,14,15,7,8,9,10)],blue_col[c(6,1)],blue_col[2]),shuffle=T)
DimPlot(p_integr_ESC, reduction = "umap", group.by = "celltypes",shape.by = "ptsharp", cols = c(red_col[1],blue_col[c(3,4,5,11,12,13,14,15,7,8,9,10)],blue_col[c(6,1)],blue_col[2]),shuffle=T)
DimPlot(p_integr_ESC, reduction = "tsne", group.by = "celltypes",shape.by = "ptsharp", cols = c(red_col[1],blue_col[c(3,4,5,11,12,13,14,15,7,8,9,10)],blue_col[c(6,1)],blue_col[2]),shuffle=T)
dev.off()

p_integr_ESC_umap <- FindNeighbors(object = p_integr_ESC, reduction = "umap", dims = 1:2, k.param = 30)
p_integr_ESC_umap <- FindClusters(object = p_integr_ESC_umap, resolution = 0.15)
clusters <- paste("C",as.numeric(p_integr_ESC_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_integr_ESC_umap@assays$RNA[,])
write.table(clusters,file="p_integr_ESC_umap_ClusteringResult.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_integr_ESC_umap_DimensinalPlot_ColoredByClusters.pdf",width=5.6,height=5)
DimPlot(object = p_integr_ESC_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_integr_ESC_umap) <- "RNA"
pdf("p_integr_ESC_pca_DimensinalPlot_ColoredByEPI.pdf",width=5,height=5)
FeaturePlot(p_integr_ESC, features = EPI[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "pca", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_ESC, features = EPI[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "pca", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_ESC, features = EPI[3], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "pca", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_integr_ESC, features = EPI, cols = c("lightgrey",cccol[2]))
dev.off()

DefaultAssay(p_integr_ESC_umap) <- "RNA"
pdf("p_integr_ESC_umap_DimensinalPlot_ColoredByEPI.pdf",width=5,height=5)
FeaturePlot(p_integr_ESC_umap, features = EPI[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_ESC_umap, features = EPI[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_ESC_umap, features = EPI[3], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_integr_ESC_umap, features = EPI, cols = c("lightgrey",cccol[2]))
dev.off()
pdf("p_integr_ESC_umap_DimensinalPlot_ColoredByTE.pdf",width=5,height=5)
FeaturePlot(p_integr_ESC_umap, features = TE[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_ESC_umap, features = TE[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_integr_ESC_umap, features = TE, cols = c("lightgrey",cccol[2]))
dev.off()
pdf("p_integr_ESC_umap_DimensinalPlot_ColoredByPE.pdf",width=5,height=5)
FeaturePlot(p_integr_ESC_umap, features = PE[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_integr_ESC_umap, features = PE[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_integr_ESC_umap, features = PE, cols = c("lightgrey",cccol[2]))
dev.off()

EPI_score <- apply(p_integr_ESC_umap@assays$RNA[EPI,], 2, mean, na.rm=T)
TE_score <- apply(p_integr_ESC_umap@assays$RNA[TE,], 2, mean, na.rm=T)
PE_score <- apply(p_integr_ESC_umap@assays$RNA[PE,], 2, mean, na.rm=T)
p_integr_ESC_umap@meta.data$EPI_score <- EPI_score
p_integr_ESC_umap@meta.data$TE_score <- TE_score
p_integr_ESC_umap@meta.data$PE_score <- PE_score
pdf("p_integr_ESC_umap_DimensinalPlot_ColoredByCellTypeFeatures.pdf",width=5,height=5)
plot(FeaturePlot(p_integr_ESC_umap, features = "EPI_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
plot(FeaturePlot(p_integr_ESC_umap, features = "TE_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
plot(FeaturePlot(p_integr_ESC_umap, features = "PE_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
DotPlot(object = p_integr_ESC_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

# p_integr_ESC_umap_anno <- RenameIdents(p_integr_ESC_umap,
#                           `0` = "TE",
#                           `1` = "TE",
#                           `2` = "TE",
#                           `3` = "TE",
#                           `4` = "TE",
#                           `5` = "EPI",
#                           `6` = "TE",
#                           `7` = "PE",
#                           `8` = "Undefined1",
#                           `9` = "Undefined2",
#                           `10` = "PE",
#                           `11` = "Undefined3",
#                           `12` = "TE")
# pdf("p_integr_ESC_umap_DimensinalPlot_ColoredByCellTypes.pdf",width=6.2,height=5)
# DimPlot(object = p_integr_ESC_umap_anno, reduction = "umap", cols = cccol, label = TRUE, label.color = "white", label.box = TRUE, repel = TRUE)
# DimPlot(object = p_integr_ESC_umap_anno, reduction = "umap", cols = cccol, repel = TRUE)
# dev.off()
# saveRDS(p_integr_ESC_umap_anno, "p_integr_ESC_umap_IntergratedMatrixAnnotated.rds")

# pdf("p_integr_ESC_umap_DotPlot_ColoredByCellTypes.pdf",width=8,height=4)
# DotPlot(object = p_integr_ESC_umap_anno, features = c(EPI, TE, PE), cols = c("white",cccol[1]))
# dev.off()

# anno <- p_integr_ESC_umap_anno@active.ident
# write.table(anno,file="p_integr_ESC_umap_anno_ClusteringResult.txt",quote=F,sep="\t",col.names=F,row.names=T)

# # markers
# DefaultAssay(p_integr_ESC_umap_anno) <- "RNA"
# p_integr_ESC_umap_anno_MarkerTest <- FindAllMarkers(object = p_integr_ESC_umap_anno, logfc.threshold = 0, min.pct = 0, return.thresh = 0.01)
# # min.pct    only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
# # min.diff.pct    only test genes that show a minimum difference in the fraction of detection between the two groups
# # logfc.threshold
# # return.thresh    Only return markers that have a p-value < return.thresh
# write.table(p_integr_ESC_umap_anno_MarkerTest,file="p_integr_ESC_umap_anno_MarkerTest.txt",quote=F,sep="\t",col.names=T,row.names=T)

# pct_cutoff <- 0.3
# logfc_cutoff <- 0.8
# pdf("p_integr_ESC_umap_anno_ClusterMarkers.VolcanoPlot.pdf",width=4*3,height=4*2)
# par(mfrow=c(2,3),mar=c(6,6,4,4))
# for (i in names(table(p_integr_ESC_umap_anno_MarkerTest[,"cluster"]))){
#     tmp_data <- p_integr_ESC_umap_anno_MarkerTest[which(p_integr_ESC_umap_anno_MarkerTest[,"cluster"]==i),]
#     x_pctdiff <- tmp_data[,"pct.1"] - tmp_data[,"pct.2"]; names(x_pctdiff) <- tmp_data[,"gene"]
#     y_logfc <- tmp_data[,"avg_log2FC"]; names(y_logfc) <- tmp_data[,"gene"]
#     up_genes <- tmp_data[,"gene"][which((x_pctdiff > pct_cutoff) & (y_logfc > logfc_cutoff))]
#     dn_genes <- tmp_data[,"gene"][which((x_pctdiff < -pct_cutoff) & (y_logfc < -logfc_cutoff))]
#     ns_genes <- setdiff(tmp_data[,"gene"], c(up_genes, dn_genes))

#     plot(x_pctdiff, y_logfc, main=paste("",i,sep=""), xlab="Delta Percentage Diff (%)", ylab="log (Fold Change)", type="p", pch=".", col=cccol[8],xlim=c(-1,1),ylim=c(-4,4.2)); box(lwd=2)
#     abline(v=0, h=0, lty=2)
#     if (length(up_genes)>0){
#         points(x_pctdiff[up_genes], y_logfc[up_genes], pch=16, col=cccol[1])
#         text(x_pctdiff[up_genes]+0.15, y_logfc[up_genes], labels=up_genes, cex=0.5)
#     }
#     if (length(dn_genes)>0){
#         points(x_pctdiff[dn_genes], y_logfc[dn_genes], pch=16, col=cccol[2])
#         text(x_pctdiff[dn_genes]+0.15, y_logfc[dn_genes], labels=dn_genes, cex=0.5)
#     }
# }

# par(mfrow=c(2,3),mar=c(6,6,4,4))
# for (i in names(table(p_integr_ESC_umap_anno_MarkerTest[,"cluster"]))){
#     tmp_data <- p_integr_ESC_umap_anno_MarkerTest[which(p_integr_ESC_umap_anno_MarkerTest[,"cluster"]==i),]
#     x_pctdiff <- tmp_data[,"pct.1"] - tmp_data[,"pct.2"]; names(x_pctdiff) <- tmp_data[,"gene"]
#     y_logfc <- tmp_data[,"avg_log2FC"]; names(y_logfc) <- tmp_data[,"gene"]
#     up_genes <- tmp_data[,"gene"][which((x_pctdiff > pct_cutoff) & (y_logfc > logfc_cutoff))]
#     dn_genes <- tmp_data[,"gene"][which((x_pctdiff < -pct_cutoff) & (y_logfc < -logfc_cutoff))]
#     ns_genes <- setdiff(tmp_data[,"gene"], c(up_genes, dn_genes))

#     plot(x_pctdiff, y_logfc, main=paste("",i,sep=""), xlab="Delta Percentage Diff (%)", ylab="log (Fold Change)", type="p", pch=".", col=cccol[8]); box(lwd=2)
#     abline(v=0, h=0, lty=2)
#     if (length(up_genes)>0){
#         points(x_pctdiff[up_genes], y_logfc[up_genes], pch=16, col=cccol[1])
#         # text(x_pctdiff[up_genes]+0.15, y_logfc[up_genes], labels=up_genes, cex=0.5)
#     }
#     if (length(dn_genes)>0){
#         points(x_pctdiff[dn_genes], y_logfc[dn_genes], pch=16, col=cccol[2])
#         # text(x_pctdiff[dn_genes]+0.15, y_logfc[dn_genes], labels=dn_genes, cex=0.5)
#     }
# }
# dev.off()

# p_integr_ESC_umap_anno_Markers <- FindAllMarkers(object = p_integr_ESC_umap_anno, min.pct = 0, logfc.threshold = logfc_cutoff, min.diff.pct = pct_cutoff, return.thresh = 0.01)
# write.table(p_integr_ESC_umap_anno_Markers,file="p_integr_ESC_umap_anno_Markers.txt",quote=F,sep="\t",col.names=T,row.names=T)
# for (i in names(table(p_integr_ESC_umap_anno_Markers[,"cluster"]))){
#     tmp_data <- p_integr_ESC_umap_anno_Markers[which(p_integr_ESC_umap_anno_Markers[,"cluster"]==i),]
#     write.table(tmp_data[,"gene"],file=paste("p_integr_ESC_umap_anno_Markers_",i,".txt",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
# }
# cluster_markers <- c()
# for (i in names(table(p_integr_ESC_umap_anno_MarkerTest[,"cluster"]))){
#     cluster_markers <- c(cluster_markers, p_integr_ESC_umap_anno_Markers[which(p_integr_ESC_umap_anno_Markers[,"cluster"]==i & p_integr_ESC_umap_anno_Markers[,"avg_log2FC"] >0 ),"gene"][1])
# }
# pdf("p_integr_ESC_umap_anno_DimensinalPlot_ColoredByClusterMarkers.pdf",width=5*3,height=5*2)
# plot(FeaturePlot(p_integr_ESC_umap_anno, features = cluster_markers, min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "umap", pt.size = 0.1, ncol = 3, keep.scale = "all") & NoAxes() & NoLegend())
# DotPlot(object = p_integr_ESC_umap_anno, features = cluster_markers, cols = c(cccol[2],cccol[1]), dot.scale = 9)
# dev.off()


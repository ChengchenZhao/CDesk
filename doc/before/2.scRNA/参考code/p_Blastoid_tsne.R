###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#f7acbc","#deab8a","#817936","#444693","#ef5b9c","#fedcbd","#7f7522","#2b4490","#f47920","#80752c","#2a5caa","#f05b72","#905a3d","#87843b","#224b8f","#f15b6c","#8f4b2e","#726930","#003a6c","#f8aba6","#87481f","#454926","#52b6a","#f69c9f","#5f3c23","#2e3a1f","#426ab3","#f58f98","#6b473c","#4d4f36","#46485f","#ca8687","#faa755","#b7ba6b","#4e72b8","#f391a9","#fab27b","#b2d235","#181d4b","#bd6758","#f58220","#5c7a29","#1a2933","#d71345","#843900","#bed742","#121a2a","#d64f44","#905d1d","#7fb80e","#0c212b","#d93a49","#8a5d19","#a3cf62","#6a6da9","#b3424a","#8c531b","#769149","#585eaa","#c76968","#826858","#6d8346","#494e8f","#bb505d","#64492b","#78a355","#afb4db","#987165","#ae6642","#abc88b","#9b95c9","#ac6767","#56452d","#74905d","#6950a1","#973c3f","#96582a","#cde6c7","#6f60aa","#b22c46","#705628","#1d953f","#867892","#a7324a","#4a3113","#77ac98","#918597","#aa363d","#412f1f","#007d65","#6f6d85","#ed1941","#845538","#84bf96","#594c6d","#f26522","#8e7437","#45b97c","#694d9f","#d2553d","#69541b","#225a1f","#6f599c","#b4534b","#d5c59f","#367459","#8552a1","#ef4136","#cd9a5b","#007947","#543044","#c63c26","#cd9a5b","#40835e","#63434f","#f3715c","#b36d41","#2b6447","#7d5886","#a7573b","#df9464","#005831","#401c44","#aa2116","#b76f40","#006c54","#472d56","#b64533","#ad8b3d","#375830","#45224a","#b54334","#dea32c","#274d3d","#411445","#853f04","#d1923f","#375830","#4b2f3d","#840228","#c88400","#27342b","#402e4c","#7a1723","#c37e00","#65c294","#c77eb5","#a03939","#c37e00","#73b9a2","#ea66a6","#8a2e3b","#e0861a","#72baa7","#f173ac","#8e453f","#ffce7b","#005344","#8f4b4a","#fcaf17","#122e29","#892f1b","#ba8448","#293047","#f6f5ec","#6b2c25","#896a45","#00ae9d","#d9d6c3","#733a31","#76624c","#508a88","#d1c7b7","#54211d","#6d5826","#70a19f","#f2eada","#78331e","#ffc20e","#50b7c1","#d3d7d4","#53261f","#fdb933","#00a6ac","#999d9c","#f15a22","#d3c6a6","#78cdd1","#a1a3a6","#b4533c","#c7a252","#008792","#9d9087","#84331f","#dec674","#94d6da","#8a8c8e","#f47a55","#b69968","#afdfe4","#74787c","#f15a22","#c1a173","#5e7c85","#7c8577","#f3704b","#dbce8f","#76becc","#72777b","#da765b","#ffd400","#90d7ec","#77787b","#c85d44","#ffd400","#009ad6","#4f5555","#ae5039","#ffe600","#145b7d","#6c4c49","#6a3427","#f0dc70","#11264f","#563624","#8f4b38","#fcf16e","#7bbfea","#3e4145","#8e3e1f","#decb00","#33a3dc","#3c3645","#f36c21","#cbc547","#228fbd","#464547","#b4532a","#6e6b41","#2468a2","#130c0e","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
cccol05 <- c("#CE001305","#16557A05","#C7A60905","#87C23205","#00879205","#A14C9405","#15A08C05","#8B7E7505","#1E7CAF05","#EA425F05","#46489A05","#E5003305","#0F231F05","#1187CD05")
options(scipen = 5)

library(Seurat)
library(hdf5r)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
p_Blastoid_1 <- Read10X_h5("/mnt/liudong/task/20230724_pig_sc/h5/InHouse/p_Blastoid_1.filter.h5")
p_Blastoid_2 <- Read10X_h5("/mnt/liudong/task/20230724_pig_sc/h5/InHouse/p_Blastoid_2.filter.h5")

obj_p_Blastoid_1 <- CreateSeuratObject(counts=p_Blastoid_1)
obj_p_Blastoid_2 <- CreateSeuratObject(counts=p_Blastoid_2)
obj_p_Blastoid_1 <- obj_p_Blastoid_1[, obj_p_Blastoid_1@meta.data$nCount_RNA >= 11027]
obj_p_Blastoid_2 <- obj_p_Blastoid_2[, obj_p_Blastoid_2@meta.data$nCount_RNA >= 9684]

# normalization
obj_p_Blastoid_1 <- NormalizeData(object = obj_p_Blastoid_1, normalization.method = "LogNormalize")
obj_p_Blastoid_2 <- NormalizeData(object = obj_p_Blastoid_2, normalization.method = "LogNormalize")

obj_p_Blastoid_1 <- FindVariableFeatures(object = obj_p_Blastoid_1, selection.method = "vst", nfeatures = 4000)
obj_p_Blastoid_2 <- FindVariableFeatures(object = obj_p_Blastoid_2, selection.method = "vst", nfeatures = 4000)

# marker genes
all_genes <- row.names(obj_p_Blastoid_1)
EPI <- cbind("SOX2", "POU5F1")
TE <- cbind("GATA3", "GATA2")
PE <- cbind("GATA6", "GATA4")
EPI <- intersect(all_genes,EPI)
TE <- intersect(all_genes,TE)
PE <- intersect(all_genes,PE)

# identify anchors
p_Blastoid_list <- list(obj_p_Blastoid_1,obj_p_Blastoid_2)
p_Blastoid_features <- SelectIntegrationFeatures(object.list = p_Blastoid_list, nfeatures = 4000)
p_Blastoid_anchors <- FindIntegrationAnchors(object.list = p_Blastoid_list, anchor.features = p_Blastoid_features)

# integrated (or batch-corrected)
p_Blastoid <- IntegrateData(anchorset = p_Blastoid_anchors, normalization.method = "LogNormalize", features.to.integrate = all_genes)

samples <- c(rep("p_Blastoid_1",nrow(obj_p_Blastoid_1@meta.data)),rep("p_Blastoid_2",nrow(obj_p_Blastoid_2@meta.data)))
p_Blastoid@meta.data$samples <- samples
write.table(file="p_Blastoid_IntergratedMatrix.txt",p_Blastoid@assays$integrated[,],sep="\t",col.names=T,row.names=T,quote=F)
write.table(file="p_Blastoid_RawCountsMatrix.txt",p_Blastoid@assays$RNA[,],sep="\t",col.names=T,row.names=T,quote=F)

p_Blastoid <- ScaleData(object = p_Blastoid)
p_Blastoid <- RunPCA(object = p_Blastoid)
p_Blastoid <- RunTSNE(object = p_Blastoid, reduction = "pca", dims = 1:30)
write.table(p_Blastoid@reductions$tsne@cell.embeddings,file="p_Blastoid_DimensinalResult_tSNE.txt",quote=F,sep="\t",col.names=T,row.names=T)

# pdf("p_Blastoid_DimensinalPlot_ColoredBySamples.pdf",width=6.4,height=5)
# # DimPlot(p_Blastoid, reduction = "umap", group.by = "samples", cols = cccol)
# DimPlot(p_Blastoid, reduction = "tsne", group.by = "samples", cols = cccol)
# dev.off()

###################################################################################################
################################                tsne               ################################
###################################################################################################
p_Blastoid_tsne <- FindNeighbors(object = p_Blastoid, reduction = "tsne", dims = 1:2, k.param = 30)
p_Blastoid_tsne <- FindClusters(object = p_Blastoid_tsne, resolution = 0.1)
clusters <- paste("C",as.numeric(p_Blastoid_tsne@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_tsne@assays$RNA[,])
write.table(clusters,file="p_Blastoid_tsne_ClusteringResult.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_tsne_DimensinalPlot_ColoredByClusters.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_tsne, reduction = "tsne", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_tsne) <- "RNA"
pdf("p_Blastoid_tsne_DimensinalPlot_ColoredByEPI.pdf",width=5.4,height=5)
FeaturePlot(p_Blastoid_tsne, features = EPI[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_Blastoid_tsne, features = EPI[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_Blastoid_tsne, features = EPI, cols = c("lightgrey",cccol[2]))
dev.off()
pdf("p_Blastoid_tsne_DimensinalPlot_ColoredByTE.pdf",width=5.4,height=5)
FeaturePlot(p_Blastoid_tsne, features = TE[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_Blastoid_tsne, features = TE[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_Blastoid_tsne, features = TE, cols = c("lightgrey",cccol[2]))
dev.off()
pdf("p_Blastoid_tsne_DimensinalPlot_ColoredByPE.pdf",width=5.4,height=5)
FeaturePlot(p_Blastoid_tsne, features = PE[1], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
FeaturePlot(p_Blastoid_tsne, features = PE[2], min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all")
DotPlot(object = p_Blastoid_tsne, features = PE, cols = c("lightgrey",cccol[2]))
dev.off()

EPI_score <- apply(p_Blastoid_tsne@assays$RNA[EPI,], 2, mean, na.rm=T)
TE_score <- apply(p_Blastoid_tsne@assays$RNA[TE,], 2, mean, na.rm=T)
PE_score <- apply(p_Blastoid_tsne@assays$RNA[PE,], 2, mean, na.rm=T)
p_Blastoid_tsne@meta.data$EPI_score <- EPI_score
p_Blastoid_tsne@meta.data$TE_score <- TE_score
p_Blastoid_tsne@meta.data$PE_score <- PE_score
pdf("p_Blastoid_tsne_DimensinalPlot_ColoredByCellTypeFeatures.pdf",width=5,height=5)
plot(FeaturePlot(p_Blastoid_tsne, features = "EPI_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
plot(FeaturePlot(p_Blastoid_tsne, features = "TE_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
plot(FeaturePlot(p_Blastoid_tsne, features = "PE_score", min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, keep.scale = "all") & NoAxes() & NoLegend())
DotPlot(object = p_Blastoid_tsne, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()


p_Blastoid_tsne_anno <- RenameIdents(p_Blastoid_tsne,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "TE",
                          `5` = "TE",
                          `6` = "PE",
                          `7` = "EPI",
                          `8` = "TE",
                          `9` = "EPI",
                          `10` = "TE",
                          `11` = "EPI_TE",
                          `12` = "TE",
                          `13` = "PE_TE",
                          `14` = "PE",
                          `15` = "Undefined",
                          `16` = "TE")
pdf("p_Blastoid_tsne_DimensinalPlot_ColoredByCellTypes.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_tsne_anno, reduction = "tsne", cols = cccol, label = TRUE, label.color = "white", label.box = TRUE, repel = TRUE)
DimPlot(object = p_Blastoid_tsne_anno, reduction = "tsne", cols = cccol, repel = TRUE)
dev.off()
saveRDS(p_Blastoid_tsne_anno, "p_Blastoid_tsne_IntergratedMatrixAnnotated.rds")

pdf("p_Blastoid_tsne_DotPlot_ColoredByCellTypes.pdf",width=8,height=4)
DotPlot(object = p_Blastoid_tsne_anno, features = c(EPI, TE, PE), cols = c("white",cccol[1]))
dev.off()

anno <- p_Blastoid_tsne_anno@active.ident
write.table(anno,file="p_Blastoid_tsne_anno_ClusteringResult.txt",quote=F,sep="\t",col.names=F,row.names=T)

# markers
DefaultAssay(p_Blastoid_tsne_anno) <- "RNA"
p_Blastoid_tsne_anno_MarkerTest <- FindAllMarkers(object = p_Blastoid_tsne_anno, logfc.threshold = 0, min.pct = 0, return.thresh = 0.01)
# min.pct    only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
# min.diff.pct    only test genes that show a minimum difference in the fraction of detection between the two groups
# logfc.threshold
# return.thresh    Only return markers that have a p-value < return.thresh
write.table(p_Blastoid_tsne_anno_MarkerTest,file="p_Blastoid_tsne_anno_MarkerTest.txt",quote=F,sep="\t",col.names=T,row.names=T)

pct_cutoff <- 0.2
logfc_cutoff <- 0.8
pdf("p_Blastoid_tsne_anno_ClusterMarkers.VolcanoPlot.pdf",width=4*3,height=4*2)
par(mfrow=c(2,3),mar=c(6,6,4,4))
for (i in names(table(p_Blastoid_tsne_anno_MarkerTest[,"cluster"]))){
    tmp_data <- p_Blastoid_tsne_anno_MarkerTest[which(p_Blastoid_tsne_anno_MarkerTest[,"cluster"]==i),]
    x_pctdiff <- tmp_data[,"pct.1"] - tmp_data[,"pct.2"]; names(x_pctdiff) <- tmp_data[,"gene"]
    y_logfc <- tmp_data[,"avg_log2FC"]; names(y_logfc) <- tmp_data[,"gene"]
    up_genes <- tmp_data[,"gene"][which((x_pctdiff > pct_cutoff) & (y_logfc > logfc_cutoff))]
    dn_genes <- tmp_data[,"gene"][which((x_pctdiff < -pct_cutoff) & (y_logfc < -logfc_cutoff))]
    ns_genes <- setdiff(tmp_data[,"gene"], c(up_genes, dn_genes))

    plot(x_pctdiff, y_logfc, main=paste("",i,sep=""), xlab="Delta Percentage Diff (%)", ylab="log (Fold Change)", type="p", pch=".", col=cccol[8],xlim=c(-1,1),ylim=c(-3,4)); box(lwd=2)
    abline(v=0, h=0, lty=2)
    if (length(up_genes)>0){
        points(x_pctdiff[up_genes], y_logfc[up_genes], pch=16, col=cccol[1])
        text(x_pctdiff[up_genes]+0.15, y_logfc[up_genes], labels=up_genes, cex=0.5)
    }
    if (length(dn_genes)>0){
        points(x_pctdiff[dn_genes], y_logfc[dn_genes], pch=16, col=cccol[2])
        text(x_pctdiff[dn_genes]+0.15, y_logfc[dn_genes], labels=dn_genes, cex=0.5)
    }
}

par(mfrow=c(2,3),mar=c(6,6,4,4))
for (i in names(table(p_Blastoid_tsne_anno_MarkerTest[,"cluster"]))){
    tmp_data <- p_Blastoid_tsne_anno_MarkerTest[which(p_Blastoid_tsne_anno_MarkerTest[,"cluster"]==i),]
    x_pctdiff <- tmp_data[,"pct.1"] - tmp_data[,"pct.2"]; names(x_pctdiff) <- tmp_data[,"gene"]
    y_logfc <- tmp_data[,"avg_log2FC"]; names(y_logfc) <- tmp_data[,"gene"]
    up_genes <- tmp_data[,"gene"][which((x_pctdiff > pct_cutoff) & (y_logfc > logfc_cutoff))]
    dn_genes <- tmp_data[,"gene"][which((x_pctdiff < -pct_cutoff) & (y_logfc < -logfc_cutoff))]
    ns_genes <- setdiff(tmp_data[,"gene"], c(up_genes, dn_genes))

    plot(x_pctdiff, y_logfc, main=paste("",i,sep=""), xlab="Delta Percentage Diff (%)", ylab="log (Fold Change)", type="p", pch=".", col=cccol[8]); box(lwd=2)
    abline(v=0, h=0, lty=2)
    if (length(up_genes)>0){
        points(x_pctdiff[up_genes], y_logfc[up_genes], pch=16, col=cccol[1])
        # text(x_pctdiff[up_genes]+0.15, y_logfc[up_genes], labels=up_genes, cex=0.5)
    }
    if (length(dn_genes)>0){
        points(x_pctdiff[dn_genes], y_logfc[dn_genes], pch=16, col=cccol[2])
        # text(x_pctdiff[dn_genes]+0.15, y_logfc[dn_genes], labels=dn_genes, cex=0.5)
    }
}
dev.off()

p_Blastoid_tsne_anno_Markers <- FindAllMarkers(object = p_Blastoid_tsne_anno, min.pct = 0, logfc.threshold = logfc_cutoff, min.diff.pct = pct_cutoff, return.thresh = 0.01)
write.table(p_Blastoid_tsne_anno_Markers,file="p_Blastoid_tsne_anno_Markers.txt",quote=F,sep="\t",col.names=T,row.names=T)
for (i in names(table(p_Blastoid_tsne_anno_Markers[,"cluster"]))){
    tmp_data <- p_Blastoid_tsne_anno_Markers[which(p_Blastoid_tsne_anno_Markers[,"cluster"]==i),]
    write.table(tmp_data[,"gene"],file=paste("p_Blastoid_tsne_anno_Markers_",i,".txt",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
}
cluster_markers <- c()
for (i in names(table(p_Blastoid_tsne_anno_MarkerTest[,"cluster"]))){
    # cluster_markers <- c(cluster_markers, p_Blastoid_tsne_anno_Markers[which(p_Blastoid_tsne_anno_Markers[,"cluster"]==i & p_Blastoid_tsne_anno_Markers[,"avg_log2FC"] > 0),"gene"][1])
    tmp_matrix <- p_Blastoid_tsne_anno_Markers[which(p_Blastoid_tsne_anno_Markers[,"cluster"]==i),]
    cluster_markers <- c(cluster_markers, tmp_matrix[order(tmp_matrix[,"pct.1"] - tmp_matrix[,"pct.2"], decreasing=T)[1],"gene"])
}

pdf("p_Blastoid_tsne_anno_DimensinalPlot_ColoredByClusterMarkers.pdf",width=5*3,height=5*2)
plot(FeaturePlot(p_Blastoid_tsne_anno, features = cluster_markers, min.cutoff = "q1", max.cutoff = "q99", col = c("#DBD8D650","#218C8C","#5DC863","#F7E61F"), reduction = "tsne", pt.size = 0.1, ncol = 3, keep.scale = "all") & NoAxes() & NoLegend())
DotPlot(object = p_Blastoid_tsne_anno, features = cluster_markers, cols = c(cccol[2],cccol[1]), dot.scale = 9)
dev.off()
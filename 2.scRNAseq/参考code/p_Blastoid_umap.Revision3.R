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
library(gplots)

# plot_col <- cccol[c(1,3,6,8)]
# names(plot_col) <- c("TE","EPI","PE","other")

jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
p_Blastoid_umap <- readRDS("p_Blastoid_umap_IntergratedMatrixAnnotated.rds")

# marker genes
all_genes <- row.names(p_Blastoid_umap)
EPI <- c("SOX2", "POU5F1")
TE <- c("GATA3", "GATA2")
PE <- c("GATA6", "GATA4")
EPI <- intersect(all_genes,EPI)
TE <- intersect(all_genes,TE)
PE <- intersect(all_genes,PE)

EPI_score <- apply(p_Blastoid_umap@assays$RNA[EPI,], 2, mean, na.rm=T)
TE_score <- apply(p_Blastoid_umap@assays$RNA[TE,], 2, mean, na.rm=T)
PE_score <- apply(p_Blastoid_umap@assays$RNA[PE,], 2, mean, na.rm=T)
p_Blastoid_umap@meta.data$EPI_score <- EPI_score
p_Blastoid_umap@meta.data$TE_score <- TE_score
p_Blastoid_umap@meta.data$PE_score <- PE_score

pct_exp_cutoff <- 50
avg_exp_scaled_cutoff <- 0.10

CellTypeAssignment <- function(ratio, avg_exp, name){
    if (sum(ratio) == 0){return("other")}
    if (sum(ratio) == 1){return(name[which(ratio==TRUE)])}
    if (sum(ratio) > 1){return(name[which.max(avg_exp)])}
}


# ###################################################################################################
# ###########################        k.param 10-50    resolution 0.05~0.25       ####################
# ###################################################################################################

# # for (each_k in seq(10,50,5)){
# for (each_k in c(70)){
#     for (each_r in seq(0.05,0.25,0.01)){
#         DefaultAssay(p_Blastoid_umap) <- "integrated"
#         p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = each_k)
#         p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = each_r)
#         pdf(paste("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.",each_r,".",each_k,".pdf",sep=""),width=5.6,height=5)
#         plot(DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE))
#         dev.off()
#     }
# }

# for (each_s in seq(1,199,18)){
#     for (each_k in c(30,40,50)){
#         for (each_r in c(0.1,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.2)){
#             DefaultAssay(p_Blastoid_umap) <- "integrated"
#             p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = each_k)
#             p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = each_r, random.seed = each_s)
#             pdf(paste("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.",each_s,".",each_r,".",each_k,".pdf",sep=""),width=5.6,height=5)
#             plot(DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE))
#             dev.off()
#         }
#     }
# }

# for (each_s in c(seq(2,17,1),seq(20,36,1),seq(38,54,1))){
#     for (each_k in c(30,40,50)){
#         for (each_r in c(0.12)){
#             DefaultAssay(p_Blastoid_umap) <- "integrated"
#             p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = each_k)
#             p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = each_r, random.seed = each_s)
#             pdf(paste("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.",each_s,".",each_r,".",each_k,".pdf",sep=""),width=5.6,height=5)
#             plot(DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE))
#             dev.off()
#         }
#     }
# }
###################################################################################################
###########################        k.param 30    resolution 0.15       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 30)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.15)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.15.30.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.15.30.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.15.30.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.15.30.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.15.30.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.15.30.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "TE",
                          `5` = "EPI",
                          `6` = "TE",
                          `7` = "PE",
                          `8` = "other",
                          `9` = "other",
                          `10` = "PE",
                          `11` = "other",
                          `12` = "TE")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.15.30.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.15.30.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 40    resolution 0.15       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 40)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.15, random.seed = 145)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.15.40.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.15.40.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.15.40.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.15.40.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.15.40.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.15.40.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "TE",
                          `6` = "PE",
                          `7` = "TE",
                          `8` = "other",
                          `9` = "other",
                          `10` = "PE",
                          `11` = "other")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.15.40.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.15.40.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 50    resolution 0.15       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 50)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.15, random.seed = 145)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.15.50.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.15.50.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.15.50.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.15.50.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.15.50.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.15.50.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "PE",
                          `6` = "TE",
                          `7` = "other",
                          `8` = "other",
                          `9` = "PE",
                          `10` = "other")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.15.50.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.15.50.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 30    resolution 0.12       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 30)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.12, random.seed = 109)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.12.30.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.12.30.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.12.30.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.12.30.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.12.30.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.12.30.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "PE",
                          `6` = "TE",
                          `7` = "other",
                          `8` = "other",
                          `9` = "PE",
                          `10` = "other",
                          `11` = "TE")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.12.30.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.12.30.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 40    resolution 0.12       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 40)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.12, random.seed = 109)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.12.40.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.12.40.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.12.40.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.12.40.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.12.40.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.12.40.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "PE",
                          `6` = "other",
                          `7` = "other",
                          `8` = "PE",
                          `9` = "other")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.12.40.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.12.40.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 50    resolution 0.12       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 50)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.12, random.seed = 109)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.12.50.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.12.50.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.12.50.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.12.50.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.12.50.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.12.50.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "PE",
                          `6` = "other",
                          `7` = "PE",
                          `8` = "other")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.12.50.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.12.50.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 30    resolution 0.18       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 30)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.18, random.seed = 109)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.18.30.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.18.30.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.18.30.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.18.30.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.18.30.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.18.30.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "TE",
                          `5` = "TE",
                          `6` = "EPI",
                          `7` = "TE",
                          `8` = "PE",
                          `9` = "EPI",
                          `10` = "TE",
                          `11` = "other",
                          `12` = "TE",
                          `13` = "PE",
                          `14` = "other",
                          `15` = "PE",
                          `16` = "other",
                          `17` = "TE")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.18.30.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.18.30.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 40    resolution 0.18       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 40)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.18, random.seed = 109)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.18.40.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.18.40.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.18.40.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.18.40.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.18.40.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.18.40.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "TE",
                          `6` = "PE",
                          `7` = "TE",
                          `8` = "other",
                          `9` = "other",
                          `10` = "PE",
                          `11` = "other")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.18.40.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.18.40.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 50    resolution 0.18       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 50)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.18, random.seed = 109)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.18.50.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.18.50.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.18.50.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.18.50.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.18.50.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.18.50.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "TE",
                          `6` = "PE",
                          `7` = "TE",
                          `8` = "other",
                          `9` = "other",
                          `10` = "PE",
                          `11` = "other")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.18.50.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.18.50.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 30    resolution 0.1       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 30)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.1, random.seed = 109)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.1.30.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.1.30.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.1.30.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.1.30.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.1.30.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.1.30.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "PE",
                          `6` = "other",
                          `7` = "other",
                          `8` = "PE",
                          `9` = "other")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.1.30.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.1.30.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 40    resolution 0.1       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 40)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.1, random.seed = 109)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.1.40.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.1.40.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.1.40.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.1.40.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.1.40.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.1.40.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "PE",
                          `6` = "other",
                          `7` = "PE",
                          `8` = "other")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.1.40.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.1.40.txt",quote=F,sep="\t",col.names=F,row.names=T)


###################################################################################################
###########################        k.param 50    resolution 0.1       ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 50)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.1, random.seed = 109)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.1.50.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.1.50.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.1.50.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.1.50.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.1.50.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.1.50.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "EPI",
                          `4` = "PE",
                          `5` = "TE",
                          `6` = "other",
                          `7` = "other")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.1.50.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.1.50.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 30    resolution 0.2        ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 30)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.2, random.seed = 127)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.2.30.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.2.30.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.2.30.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.2.30.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.2.30.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.2.30.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "PE",
                          `5` = "TE",
                          `6` = "TE",
                          `7` = "EPI",
                          `8` = "TE",
                          `9` = "EPI",
                          `10` = "TE",
                          `11` = "other",
                          `12` = "other",
                          `13` = "PE",
                          `14` = "other",
                          `15` = "TE")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.2.30.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,6,3,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.2.30.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 40    resolution 0.2        ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 40)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.2, random.seed = 127)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.2.40.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.2.40.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.2.40.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.2.40.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.2.40.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.2.40.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "TE",
                          `6` = "TE",
                          `7` = "PE",
                          `8` = "TE",
                          `9` = "other",
                          `10` = "PE",
                          `11` = "other",
                          `12` = "PE",
                          `13` = "other",
                          `14` = "TE")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.2.40.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.2.40.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################        k.param 50    resolution 0.2        ############################
###################################################################################################
DefaultAssay(p_Blastoid_umap) <- "integrated"
p_Blastoid_umap <- FindNeighbors(object = p_Blastoid_umap, reduction = "umap", dims = 1:2, k.param = 50)
p_Blastoid_umap <- FindClusters(object = p_Blastoid_umap, resolution = 0.2, random.seed = 127)
clusters <- paste("C",as.numeric(p_Blastoid_umap@meta.data$seurat_clusters)-1,sep="")
names(clusters) <- colnames(p_Blastoid_umap@assays$RNA[,])
write.table(clusters,file="p_Blastoid_umap_ClusteringResult.0.2.50.txt",quote=F,sep="\t",col.names=F,row.names=T)

pdf("p_Blastoid_umap_DimensinalPlot_ColoredByClusters.0.2.50.pdf",width=5.6,height=5)
DimPlot(object = p_Blastoid_umap, reduction = "umap", cols = cccol, label = TRUE, repel = TRUE)
dev.off()

DefaultAssay(p_Blastoid_umap) <- "RNA"
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypeFeatures.0.2.50.pdf",width=5,height=5)
DotPlot(object = p_Blastoid_umap, features = c("EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
dev.off()

dp <- DotPlot(object = p_Blastoid_umap, features = c(PE, EPI, TE,"EPI_score","TE_score","PE_score"), cols = c("lightgrey",cccol[2]))
pct_exp <- matrix(dp$data$pct.exp, c(length(dp$data$pct.exp)/9,9), byrow=T)
avg_exp_scaled <- matrix(dp$data$avg.exp.scaled, c(length(dp$data$pct.exp)/9,9), byrow=T)
colnames(pct_exp) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(pct_exp) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(pct_exp, file="pct_exp.0.2.50.txt", quote=F, sep="\t", row.names=T, col.names=T)
colnames(avg_exp_scaled) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(avg_exp_scaled) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")
write.table(avg_exp_scaled, file="avg_exp_scaled.0.2.50.txt", quote=F, sep="\t", row.names=T, col.names=T)

rationale <- pct_exp > pct_exp_cutoff & (avg_exp_scaled > avg_exp_scaled_cutoff)
colnames(rationale) <- c(PE, EPI, TE,"EPI_score","TE_score","PE_score")
rownames(rationale) <- paste("c", seq(length(dp$data$pct.exp)/9)-1, sep="")

celltype <- c()
for (i in seq(nrow(rationale))){
    celltype <- c(celltype, CellTypeAssignment(rationale[i, c("EPI_score","TE_score","PE_score")], avg_exp_scaled[i, c("EPI_score","TE_score","PE_score")], c("EPI","TE","PE")))
}
rationale <- cbind(celltype, rationale)
colnames(rationale) <- c("AssignCellType",PE, EPI, TE,"EPI_score","TE_score","PE_score")
write.table(rationale, file="rationale.0.2.50.txt", quote=F, sep="\t", row.names=T, col.names=T)
p_Blastoid_umap_tmp <- RenameIdents(p_Blastoid_umap,
                          `0` = "TE",
                          `1` = "TE",
                          `2` = "TE",
                          `3` = "TE",
                          `4` = "EPI",
                          `5` = "TE",
                          `6` = "TE",
                          `7` = "PE",
                          `8` = "other",
                          `9` = "other",
                          `10` = "PE",
                          `11` = "other")
pdf("p_Blastoid_umap_DimensinalPlot_ColoredByCellTypes.0.2.50.pdf",width=6.2,height=5)
DimPlot(object = p_Blastoid_umap_tmp, reduction = "umap", cols = cccol[c(1,3,6,8)], repel = TRUE)
dev.off()
anno <- p_Blastoid_umap_tmp@active.ident
write.table(anno,file="p_Blastoid_umap_anno_ClusteringResult.0.2.50.txt",quote=F,sep="\t",col.names=F,row.names=T)

###################################################################################################
###########################           Jaccard index calculation        ############################
###################################################################################################
cr_01_30 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.1.30.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.1.30.txt")[,2]),sep="_")
cr_01_40 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.1.40.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.1.40.txt")[,2]),sep="_")
cr_01_50 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.1.50.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.1.50.txt")[,2]),sep="_")
cr_012_30 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.12.30.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.12.30.txt")[,2]),sep="_")
cr_012_40 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.12.40.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.12.40.txt")[,2]),sep="_")
cr_012_50 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.12.50.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.12.50.txt")[,2]),sep="_")
cr_015_30 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.15.30.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.15.30.txt")[,2]),sep="_")
cr_015_40 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.15.40.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.15.40.txt")[,2]),sep="_")
cr_015_50 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.15.50.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.15.50.txt")[,2]),sep="_")
cr_018_30 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.18.30.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.18.30.txt")[,2]),sep="_")
cr_018_40 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.18.40.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.18.40.txt")[,2]),sep="_")
cr_018_50 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.18.50.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.18.50.txt")[,2]),sep="_")
cr_02_30 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.2.30.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.2.30.txt")[,2]),sep="_")
cr_02_40 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.2.40.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.2.40.txt")[,2]),sep="_")
cr_02_50 <- paste(as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.2.50.txt")[,1]),as.vector(read.table("p_Blastoid_umap_anno_ClusteringResult.0.2.50.txt")[,2]),sep="_")

data_matrix <- cbind(cr_01_30,cr_01_40,cr_01_50,cr_012_30,cr_012_40,cr_012_50,cr_015_30,cr_015_40,cr_015_50,cr_018_30,cr_018_40,cr_018_50,cr_02_30,cr_02_40,cr_02_50)
colnames(data_matrix) <- c("cr_01_30","cr_01_40","cr_01_50","cr_012_30","cr_012_40","cr_012_50","cr_015_30","cr_015_40","cr_015_50","cr_018_30","cr_018_40","cr_018_50","cr_02_30","cr_02_40","cr_02_50")
jaccindex <- c()
for (i in colnames(data_matrix)){
    tmp_ji <- c()
    for (j in colnames(data_matrix)){
        # print(c(i,j,jaccard(data_matrix[,i], data_matrix[,j])))
        tmp_ji <- c(tmp_ji, jaccard(data_matrix[,i], data_matrix[,j]))
    }
    jaccindex <- rbind(jaccindex, tmp_ji)
}
colnames(jaccindex) <- c("res=0.1,k=30","res=0.1,k=40","res=0.1,k=50","res=0.12,k=30","res=0.12,k=40","res=0.12,k=50","res=0.15,k=30","res=0.15,k=40","res=0.15,k=50","res=0.18,k=30","res=0.18,k=40","res=0.18,k=50","res=0.2,k=30","res=0.2,k=40","res=0.2,k=50")
rownames(jaccindex) <- c("res=0.1,k=30","res=0.1,k=40","res=0.1,k=50","res=0.12,k=30","res=0.12,k=40","res=0.12,k=50","res=0.15,k=30","res=0.15,k=40","res=0.15,k=50","res=0.18,k=30","res=0.18,k=40","res=0.18,k=50","res=0.2,k=30","res=0.2,k=40","res=0.2,k=50")


pdf("JaccardIndex.Heatmap.pdf",width=6,height=5)
plot_data <- jaccindex
layout(matrix(c(rep(1,12),rep(2,3)),nrow=1,ncol=15))
all_val <- as.matrix(plot_data) # using same scale bar
zmax <- max(na.omit(all_val))
zmin <- min(na.omit(all_val))
ColorRamp <- colorRampPalette(c("#605571","#3A518A","#218C8C","#5DC863","#F7E61F"), bias=1)(10000)
# ColorRamp <- colorRampPalette(c("lightblue","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence

tmp_matrix <- plot_data
tmp_min <- as.numeric(max(zmin,min(tmp_matrix,na.rm=T)))
tmp_max <- as.numeric(min(zmax,max(tmp_matrix,na.rm=T)))
tmp_matrix[tmp_matrix<tmp_min] <- tmp_min
tmp_matrix[tmp_matrix>tmp_max] <- tmp_max
tmp_LeftRange <- round((10000/(zmax-zmin))*(tmp_min-zmin))
tmp_RightRange <- round((10000/(zmax-zmin))*(tmp_max-zmin))
par(mar=c(6,6,2,2))
image(1:ncol(tmp_matrix), 1:nrow(tmp_matrix), t(tmp_matrix), axes=F, col=ColorRamp[tmp_LeftRange:tmp_RightRange], xlab="", main="JGES",ylab="");box(lwd=2)
axis(side=1,1:ncol(tmp_matrix),labels=colnames(data_matrix),las=2)
axis(side=2,1:nrow(tmp_matrix),labels=colnames(data_matrix),las=2)


par(mar=c(6,4,2,2))
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="log2(fpkm+1)",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=2,c(zmin,round((zmax+zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()

pdf("JaccardIndex.hCluster.Heatmap.pdf",width=5,height=5)
heatmap.2(jaccindex,main="",col=ColorRamp,key=F,trace="none",margins=c(8,8), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
dev.off()
pdf("JaccardIndex.hCluster.Heatmap.legend.pdf",width=3,height=1.8)
par(mar=c(4,2,4,2))
plotMatrix <- jaccindex
all_val <- as.matrix(plotMatrix) # using same scale bar
zmax <- max(na.omit(all_val))
zmin <- min(na.omit(all_val))
ColorLevels <- seq(to=zmax,from=zmin, length=1000)   #number sequence
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="Jaccard Index",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,3),round((zmax+zmin)/2,3),round(zmax,3)))
dev.off()


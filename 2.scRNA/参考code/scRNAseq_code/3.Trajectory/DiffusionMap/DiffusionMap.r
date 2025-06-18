###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#f7acbc","#deab8a","#817936","#444693","#ef5b9c","#fedcbd","#7f7522","#2b4490","#f47920","#80752c","#2a5caa","#f05b72","#905a3d","#87843b","#224b8f","#f15b6c","#8f4b2e","#726930","#003a6c","#f8aba6","#87481f","#454926","#102b6a","#f69c9f","#5f3c23","#2e3a1f","#426ab3","#f58f98","#6b473c","#4d4f36","#46485f","#ca8687","#faa755","#b7ba6b","#4e72b8","#f391a9","#fab27b","#b2d235","#181d4b","#bd6758","#f58220","#5c7a29","#1a2933","#d71345","#843900","#bed742","#121a2a","#d64f44","#905d1d","#7fb80e","#0c212b","#d93a49","#8a5d19","#a3cf62","#6a6da9","#b3424a","#8c531b","#769149","#585eaa","#c76968","#826858","#6d8346","#494e8f","#bb505d","#64492b","#78a355","#afb4db","#987165","#ae6642","#abc88b","#9b95c9","#ac6767","#56452d","#74905d","#6950a1","#973c3f","#96582a","#cde6c7","#6f60aa","#b22c46","#705628","#1d953f","#867892","#a7324a","#4a3113","#77ac98","#918597","#aa363d","#412f1f","#007d65","#6f6d85","#ed1941","#845538","#84bf96","#594c6d","#f26522","#8e7437","#45b97c","#694d9f","#d2553d","#69541b","#225a1f","#6f599c","#b4534b","#d5c59f","#367459","#8552a1","#ef4136","#cd9a5b","#007947","#543044","#c63c26","#cd9a5b","#40835e","#63434f","#f3715c","#b36d41","#2b6447","#7d5886","#a7573b","#df9464","#005831","#401c44","#aa2116","#b76f40","#006c54","#472d56","#b64533","#ad8b3d","#375830","#45224a","#b54334","#dea32c","#274d3d","#411445","#853f04","#d1923f","#375830","#4b2f3d","#840228","#c88400","#27342b","#402e4c","#7a1723","#c37e00","#65c294","#c77eb5","#a03939","#c37e00","#73b9a2","#ea66a6","#8a2e3b","#e0861a","#72baa7","#f173ac","#8e453f","#ffce7b","#005344","#8f4b4a","#fcaf17","#122e29","#892f1b","#ba8448","#293047","#f6f5ec","#6b2c25","#896a45","#00ae9d","#d9d6c3","#733a31","#76624c","#508a88","#d1c7b7","#54211d","#6d5826","#70a19f","#f2eada","#78331e","#ffc20e","#50b7c1","#d3d7d4","#53261f","#fdb933","#00a6ac","#999d9c","#f15a22","#d3c6a6","#78cdd1","#a1a3a6","#b4533c","#c7a252","#008792","#9d9087","#84331f","#dec674","#94d6da","#8a8c8e","#f47a55","#b69968","#afdfe4","#74787c","#f15a22","#c1a173","#5e7c85","#7c8577","#f3704b","#dbce8f","#76becc","#72777b","#da765b","#ffd400","#90d7ec","#77787b","#c85d44","#ffd400","#009ad6","#4f5555","#ae5039","#ffe600","#145b7d","#6c4c49","#6a3427","#f0dc70","#11264f","#563624","#8f4b38","#fcf16e","#7bbfea","#3e4145","#8e3e1f","#decb00","#33a3dc","#3c3645","#f36c21","#cbc547","#228fbd","#464547","#b4532a","#6e6b41","#2468a2","#130c0e","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")

cccol01 <- paste(cccol,"01",sep="")
cccol05 <- paste(cccol,"05",sep="")
cccol30 <- paste(cccol,"30",sep="")
cccol50 <- paste(cccol,"50",sep="")
cccol80 <- paste(cccol,"80",sep="")

ccpch <- c(3,4,7,8,15,17,18,19)
red_col <- colorRampPalette(c(cccol[1],cccol[4]), bias=1)(16)
blue_col <- colorRampPalette(c(cccol[2],cccol[4]), bias=1)(16)
grey_col <- colorRampPalette(c(cccol[8],"black"), bias=1)(16)

# options(scipen = 200)

library("Seurat")
# library("diffusionMap")
# library("devtools")
# install_github("theislab/destiny")
library("destiny")
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
PGC_KL <- readRDS("PGC.rds")
PGC_XM <- readRDS("cPGC_WG_PGC.rds")

PGC_KL <- NormalizeData(object = PGC_KL, normalization.method = "LogNormalize")
PGC_XM <- NormalizeData(object = PGC_XM, normalization.method = "LogNormalize")
PGC_KL <- FindVariableFeatures(object = PGC_KL, selection.method = "vst", nfeatures = 1200)
PGC_XM <- FindVariableFeatures(object = PGC_XM, selection.method = "vst", nfeatures = 1200)

all_genes <- intersect(row.names(PGC_KL),row.names(PGC_XM))

# identify anchors
chk_list <- list(PGC_KL, PGC_XM)
chk_features <- SelectIntegrationFeatures(object.list = chk_list, nfeatures = 1200)
chk_anchors <- FindIntegrationAnchors(object.list = chk_list, anchor.features = chk_features)

# integrated (or ‘batch-corrected’)
data_integrated <- IntegrateData(anchorset = chk_anchors, normalization.method = "LogNormalize", features.to.integrate = all_genes)

samples <- c(rep("PGC_KL",ncol(PGC_KL)),rep("PGC_XM",ncol(PGC_XM)))
clusters <- c(PGC_KL@meta.data$seurat_clusters, PGC_XM@meta.data$seurat_clusters)
stages <- c(PGC_KL@meta.data$stages, PGC_XM@meta.data$stages)

data_integrated@meta.data$samples <- samples
data_integrated@meta.data$clusters <- clusters
data_integrated@meta.data$stages <- stages

# saveRDS(data_integrated, "PGC_IntergratedMatrix.rds")
# write.table(file="PGC_IntergratedMatrix.txt",data_integrated@assays$integrated[,],sep="\t",col.names=T,row.names=T,quote=F)

###################################################################################################
################################                PLOT               ################################
###################################################################################################
# COLS1 <- cccol50[clusters]
COLS1 <- rep(1,length(clusters))
COLS1[which(clusters==13)] <- cccol80[13]
COLS1[which(clusters==14)] <- cccol80[14]
COLS1[which(clusters==15)] <- cccol80[15]
COLS1[which(clusters==24)] <- cccol80[24]
COLS1[which(clusters==26)] <- cccol80[26]
COLS1[which(clusters==18)] <- cccol[1]

COLS2 <- rep(1,length(stages)); names(COLS2) <- c(colnames(PGC_KL),colnames(PGC_XM))
COLS2[which(stages=="01.EGKX_0h_WE")] <- red_col[1]
COLS2[which(stages=="02.EGKXI_1h_WE")] <- red_col[2]
COLS2[which(stages=="03.EGKXI_2h_WE")] <- red_col[3]
COLS2[which(stages=="04.EGKXII_3h_WE")] <- red_col[4]
COLS2[which(stages=="05.EGKXII_4h_WE")] <- red_col[5]
COLS2[which(stages=="06.EGKXIII_5h_WE")] <- red_col[6]
COLS2[which(stages=="07.HH2_6h_WE")] <- red_col[7]
COLS2[which(stages=="08.HH2a_10h_WE")] <- red_col[8]
COLS2[which(stages=="09.HH3_13h_WE")] <- red_col[9]
COLS2[which(stages=="10.HH3a_16h_WE")] <- red_col[10]
COLS2[which(stages=="11.HH4_19h_WE")] <- red_col[11]
COLS2[which(stages=="12.HH5_22h_WE")] <- red_col[12]
COLS2[which(stages=="13.HH7_23h_WE")] <- red_col[13]
COLS2[which(stages=="14.HH8_26h_WE")] <- red_col[14]
COLS2[which(stages=="15.HH9_29h_WE")] <- red_col[15]
COLS2[which(stages=="16.HH10_33h_WE")] <- red_col[16]
COLS2[which(clusters==13)] <- cccol80[13]
COLS2[which(clusters==14)] <- cccol80[14]
COLS2[which(clusters==15)] <- cccol80[15]
COLS2[which(clusters==24)] <- cccol80[24]
COLS2[which(clusters==26)] <- cccol80[26]

# before intergrated
data_matrix <- t(cbind(as.matrix(PGC_KL@assays$RNA[all_genes,]),as.matrix(PGC_XM@assays$RNA[all_genes,])))

sigmas <- find_sigmas(data_matrix, verbose = FALSE)
dm1 <- DiffusionMap(data_matrix, sigmas, verbose = FALSE, n_pcs = 50, n_eigs = 10)
pdf("DiffusionMap.BeforeIntergrated.pdf",width=6,height=6)
plot(sigmas)
plot(dm1@eigenvectors[,"DC1"],dm1@eigenvectors[,"DC2"], col=COLS1, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC1"],dm1@eigenvectors[,"DC3"], col=COLS1, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC2"],dm1@eigenvectors[,"DC3"], col=COLS1, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC1"],dm1@eigenvectors[,"DC2"], col=COLS2, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC1"],dm1@eigenvectors[,"DC3"], col=COLS2, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC2"],dm1@eigenvectors[,"DC3"], col=COLS2, pch=16, cex=0.3)
dev.off()


# after intergrated
data_matrix <- t(as.matrix(data_integrated@assays$integrated[all_genes,]))

sigmas <- find_sigmas(data_matrix, verbose = FALSE)
dm1 <- DiffusionMap(data_matrix, sigmas, verbose = FALSE, n_pcs = 50, n_eigs = 10)
pdf("DiffusionMap.AfterIntergrated.pdf",width=6,height=6)
plot(sigmas)
plot(dm1, col=COLS1)
plot(dm1, col=COLS2)
plot(dm1@eigenvectors[,"DC1"],dm1@eigenvectors[,"DC2"], col=COLS1, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC1"],dm1@eigenvectors[,"DC3"], col=COLS1, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC2"],dm1@eigenvectors[,"DC3"], col=COLS1, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC1"],dm1@eigenvectors[,"DC2"], col=COLS2, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC1"],dm1@eigenvectors[,"DC3"], col=COLS2, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC2"],dm1@eigenvectors[,"DC3"], col=COLS2, pch=16, cex=0.3)
dev.off()

pdf("DiffusionMap.legend.pdf",width=6,height=6)
plot(1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("center",c("C12","C13","C14","C23","C25"),col=cccol80[c(13,14,15,24,26)],pch=16,bty="n")
plot(1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("center",c("01.EGKX_0h_WE","02.EGKXI_1h_WE","03.EGKXI_2h_WE","04.EGKXII_3h_WE","05.EGKXII_4h_WE","06.EGKXIII_5h_WE","07.HH2_6h_WE","08.HH2a_10h_WE","09.HH3_13h_WE","10.HH3a_16h_WE","11.HH4_19h_WE","12.HH5_22h_WE","13.HH7_23h_WE","14.HH8_26h_WE","15.HH9_29h_WE","16.HH10_33h_WE"),col=red_col,ncol=2,pch=16,bty="n")
dev.off()

# only late data
clusters_XM <- c(PGC_XM@meta.data$seurat_clusters)
COLS_XM <- rep(1,length(clusters_XM))
COLS_XM[which(clusters_XM==13)] <- cccol80[13]
COLS_XM[which(clusters_XM==14)] <- cccol80[14]
COLS_XM[which(clusters_XM==15)] <- cccol80[15]
COLS_XM[which(clusters_XM==24)] <- cccol80[24]
COLS_XM[which(clusters_XM==26)] <- cccol80[26]

data_matrix <- t(as.matrix(PGC_XM@assays$RNA[all_genes,]))

sigmas <- find_sigmas(data_matrix, verbose = FALSE)
dm1 <- DiffusionMap(data_matrix, sigmas, verbose = FALSE, n_pcs = 50, n_eigs = 10)
pdf("DiffusionMap.XM.pdf",width=6,height=6)
plot(sigmas)
plot(dm1, col=COLS_XM)
plot(dm1@eigenvectors[,"DC1"],dm1@eigenvectors[,"DC2"], col=COLS_XM, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC1"],dm1@eigenvectors[,"DC3"], col=COLS_XM, pch=16, cex=0.3)
plot(dm1@eigenvectors[,"DC2"],dm1@eigenvectors[,"DC3"], col=COLS_XM, pch=16, cex=0.3)
dev.off()

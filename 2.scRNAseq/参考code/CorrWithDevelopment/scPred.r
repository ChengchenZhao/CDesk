###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
allcol <- c("#f7acbc","#deab8a","#817936","#444693","#ef5b9c","#fedcbd","#7f7522","#2b4490","#f47920","#80752c","#2a5caa","#f05b72","#905a3d","#87843b","#224b8f","#f15b6c","#8f4b2e","#726930","#003a6c","#f8aba6","#87481f","#454926","#102b6a","#f69c9f","#5f3c23","#2e3a1f","#426ab3","#f58f98","#6b473c","#4d4f36","#46485f","#ca8687","#faa755","#b7ba6b","#4e72b8","#f391a9","#fab27b","#b2d235","#181d4b","#bd6758","#f58220","#5c7a29","#1a2933","#d71345","#843900","#bed742","#121a2a","#d64f44","#905d1d","#7fb80e","#0c212b","#d93a49","#8a5d19","#a3cf62","#6a6da9","#b3424a","#8c531b","#769149","#585eaa","#c76968","#826858","#6d8346","#494e8f","#bb505d","#64492b","#78a355","#afb4db","#987165","#ae6642","#abc88b","#9b95c9","#ac6767","#56452d","#74905d","#6950a1","#973c3f","#96582a","#cde6c7","#6f60aa","#b22c46","#705628","#1d953f","#867892","#a7324a","#4a3113","#77ac98","#918597","#aa363d","#412f1f","#007d65","#6f6d85","#ed1941","#845538","#84bf96","#594c6d","#f26522","#8e7437","#45b97c","#694d9f","#d2553d","#69541b","#225a1f","#6f599c","#b4534b","#d5c59f","#367459","#8552a1","#ef4136","#cd9a5b","#007947","#543044","#c63c26","#cd9a5b","#40835e","#63434f","#f3715c","#b36d41","#2b6447","#7d5886","#a7573b","#df9464","#005831","#401c44","#aa2116","#b76f40","#006c54","#472d56","#b64533","#ad8b3d","#375830","#45224a","#b54334","#dea32c","#274d3d","#411445","#853f04","#d1923f","#375830","#4b2f3d","#840228","#c88400","#27342b","#402e4c","#7a1723","#c37e00","#65c294","#c77eb5","#a03939","#c37e00","#73b9a2","#ea66a6","#8a2e3b","#e0861a","#72baa7","#f173ac","#8e453f","#ffce7b","#005344","#8f4b4a","#fcaf17","#122e29","#892f1b","#ba8448","#293047","#f6f5ec","#6b2c25","#896a45","#00ae9d","#d9d6c3","#733a31","#76624c","#508a88","#d1c7b7","#54211d","#6d5826","#70a19f","#f2eada","#78331e","#ffc20e","#50b7c1","#d3d7d4","#53261f","#fdb933","#00a6ac","#999d9c","#f15a22","#d3c6a6","#78cdd1","#a1a3a6","#b4533c","#c7a252","#008792","#9d9087","#84331f","#dec674","#94d6da","#8a8c8e","#f47a55","#b69968","#afdfe4","#74787c","#f15a22","#c1a173","#5e7c85","#7c8577","#f3704b","#dbce8f","#76becc","#72777b","#da765b","#ffd400","#90d7ec","#77787b","#c85d44","#ffd400","#009ad6","#4f5555","#ae5039","#ffe600","#145b7d","#6c4c49","#6a3427","#f0dc70","#11264f","#563624","#8f4b38","#fcf16e","#7bbfea","#3e4145","#8e3e1f","#decb00","#33a3dc","#3c3645","#f36c21","#cbc547","#228fbd","#464547","#b4532a","#6e6b41","#2468a2","#130c0e","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
cccol05 <- c("#CE001305","#16557A05","#C7A60905","#87C23205","#00879205","#A14C9405","#15A08C05","#8B7E7505","#1E7CAF05","#EA425F05","#46489A05","#E5003305","#0F231F05","#1187CD05")
options(scipen = 5)
red_col <- colorRampPalette(c(cccol[1],cccol[8],cccol[3],cccol[6]), bias=1)(4)
blue_col <- c(colorRampPalette(c(cccol[4],cccol[7]), bias=1)(5),colorRampPalette(c(cccol[2],cccol[11],cccol[9],cccol[14]), bias=1)(5),colorRampPalette(c(cccol[16],cccol[15]), bias=1)(5))
grey_col <- colorRampPalette(c(cccol[8],"black"), bias=1)(10)

library("Seurat")
library("hdf5r")
library("scPred")
library("magrittr")

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################

## ref
# ref <- readRDS("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/Revision2/p_integr_revsub_tsne_anno.rds")
# ref <- readRDS("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/Revision2/p_integr_revsub_umap_anno.rds")

ref <- readRDS("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/Revision/E789/p_integr_dev_sub_IntergratedCluster.rds")
ref <- RenameIdents(ref,
                    `0` = "UD",
                    `1` = "PE",
                    `2` = "TE",
                    `3` = "TE",
                    `4` = "EPI")

# clusters <- paste("C",as.numeric(p_integr_dev_sub@meta.data$seurat_clusters)-1,sep="")
# sub_cell 
# celltypes[names(which(clusters=="C0" | clusters=="C1" | clusters=="C2"))] = "d_PE"
# ref <- subset()

## query
qry <- readRDS("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/p_Blastoid_umap_IntergratedMatrixAnnotated.rds")

pdf("scPred.Clusters.pdf", width=5, height=5)
# DimPlot(object = ref, reduction = "tsne", cols = cccol[c(3,8,1,6)], label = TRUE, repel = TRUE)
DimPlot(object = ref, reduction = "umap", cols = cccol[c(6,3,1)], label = TRUE, repel = TRUE)
DimPlot(object = qry, reduction = "umap", cols = cccol[c(6,1,3,8)], label = TRUE, repel = TRUE)
dev.off()

common_genes <- intersect(rownames(ref),rownames(ref))

ref@meta.data$cell_type <- ref@active.ident
qry@meta.data$cell_type <- qry@active.ident

DefaultAssay(ref) <- "integrated"
DefaultAssay(qry) <- "integrated"

meta_data <- read.table("/mnt/zhaochengchen/Work/OTHER/00.ZCC/1.pIPS/CorrWithDevelopment/meta.txt",header=T,row.names=2,sep="\t")
E7_LB_cells <- intersect(colnames(ref),rownames(meta_data))[which(meta_data[intersect(colnames(ref),rownames(meta_data)),2] == "E7_LB")]
E8_HB_cells <- intersect(colnames(ref),rownames(meta_data))[which(meta_data[intersect(colnames(ref),rownames(meta_data)),2] == "E8_HB")]
E9_EBi_cells <- intersect(colnames(ref),rownames(meta_data))[which(meta_data[intersect(colnames(ref),rownames(meta_data)),2] == "E9_EBi")]

ref_E7 <- subset(ref, cells = E7_LB_cells)
ref_E8 <- subset(ref, cells = E8_HB_cells)
ref_E9 <- subset(ref, cells = E9_EBi_cells)
###################################################################################################
################################             CALU & PLOT           ################################
###################################################################################################
# getFeatureSpace will create a scPred object stored in the @misc slot. This object will contained all required information to classify cells. See ?getFeatureSpace help documentation.
ref <- getFeatureSpace(ref, "cell_type")

# Secondly, we train the classifiers for each cell using the trainModel function. By default, scPred will use a support vector machine with a radial kernel.
# ref <- trainModel(ref)
ref <- trainModel(ref, model = "glm")
# ref <- trainModel(ref, model = "mda")

# Training probabilities for each cell in the reference data can be accessed using the get_probabilities method:
get_probabilities(ref) %>% head()

# We can use the get_scpred method to retrieve the scPred object from the Seurat object. Printing a scPred object will show for each cell type:

# The number of cells
# Number of features used to train the model
# Prediction model (method)
# Performance metrics
get_scpred(ref)

# To visualize the performance for each cell type we can use the plot_probabilities function:
pdf("scPred.Probabilities_svmRadial.pdf", width=5, height=5)
plot_probabilities(ref)
dev.off()

# # A different model can be specified using the model parameter and providing the method value from caret (e.g. mda for a mixture discriminant analysis using the mda package). Additionally, if only an mda model wants to be applied to a subset of cells, we can specify this using the reclassify parameter. In this case, we want to train different models for “cMono” and “ncMono” to improve their classification performance:
# ref <- trainModel(ref, model = "mda", reclassify = c("cMono", "ncMono"))
# # The code above trains a mixture discriminant analysis for two cell types and preserves the previous support vector machines for the remaing cell types.
# # We can observe a change in the sensitivity for “cMono” and “ncMono”:
# get_scpred(ref)

# # and also verify that higher probabilities for these cell types were obtained by plotting the training probabilities again:
# pdf("scPred.Probabilities_svmRadial2.pdf", width=5, height=5)
# plot_probabilities(ref)
# dev.off()


# # First, let’s normalize the query dataset (cells to be classfied).
# qry <- NormalizeData(qry)
# Finally, we ca classify the cells from the query data using the scPredict function. The first argument corresponds to the query object and the second to the reference object (with a scPred model trained already).

# scPred now uses Harmony to align the query data onto the training low-dimensional space used as reference. Once the data is aligned, cells are classified using the pre-trained models.
# scPredict will return the query dataset. Make sure the left-side value of the <- operator corresponds to the query data.
qry <- scPredict(qry, ref, recompute_alignment=T, threshold=0.5)

pdf("scPred.scpredPrediction.pdf", width=5, height=5)
DimPlot(qry, group.by = "scpred_prediction", reduction = "scpred", label = TRUE, repel = TRUE)
DimPlot(qry, group.by = "scpred_prediction", reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

# We can also run UMAP using the aligned data as an input
qry <- RunUMAP(qry, reduction = "scpred", dims = 1:30)
pdf("scPred.scpredPredictionUMAP.pdf", width=5, height=5)
DimPlot(qry, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
dev.off()

pred_numb <- crossTab(qry, "cell_type", "scpred_prediction")
pred_prop <- crossTab(qry, "cell_type", "scpred_prediction", output = "prop")
# write.table(pred_numb, file="pred_numb.txt", quote=F, sep="\t")
# write.table(pred_prop, file="pred_prop.txt", quote=F, sep="\t")
pdf("PredProp.PieChart.pdf", height=3, width=3*5)
par(mfrow=c(1,5))
pie(pred_prop[,1], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,1]*100,1),"%",sep=""))
pie(pred_prop[,2], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,2]*100,1),"%",sep=""))
pie(pred_prop[,3], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,3]*100,1),"%",sep=""))
pie(pred_prop[,4], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,4]*100,1),"%",sep=""))

plot(1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("center",c("d_EPI","d_PE","d_TE","d_UD","unassigned"),col=cccol[c(3,6,1,8,13)],pch=15,bty="n")
dev.off()

############ E7
ref_E7 <- getFeatureSpace(ref_E7, "cell_type")
ref_E7 <- trainModel(ref_E7, model = "glm")
get_scpred(ref_E7)
qry <- scPredict(qry, ref_E7, recompute_alignment=T, threshold=0.5)
pred_numb <- crossTab(qry, "cell_type", "scpred_prediction")
pred_prop <- crossTab(qry, "cell_type", "scpred_prediction", output = "prop")
# write.table(pred_numb, file="E7_pred_numb.txt", quote=F, sep="\t")
# write.table(pred_prop, file="E7_pred_prop.txt", quote=F, sep="\t")
pdf("PredProp.PieChart.E7.pdf", height=3, width=3*4)
par(mfrow=c(1,4))
pie(pred_prop[,1], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,1]*100,1),"%",sep=""))
pie(pred_prop[,2], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,2]*100,1),"%",sep=""))
pie(pred_prop[,3], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,3]*100,1),"%",sep=""))
pie(pred_prop[,4], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,4]*100,1),"%",sep=""))
dev.off()

############ E8
ref_E8 <- getFeatureSpace(ref_E8, "cell_type")
ref_E8 <- trainModel(ref_E8, model = "glm")
get_scpred(ref_E8)
qry <- scPredict(qry, ref_E8, recompute_alignment=T, threshold=0.5)
pred_numb <- crossTab(qry, "cell_type", "scpred_prediction")
pred_prop <- crossTab(qry, "cell_type", "scpred_prediction", output = "prop")
# write.table(pred_numb, file="E8_pred_numb.txt", quote=F, sep="\t")
# write.table(pred_prop, file="E8_pred_prop.txt", quote=F, sep="\t")
pdf("PredProp.PieChart.E8.pdf", height=3, width=3*4)
par(mfrow=c(1,4))
pie(pred_prop[,1], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,1]*100,1),"%",sep=""))
pie(pred_prop[,2], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,2]*100,1),"%",sep=""))
pie(pred_prop[,3], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,3]*100,1),"%",sep=""))
pie(pred_prop[,4], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,4]*100,1),"%",sep=""))
dev.off()

############ E9
ref_E9 <- getFeatureSpace(ref_E9, "cell_type")
ref_E9 <- trainModel(ref_E9, model = "glm")
get_scpred(ref_E9)
qry <- scPredict(qry, ref_E9, recompute_alignment=T, threshold=0.5)
pred_numb <- crossTab(qry, "cell_type", "scpred_prediction")
pred_prop <- crossTab(qry, "cell_type", "scpred_prediction", output = "prop")
# write.table(pred_numb, file="E9_pred_numb.txt", quote=F, sep="\t")
# write.table(pred_prop, file="E9_pred_prop.txt", quote=F, sep="\t")
pdf("PredProp.PieChart.E9.pdf", height=3, width=3*4)
par(mfrow=c(1,4))
pie(pred_prop[,1], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,1]*100,1),"%",sep=""))
pie(pred_prop[,2], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,2]*100,1),"%",sep=""))
pie(pred_prop[,3], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,3]*100,1),"%",sep=""))
pie(pred_prop[,4], col=cccol[c(3,6,1,8,13)], border="white", labels=paste(round(pred_prop[,4]*100,1),"%",sep=""))
dev.off()

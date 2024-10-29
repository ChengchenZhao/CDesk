
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
setwd("/mnt/lvyawen/Works/works1/ChickEmbryo/4.figure/5.test/cluster_PGC_diffused")
set.seed(1118)
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
cccol01 <- paste(cccol,"01",sep="")
cccol05 <- paste(cccol,"05",sep="")
cccol30 <- paste(cccol,"30",sep="")
cccol50 <- paste(cccol,"50",sep="")
cccol80 <- paste(cccol,"80",sep="")
allcol <- c("#f7acbc","#deab8a","#817936","#444693","#ef5b9c","#fedcbd","#7f7522","#2b4490","#f47920","#80752c","#2a5caa","#f05b72","#905a3d","#87843b","#224b8f","#f15b6c","#8f4b2e","#726930","#003a6c","#f8aba6","#87481f","#454926","#102b6a","#f69c9f","#5f3c23","#2e3a1f","#426ab3","#f58f98","#6b473c","#4d4f36","#46485f","#ca8687","#faa755","#b7ba6b","#4e72b8","#f391a9","#fab27b","#b2d235","#181d4b","#bd6758","#f58220","#5c7a29","#1a2933","#d71345","#843900","#bed742","#121a2a","#d64f44","#905d1d","#7fb80e","#0c212b","#d93a49","#8a5d19","#a3cf62","#6a6da9","#b3424a","#8c531b","#769149","#585eaa","#c76968","#826858","#6d8346","#494e8f","#bb505d","#64492b","#78a355","#afb4db","#987165","#ae6642","#abc88b","#9b95c9","#ac6767","#56452d","#74905d","#6950a1","#973c3f","#96582a","#cde6c7","#6f60aa","#b22c46","#705628","#1d953f","#867892","#a7324a","#4a3113","#77ac98","#918597","#aa363d","#412f1f","#007d65","#6f6d85","#ed1941","#845538","#84bf96","#594c6d","#f26522","#8e7437","#45b97c","#694d9f","#d2553d","#69541b","#225a1f","#6f599c","#b4534b","#d5c59f","#367459","#8552a1","#ef4136","#cd9a5b","#007947","#543044","#c63c26","#cd9a5b","#40835e","#63434f","#f3715c","#b36d41","#2b6447","#7d5886","#a7573b","#df9464","#005831","#401c44","#aa2116","#b76f40","#006c54","#472d56","#b64533","#ad8b3d","#375830","#45224a","#b54334","#dea32c","#274d3d","#411445","#853f04","#d1923f","#375830","#4b2f3d","#840228","#c88400","#27342b","#402e4c","#7a1723","#c37e00","#65c294","#c77eb5","#a03939","#c37e00","#73b9a2","#ea66a6","#8a2e3b","#e0861a","#72baa7","#f173ac","#8e453f","#ffce7b","#005344","#8f4b4a","#fcaf17","#122e29","#892f1b","#ba8448","#293047","#f6f5ec","#6b2c25","#896a45","#00ae9d","#d9d6c3","#733a31","#76624c","#508a88","#d1c7b7","#54211d","#6d5826","#70a19f","#f2eada","#78331e","#ffc20e","#50b7c1","#d3d7d4","#53261f","#fdb933","#00a6ac","#999d9c","#f15a22","#d3c6a6","#78cdd1","#a1a3a6","#b4533c","#c7a252","#008792","#9d9087","#84331f","#dec674","#94d6da","#8a8c8e","#f47a55","#b69968","#afdfe4","#74787c","#f15a22","#c1a173","#5e7c85","#7c8577","#f3704b","#dbce8f","#76becc","#72777b","#da765b","#ffd400","#90d7ec","#77787b","#c85d44","#ffd400","#009ad6","#4f5555","#ae5039","#ffe600","#145b7d","#6c4c49","#6a3427","#f0dc70","#11264f","#563624","#8f4b38","#fcf16e","#7bbfea","#3e4145","#8e3e1f","#decb00","#33a3dc","#3c3645","#f36c21","#cbc547","#228fbd","#464547","#b4532a","#6e6b41","#2468a2","#130c0e","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")

# 加载必要的库
library(amap)  # 用于KMeans函数
library(dplyr)  
library(ggplot2)  
library(cluster)  

# 设置科学计数法的显示选项
options(scipen = 200)

# Z-Score标准化函数
ZScore <- function(x) {
  x <- as.numeric(x)
  return((x - mean(x)) / sqrt(var(x)))
}

# 读取数据的函数
read_data <- function(fpkm_file) {
  # 读取主数据集
  fpkm <- read.table(fpkm_file, header = TRUE, row.names = 1, sep = "\t")
  all_genes <- rownames(fpkm)
  all_samples <- colnames(fpkm)
  logfpkm <- log2(fpkm + 1)
  
  return(list(fpkm = fpkm, logfpkm = logfpkm, all_genes = all_genes, all_samples = all_samples))
}

# Z-Score归一化函数
normalize_data <- function(logfpkm) {
  zscore <- apply(logfpkm, 2, ZScore)
  return(zscore)
}

# 运行KMeans聚类并导出结果的函数
run_kmeans <- function(zscore, genes_to_cluster, k = 12, seed = 999) {
  set.seed(seed)
  clusterData <- zscore[genes_to_cluster, ]
  km <- Kmeans(clusterData, k, method = "correlation", iter.max = 100, nstart = 10)
  
  # 保存聚类的基因
  for (id in seq(k)) {
    modGenes = names(which(km$cluster == id))
    fileName = paste("Zscore.cluster_", id, "_genes.txt", sep = "")
    write.table(modGenes, file = fileName, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  return(km)
}

# 绘制聚类结果的函数
plot_clusters <- function(km, plotData, all_samples, output_prefix = "Fig") {
  pdf(file = paste0(output_prefix, ".ExpClustersLinePlot.OnePath.ZScore.pdf"), width = 5.55, height = 3)
  selected_cluster <- seq(length(unique(km$cluster)))
  clusterCol <- cccol[1:length(unique(km$cluster))]
  
  for (each_i in seq(length(selected_cluster))) {
    each <- selected_cluster[each_i]
    modGenes <- names(which(km$cluster == each))
    v1 <- apply(plotData[modGenes, ], 2, mean)
    n <- length(modGenes)
    sd <- apply(plotData[modGenes, ], 2, sd)
    
    # 计算置信区间
    v2 <- v1 - sd/sqrt(n)*qt(1 - 0.05/2, n - 1)
    v3 <- v1 + sd/sqrt(n)*qt(1 - 0.05/2, n - 1)
    
    par(mar = c(6.4, 4, 2, 2))
    plot(v1, lwd = 3, type = "l", col = clusterCol[each_i], ylim = c(min(v2), max(v3)), xlab = NA, ylab = "Zscore", xaxt = "n", main = paste("C", each, "(n=", length(modGenes), ")", sep = ""))
    axis(side = 1, 1:length(all_samples), all_samples, las = 2)
    polygon(c(1, 1:length(all_samples), length(all_samples):2), c(v2[1], v3, v2[length(all_samples):2]), col = adjustcolor("grey", alpha.f = 0.4), border = NA)
  }
  dev.off()
}

# Step 1: Load Data
result <- read_data("fpkm.txt")  # 替换为你的数据集路径

# Step 2: Data Preprocessing
data_filtered <- result$logfpkm %>%
  filter_if(is.numeric, all_vars(. > 1))  # 示例过滤，移除低表达

# Step 3: Identify Highly Variable Genes
highly_var_genes <- rowVars(as.matrix(data_filtered))
top_genes <- order(highly_var_genes, decreasing = TRUE)[1:2000]  # 选择前200个高变异基因
data_top <- data_filtered[top_genes, ]

# Step 4: Dimensionality Reduction with PCA
pca <- prcomp(t(data_top), center = TRUE, scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], col = "blue", xlab = "PC1", ylab = "PC2", main = "PCA of Top Genes")

# Step 5: Clustering with K-means
set.seed(123)  # 设置随机种子以保证可重复性
km <- run_kmeans(normalize_data(data_top), top_genes, k = 3)  # 根据数据调整聚类数
data_filtered$Cluster <- factor(km$cluster)

# Step 6: Visualize Clusters
ggplot(data_filtered, aes(x = pca$x[,1], y = pca$x[,2], color = Cluster)) +
  geom_point() +
  labs(title = "Gene Clusters in PCA Space", x = "PC1", y = "PC2") +
  theme_minimal()

# Step 7: Save Clustering Results
write.csv(data_filtered, "clustered_data.csv", row.names = FALSE)

# Step 8: 绘制聚类结果
plot_clusters(km, data_top, result$all_samples)


  File "/mnt/runtask/CDesk/src/Tasks/BulkRNAseqBW/__init__.py", line 75, in run
    self.exec_cmd(" ".join(cmd), logger)
  File "CDesk_base.py", line 428, in CDesk_base.taskTemplate.exec_cmd
KeyboardInterruptdd


# 示例用法
# result <- read_data("fpkm.txt", "fpkm_dev.txt")
# merged_data <- merge_data(result$fpkm, result$fpkm_dev, c("mESC", "tCLC"), c("C2"), result$common_genes)
# zscore_data <- normalize_data(merged_data)
# km_result <- run_kmeans(zscore_data, genes_to_cluster = common_genes, path_samples = c("mESC", "C2"))
# plot_clusters(km_result, zscore_data, c("mESC", "tCLC", "C2"))





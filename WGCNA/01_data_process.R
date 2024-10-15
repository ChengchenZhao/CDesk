### 测序数据
rm(list = ls())
### 加载R包
library(WGCNA)
### 读取表达量数据
raw_counts <- data.table::fread('count_table.csv', data.table = F)
### 第一列变行名
rownames(raw_counts) <- raw_counts[,1]
raw_counts <- raw_counts[,-1]
### 查看数据
head(raw_counts)
dim(raw_counts)
test <- raw_counts[1:10,1:10]

###进行行列转换 t(),
datExpr0 <- as.data.frame(t(raw_counts))

###################################################################################
### 1.goodSamplesGenes 的输入条件
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

### 如果没有达标就需要筛选

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    #printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

### 力度不够,还可以再做一些筛选工作
raw_counts <- as.data.frame(t(datExpr0))

### 读入性状数据
library(RColorBrewer)
pheno <- data.table::fread('sample_metadata.csv',data.table = F)
num_conditions <- nlevels(as.factor(pheno$condition))
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
cond_colors <- pal[as.integer(as.factor(pheno$condition))]

library(gplots)
heatmap.2(cor(raw_counts), 
          RowSideColors=cond_colors,
          trace='none', 
          main='Sample correlations (raw)')

### 2.低counts数据过滤
low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)
raw_counts <- raw_counts[!low_count_mask,]
sum(low_count_mask)

### 3。log2转换
log_counts <- log2(raw_counts + 1)

### 数据变形
x = reshape::melt(as.matrix(log_counts))
colnames(x) = c('gene_id', 'sample', 'value')
library(ggplot2)
ggplot(x, aes(x=value, color=sample)) + 
  geom_density()

### 没有log化的
x1 = reshape::melt(as.matrix(raw_counts))
colnames(x1) = c('gene_id', 'sample', 'value')
ggplot(x1, aes(x=value, color=sample)) + 
  geom_density()

heatmap.2(cor(log_counts), 
          RowSideColors=cond_colors,
          trace='none', main='Sample correlations (log2-transformed)')

### 4.去除没有改变的基因
# first, let's remove any genes with _zero_ variance since these are not
# going to help us, and may cause problems with some of the models
log_counts <- log_counts[apply(log_counts, 1, var) > 0,]

### 聚类看异常样本
sampleTree = hclust(dist(t(log_counts)), method = "average")
### 画图
pdf("sample_clustering_to_detect_outliers.pdf", width = 6, height = 6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

### 提取表达量和性状信息的交集
rownames(pheno) <- pheno$sample_id
index <- intersect(pheno$sample_id,colnames(log_counts))
pheno <- pheno[index,]
datExpr = as.data.frame(t(log_counts[,index]))

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(as.numeric(as.factor(pheno$condition)), signed = FALSE)

pdf("sample_dendrogram_and_trait_heatmap.pdf", width = 8, height = 6)
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(pheno)[2], 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

save(datExpr,pheno,file = "datExpr0.Rdata")

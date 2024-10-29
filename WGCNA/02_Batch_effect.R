################################################
### batch effect
load(file = "datExpr0.Rdata")
sampleTree = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(as.numeric(as.factor(pheno$condition)), signed = FALSE)
plotDendroAndColors(sampleTree, 
                    traitColors,
                    groupLabels = names(pheno)[2], 
                    main = "Sample dendrogram and trait heatmap")

### 批次矫正

model <- model.matrix(~as.factor(pheno$condition))
library(sva)
datExpr <- ComBat(dat = as.matrix(t(datExpr)), batch = pheno$batch, mod = model)

### 作图查看
sampleTree = hclust(dist(t(datExpr)), method = "average")
traitColors = numbers2colors(as.numeric(as.factor(pheno$condition)), signed = FALSE)
pdf("sample_dendrogram_and_trait_heatmap_removeBatchEffect.pdf", width = 8, height = 6)
plotDendroAndColors(sampleTree, 
                    traitColors,
                    groupLabels = names(pheno)[2], 
                    main = "Sample dendrogram and trait heatmap-removeBatchEffect")
dev.off()

### quntil normalization
### 理解 Quntile Normalization 
###  https://mp.weixin.qq.com/s/YL3UeWzA6Ysi5b3nnTxiww
color <- as.numeric(as.factor(pheno$condition))
boxplot(datExpr,outline=FALSE, notch=T, las=2,col=color)
## 组间校正
library(limma) 
datExpr=normalizeBetweenArrays(datExpr)
pdf("boxplot.pdf", width = 8, height = 6)
boxplot(datExpr,outline=FALSE, notch=T, las=2,col=color)
dev.off()

### 热图
num_conditions <- nlevels(as.factor(pheno$condition))
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
cond_colors <- pal[as.integer(as.factor(pheno$condition))]

pdf("Sample correlations (batch effect).pdf", width = 8, height = 6)
heatmap.2(cor(datExpr), 
          RowSideColors=cond_colors,
          trace='none', 
          main='Sample correlations (batch effect)')
dev.off()

### 行列转置
datExpr <- t(datExpr)
save(datExpr,pheno,file = "remove_batch.Rdata")

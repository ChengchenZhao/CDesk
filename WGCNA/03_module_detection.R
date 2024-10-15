#####WGCNA
library(WGCNA)

### 设置多线程
enableWGCNAThreads()
### 加载上一步的数据
load(file = "remove_batch.Rdata")

### 软阈值的确定，使用函数pickSoftThreshold

powers = c(c(1:10), seq(from = 12, to=30, by=2))
length(powers)
### 这个power的长度，决定了最终返回结果的行数目

# Call the network topology analysis function
### 模块这里都是作用于表达量数据，没有性状数据的事情
### 稍微耗时间。
#pickSoftThreshold用于计算不同软阈值功率下的网络特性，帮助选择最佳的功率值。
sft = pickSoftThreshold(datExpr, powerVector = powers,networkType="signed", verbose = 5)

test <- sft$fitIndices# 包含了每个功率值对应的网络特征信息，例如无标度拟合指数和平均连接度。

### 画图呈现结果
pdf("Scale_independence_and_Mean_connectivity.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9#将字体大小设为 0.9 倍默认大小
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

# 一步法网络构建以及模块发现
# power = 12，这里设置了软阈值功率为 12。这个值通常通过之前的分析（如 pickSoftThreshold()）来确定

# TOMType 指定使用哪种类型的 TOM（拓扑重叠矩阵），这里选择 "unsigned" 表示不区分相关性的符号，正负相关都会被同等对待。若设置为 "signed"，则会区分正负相关。
# networkType 网络类型设置为 "signed"，表示将区分正相关和负相关的基因关系，负相关的基因在构建网络时将被弱化。
# reassignThreshold 设定模块重分配的阈值。如果基因与其模块的相关性小于这个阈值，则会被重新分配。设置为 0 表示不重分配基因。
# 设定模块合并的高度阈值。模块之间的聚类树距离小于 0.25 的模块将被合并，这个值越低，合并越严格。
# numericLabels 将模块标签设置为数值而非颜色。这样每个模块的标签是数字，而不是默认的颜色名称（如 "blue", "turquoise" 等）。
# pamRespectsDendro 参数的选择取决于您对模块划分的需求：如果希望保持树状结构一致性，使用 TRUE；如果希望最大化模块的内部一致性，可以使用 FALSE，允许 PAM 进行更自由的分配。
# saveTOMs 表示是否保存计算得到的 TOM 矩阵。TOM 是用于衡量基因对之间相似度的矩阵，它是网络构建的重要中间结果。
# saveTOMFileBase = "yourname"：设置保存的 TOM 文件的基础文件名。
# verbose 设定输出详细程度，值越高输出越详细

cor <- WGCNA::cor#重载了 cor() 函数
net = blockwiseModules(datExpr, power = 12,
                       TOMType = "unsigned", minModuleSize = 30,
                       networkType="signed",
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "yourname", 
                       verbose = 3)
cor<-stats::cor#stats 包中的 cor() 函数重新分配给当前环境中的变量 cor，以便在后续代码中使用默认的 cor() 函数。

### 返回的net是个列表，Module Epigene 表达数据也在里面
table(net$colors)
### 0表示的是不在模块内的基因

### 画图
pdf("Dendrogram_with_Module_Colors.pdf", width = 12, height = 9)
### 标签转为颜色
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plot(net$dendrograms[[1]])
plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

### 保存数据
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
geneTree = net$dendrograms[[1]]
save(moduleLabels, moduleColors, geneTree, file = "module_detection.RData")




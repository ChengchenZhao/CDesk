###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(monocle)
set.seed(0)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################

##### expression matrix
# raw_data <- read.table("/mnt/Storage/home/zhaochengchen/Work/2.naiveES/result/scRNAseq/0.mapping/20171227/naiveES_20171227_fpkm.txt",row.names=1,header=T)
selected_cells <- as.vector(read.table("selected_cells.txt")[,1])
raw_data <- read.table("naiveES_20171227_fpkm.txt",row.names=1,header=T)
data <- raw_data[,selected_cells]
logfpkm <- log2(data+1)
d0 <- grep("d0_", selected_cells)
d1 <- grep("d1_", selected_cells)
d3 <- grep("d3_", selected_cells)
d6 <- grep("d6_", selected_cells)
d9 <- grep("d9_", selected_cells)
d12 <- grep("d12_", selected_cells)
F10 <- grep("F10_", selected_cells)

sample_stage <- rep("",length(selected_cells))
sample_stage[d0] <- "d0"
sample_stage[d1] <- "d1"
sample_stage[d3] <- "d3"
sample_stage[d6] <- "d6"
sample_stage[d9] <- "d9"
sample_stage[d12] <- "d12"
sample_stage[F10] <- "F10"

expr_matrix <- data
sample_sheet <- data.frame(stage=as.character(sample_stage));rownames(sample_sheet) <- selected_cells
gene_annotation <- data.frame(gene_short_name=rownames(logfpkm));rownames(gene_annotation) <- rownames(logfpkm)

# negbinomial.size()	UMIs, Transcript counts from experiments with spike-ins or relative2abs(), raw read counts	Negative binomial distribution with fixed variance (which is automatically calculated by Monocle). Recommended for most users.
# negbinomial()	UMIs, Transcript counts from experiments with spike-ins or relative2abs, raw read counts	Slightly more accurate than negbinomial.size(), but much, much slower. Not recommended except for very small datasets.
# tobit()	FPKM, TPM	Tobits are truncated normal distributions. Using tobit() will tell Monocle to log-transform your data where appropriate. Do not transform it yourself.
# gaussianff()	log-transformed FPKM/TPMs, Ct values from single-cell qPCR	If you want to use Monocle on data you have already transformed to be normally distributed, you can use this function, though some Monocle features may not work well.
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
mono_data <- newCellDataSet(as.matrix(expr_matrix),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit=0.5,
                       expressionFamily=tobit())
mono_data$stage <- factor(mono_data$stage,levels=c("d0","d1","d3","d6","d9","d12","F10"))
###################################################################################################
################################                PLOT               ################################
###################################################################################################

# estimateSizeFactors() and estimateDispersions() will only work, and are only needed, if you are working with a CellDataSet with a negbinomial() or negbinomial.size() expression family.
# mono_data <- estimateSizeFactors(mono_data)
# mono_data <- estimateDispersions(mono_data)

mono_data <- detectGenes(mono_data, min_expr = 0.01)
expressed_genes <- row.names(subset(fData(mono_data), num_cells_expressed >= ncol(expr_matrix)*0.1))
mono_data <- mono_data[expressed_genes,]

###### differential expressed gene
diff_test_res <- differentialGeneTest(mono_data[expressed_genes,],
              fullModelFormulaStr = "~stage")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
mono_data <- setOrderingFilter(mono_data, ordering_genes)
mono_data <- reduceDimension(mono_data, max_components=2)
mono_data <- orderCells(mono_data)
mono_data <- orderCells(mono_data,root_state=3)

pdf('PseudotimeAnalysis.pdf',width=4,height=4)
plot_cell_trajectory(mono_data,color_by="stage",show_branch_points=F)
plot_cell_trajectory(mono_data,color_by="Pseudotime",show_branch_points=F)
dev.off()

pdf('PseudotimeAnalysis_EachBranch.pdf',width=6,height=6)
plot_cell_trajectory(mono_data,color_by="stage",show_branch_points=F) + facet_wrap(~stage, nrow = 3)
dev.off()

# gene expression in Pseudotime
selected_genes <- row.names(subset(fData(mono_data),gene_short_name %in% c("CBX5","DNMT3L")))
# selected_genes <- row.names(subset(fData(mono_data),gene_short_name %in% c("ZNF578","ZNF600")))
cds_subset <- mono_data[selected_genes,]

pdf('MarkerGeneExpressionInPseudotimeOrder.pdf',width=6,height=6)
plot_genes_in_pseudotime(cds_subset, color_by = "stage")
plot_genes_branched_pseudotime(cds_subset,branch_point = 1,color_by = "stage",ncol = 1,branch_labels = c("d0_d1","d6_F10"))
plot_genes_jitter(cds_subset, color_by = "stage")
plot_genes_jitter(cds_subset,grouping = "stage",color_by = "stage",nrow= NULL,ncol = 1)
dev.off()



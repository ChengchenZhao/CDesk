#! /usr/bin/Rscript



# Description : Perform the GO/KEGG enrichment analysis based on the pubsihed GO/KEGG gene sets
# Author      : Zheting Zhang
# Time        : 2024-11-17



# 检测并安装必要的包

packages <- c("BiocManager", "clusterProfiler", "org.Hs.eg.db", "org.Rn.eg.db", "org.Mm.eg.db", "org.Gg.eg.db", "org.Ss.eg.db", "enrichplot", "ggplot2","clusterProfiler")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "BiocManager") {
      install.packages("BiocManager")
    } else {
      BiocManager::install(pkg)
    }
  }
}

#Load packages

suppressMessages(library(dplyr))
library(org.Hs.eg.db)
library(org.Rn.eg.db)
library(org.Mm.eg.db) %>% suppressMessages() 
library(org.Gg.eg.db)
library(org.Ss.eg.db)
library(enrichplot) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()
# 定义物种数据库映射

species_dbs <- list(
  "human" = org.Hs.eg.db,
  "rat" = org.Rn.eg.db,
  "mouse" = org.Mm.eg.db,
  "chicken" = org.Gg.eg.db,
  "pig" = org.Ss.eg.db)


# 定义函数进行GO或KEGG富集分析

perform_enrichment <- function(gene_file, analysis_type, species, organism) {
  # 读取基因列表文件

  gene_list <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)
  # 选择物种数据库

  org_db <- species_dbs[[species]]


  # transform geneID

  gene.df <- bitr(geneID = gene_list$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)

  # 进行GO富集分析

  if (analysis_type == "GO") {
    enrich_result <- enrichGO(gene = gene.df$ENTREZID,
                              OrgDb = org_db,
                              ont = "ALL",
        		                  pAdjustMethod = "BH",
			                       pvalueCutoff = 0.05,
			                       qvalueCutoff = 0.5,
			                       keyType = "ENTREZID",
                              readable = TRUE)
  }
  
  # 进行KEGG富集分析

  if (analysis_type == "KEGG") {
    enrich_result <- enrichKEGG(gene = gene.df$ENTREZID,
                                organism = organism,
			                         	keyType = "kegg",
			                        	pAdjustMethod = "BH",
			                         	minGSSize = 1, 
			                         	maxGSSize = 5000,
			                         	pvalueCutoff = 0.05,
		                       		qvalueCutoff = 0.5,
                                readable = TRUE,
			                       	keyType = "ENTREZID")
  }
  
  # 输出结果表格

#  output_file <- paste0(gene_file,analysis_type, "_enrichment_results.txt")
#  write.table(as.data.frame(enrich_result), file = output_file, row.names = FALSE, sep = "\t")
  
  # 返回结果对象

  return(enrich_result)
}

# 去冗余函数
simplify_enrichment <- function(enrich_result, cutoff) {
  if (nrow(enrich_result) > 0) {
    simplified_result <- simplify(enrich_result, cutoff = cutoff)
  } else {
    simplified_result <- enrich_result
    
  }
  return(simplified_result)
}

# 从命令行参数读取输入

args <- commandArgs(trailingOnly = TRUE)

if ( length(args) < 3) {
  stop("Usage: Rscript enrichment_analysis.R <gene_file(SYMBOL)> GO <species> <similarity cutoff level[0:1]>
           or: Rscript enrichment_analysis.R <gene_file(SYMBOL)> KEGG <species>  <organism>")
}

if (args[2] == "GO" & length(args) != 4) {
  stop("Usage: Rscript enrichment_analysis.R <gene_file(SYMBOL)> GO <species> <similarity cutoff level[0:1]> ")
}

if (args[2] == "KEGG" & length(args) != 4) {
  stop("Usage: Rscript enrichment_analysis.R <gene_file(SYMBOL)> KEGG <species>  <organism> ")
}



gene_file <- args[1]
analysis_type <- args[2]
species <- args[3]
if (analysis_type == "GO") {
  cutoff <- args[4] %>% as.numeric()
}

if (analysis_type == "KEGG") {
  organism <- args[4]
}




# 检查输入是否有效

if (!file.exists(gene_file)) {
  stop("Gene list file does not exist")
}
if (!(analysis_type %in% c("GO", "KEGG"))) {
  stop("Analysis type must be GO or KEGG")
}
if (!(species %in% names(species_dbs))) {
  stop("Species must be human, rat, mouse, chicken, or pig")
}
if (analysis_type == "KEGG" & length(args)!= 6) {
  stop("Please provide organism for KEGG")
}
if (cutoff <= 0 || cutoff > 1){
  stop("cutoff value must be [0:1]")
}


# 执行富集分析

enrichment_result <- perform_enrichment(gene_file, analysis_type, species,organism)

# 合并冗余项

if (analysis_type == "GO") {
  enrichment_result <- simplify_enrichment(enrichment_result,cutoff)
}

# 层级筛选

#enrichment_result <- gofilter(enrichment_result,level = level_low:level_high)

# 打印结果
output_file <- paste0(gene_file,"_",analysis_type,cutoff,"cutoff","_enrichment_results.txt")
write.table(as.data.frame(enrichment_result), file = output_file, row.names = FALSE, sep = "\t",quote=FALSE)

cat("Thanks for using ! ^_^ , any question please contact with Bioinformatic Team of CFCL\n")

# 检测并安装必要的包

packages <- c( "ggplot2","clusterProfiler")

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
library(ggplot2) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()
# 定义物种数据库映射



  # 定义函数进行富集分析

  perform_enrichment <- function(gene_file,customer_term) {
  
  # 读取基因列表文件

  gene_list <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)
  
  # 定义自定义term

  customer_term <- read.table(customer_term,header = FALSE,stringsAsFactors = FALSE,col.names = c("Name","Gene"),sep = "\t")
  
  if (ncol(customer_term)<2) {
    stop("customer file needs two colums: functional term and gene(symbol)")
  }
  
  customer_term_id <- unique(customer_term$Name)
  
  customer_term_id <- data.frame(Term=paste0("term_",c(1:length(customer_term_id))),Name=customer_term_id)
  
  customer_term <- merge(customer_term,customer_term_id,by="Name") %>% select(Term,Gene)
  

  # 进行富集分析

    enrich_result <- enricher(gene = gene_list$V1,
        		                  pAdjustMethod = "BH",
			                        pvalueCutoff = 0.05,
			                        qvalueCutoff = 0.5,
			                        TERM2GENE = customer_term,
			                        TERM2NAME = customer_term_id)
  
  

  
  # 输出结果表格

#  output_file <- paste0(gene_file,analysis_type, "_enrichment_results.txt")
#  write.table(as.data.frame(enrich_result), file = output_file, row.names = FALSE, sep = "\t")
  
  # 返回结果对象

  return(enrich_result)
}

# 从命令行参数读取输入

args <- commandArgs(trailingOnly = TRUE)

if ( length(args) < 2) {
  stop("Usage: Rscript enrichment_analysis_customer.R <gene_file(SYMBOL)> <customer_term(term\tsymbol)>")
}


gene_file <- args[1]
customer_term <- args[2]



# 检查输入是否有效

if (!file.exists(gene_file)) {
  stop("Gene list file does not exist")
}
if (!file.exists(customer_term)) {
  stop("customer_term file does not exist")
}


# 执行富集分析

enrichment_result <- perform_enrichment(gene_file, customer_term)

# 打印结果
output_file <- paste0(gene_file,"_",customer_term,"_enrichment_results.txt")
write.table(as.data.frame(enrichment_result), file = output_file, row.names = FALSE, sep = "\t",quote=FALSE)

cat("Thanks for using ! ^_^ , any question please contact with Bioinformatic Team of CFCL\n")

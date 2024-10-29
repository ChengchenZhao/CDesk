# 加载必要的包

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}



library(ggplot2)

# 定义函数读取和处理GO/KEGG结果表格

#enrich_result <- read.table(enrichment_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# 定义函数绘制气泡图


plot_bubble <- function(enrich_result,high_color,low_color) {
  
  if (colnames(enrich_result)[1]=="ONTOLOGY") {
  
  p <- ggplot(enrich_result)+
    geom_point(aes(x=GeneRatio,y=Description,col=-log10(p.adjust),size=Count))+
    theme_bw()+
    scale_color_gradient(high = high_color,low = low_color)+
    facet_grid(ONTOLOGY~.,scales = "free_y",space = "free_y")+
    theme(text = element_text(size = 25))
  return(p)
}
  if (colnames(enrich_result)[1]!="ONTOLOGY") {
    p <- ggplot(enrich_result)+
      geom_point(aes(x=GeneRatio,y=Description,col=-log10(p.adjust),size=Count))+
      theme_bw()+
      scale_color_gradient(high = high_color,low = low_color)+
      theme(text = element_text(size = 25)) 
    return(p)
  }}
# 定义函数绘制柱状图

    
  
plot_bar <- function(enrich_result) {
  if (colnames(enrich_result)[1]=="ONTOLOGY") {
  p <- ggplot(enrich_result) +
    geom_bar(aes(x=-log10(p.adjust),y=Description,fill=ONTOLOGY),stat = "identity")+
    theme_bw() +
    facet_grid(ONTOLOGY~.,scales = "free_y",space = "free_y")+
    theme(text = element_text(size = 25))
  return(p)
}
  if (colnames(enrich_result)[1]!="ONTOLOGY") {
    p <- ggplot(enrich_result) +
      geom_bar(aes(x=-log10(p.adjust),y=Description),stat = "identity")+
      theme_bw() +
      theme(text = element_text(size = 25))
    return(p)  
  }}

# 从命令行参数读取输入

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot_enrichment.R <enrichment_file> <color for high value> <color for low value>")
}

enrichment_file <- args[1]
high_color <- args[2]
low_color <- args[3]


# 读取GO/KEGG结果表格

enrich_result <- read.table(enrichment_file,header = TRUE,sep = "\t")
enrich_result$Description <- factor(enrich_result$Description,levels = as.character(enrich_result$Description[rev(order(enrich_result$p.adjust))]))
enrich_result$GeneRatio <-  sapply(enrich_result$GeneRatio, cal <- function(x){out=eval(parse(text = x)); return(out)})


# 检查输入是否有效

if (!file.exists(enrichment_file)) {
  stop("Enrichment file does not exist")
}



# 绘制气泡图
pdf(paste0(enrichment_file,"_bubble.pdf"),width = 20,height = 10)
plot_bubble(enrich_result,high_color,low_color)
dev.off()

# 绘制柱状图

pdf(paste0(enrichment_file,"_bar.pdf"),width = 20,height = 10)
plot_bar(enrich_result)
dev.off()

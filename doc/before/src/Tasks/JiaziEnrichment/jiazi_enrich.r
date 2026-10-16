#!/usr/bin/Rscript-4.2

# 关闭警告
options(warn=-1)

# 解析参数的包
library(optparse)

# 生成命令行参数
option_list = list(
    make_option(
        c("-t", "--term_prefix"),
        default="/mnt/liudong/data/Genome/jiazi_enrich/term",
        help="enrich the gene list with this prefix's term\n                which contains the [term_prefix]2gene.txt and the [term_prefix]2name.txt,\n                default is the mouse's term\n                [term_prefix]2name.txt: contains two columns(term, name)\n                [term_prefix]2gene.txt: contains two columns(term, gene)"
    ),
    make_option(
        c("-o", "--output"),
        default="result",
        help="output name with path,\n                default is the gene_list file name"
    ),
    make_option(
        c("-p", "--pvalueCutoff"),
        default=0.2,
        type="numeric",
        help="adjusted pvalue cutoff on enrichment tests to report"
    ),
    make_option(
        c("-q", "--qvalueCutoff"),
        default=0.05,
        type="numeric",
        help="qvalue cutoff on enrichment tests to report as significant"
    ),
    make_option(
        c("-m", "--pAdjustMethod"),
        default="BH",
        help="one of \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\", \"none\""
    ),
    make_option(
        c("--minGSSize"),
        default=10,
        type="numeric",
        help="minimal size of genes annotated for testing"
    ),
    make_option(
        c("--maxGSSize"),
        default=500,
        type="numeric",
        help="maximal size of genes annotated for testing"
    )
)
parser = OptionParser(
    usage="%prog [options] gene_list",
    option_list=option_list
)

arguments = parse_args(
    parser,
    positional_arguments=c(1)
)

file = arguments$args
opt = arguments$options
term2gene = paste(opt$term_prefix, "2gene.txt", sep="")
term2name = paste(opt$term_prefix, "2name.txt", sep="")
output = opt$output
enrich_pvalueCutoff = opt$pvalueCutoff
enrich_qvalueCutoff = opt$qvalueCutoff
enrich_pAdjustMethod = opt$pAdjustMethod
enrich_minGSSize = opt$minGSSize
enrich_maxGSSize = opt$maxGSSize

if (output == "result"){
    t = strsplit(file, "/")[[1]]
    output = gsub(".txt$", "", t[length(t)])
}

# 验证参数
# sprintf有长度限制，因此分割
cat(sprintf(
    "Check the args:\n      gene_list:  %s\n      term2gene:  %s\n      term2name:  %s\n",
    file,
    term2gene,
    term2name
))
cat(sprintf(
    "         output:  %s\n   pvalueCutoff:  %s\n   qvalueCutoff:  %s\n",
    output,
    enrich_pvalueCutoff,
    enrich_qvalueCutoff
))
cat(sprintf(
    "  pAdjustMethod:  %s\n      minGSSize:  %s\n      maxGSSize:  %s\n",
    enrich_pAdjustMethod,
    enrich_minGSSize,
    enrich_maxGSSize
))

err_list = c()

if(file.access(file) == -1) {
    err_list = c(err_list, sprintf("\nThe gene list file ( %s ) does not exist!", file))
}
if(file.access(term2gene) == -1) {
    err_list = c(err_list, sprintf("\nThe term2gene file ( %s ) does not exist!", term2gene))
}
if(file.access(term2name) == -1) {
    err_list = c(err_list, sprintf("\nThe term2name file ( %s ) does not exist!", term2name))
}
if(length(err_list) > 0){
    stop(err_list)
}

cat("Enriching the gene_list with the term ...\n")
out <- file("JiaziEnrich.log", open="a")
sink(out, type="message")

library(clusterProfiler)

gene_list = readLines(file)
term2gene = read.table(term2gene, sep = "\t", header = T)
term2name = read.table(term2name,sep = "\t", header = T)

eg <- enricher(
    gene_list,
    pAdjustMethod = enrich_pAdjustMethod,
    pvalueCutoff = enrich_pvalueCutoff,
    qvalueCutoff = enrich_qvalueCutoff,
    TERM2GENE = term2gene,
    TERM2NAME = term2name, 
    minGSSize = enrich_minGSSize,
    maxGSSize = enrich_maxGSSize
)
out_gene <- c()
out_name <- c()
for (each_termindex in seq(length(eg@result$ID))){
    out_name <- c(out_name,paste(eg@result$Description[each_termindex],sep="\t"))
    out_gene <- c(out_gene, paste(strsplit(eg@result$geneID,"/")[[each_termindex]],collapse=","))
}
names(out_gene) <- out_name
# eg@result$ID           eg@result$Description  eg@result$GeneRatio    eg@result$BgRatio      eg@result$pvalue       eg@result$p.adjust     eg@result$qvalue       eg@result$geneID       eg@result$Count
out_matrix <- cbind(eg@result$pvalue,eg@result$p.adjust,eg@result$qvalue,out_gene)
rownames(out_matrix) <- out_name
colnames(out_matrix) <- c("pvalue","p.adjust","qvalue","genes")
write.table(out_matrix,file=paste0(output,".JiaziEnrich.txt"),col.names=T,row.names=T,sep="\t",quote=F)

pdf(paste0(output,".JiaziEnrich.pdf"),width=10,height=4)
dotplot(eg,showCategory = 12,font.size = 10)
dev.off()

sink(type="message")
close(out)

cat(sprintf(
    "Enrich result: %s, %s\nFinish!\n",
    paste0(output,".JiaziEnrich.txt"),
    paste0(output,".JiaziEnrich.pdf")
))

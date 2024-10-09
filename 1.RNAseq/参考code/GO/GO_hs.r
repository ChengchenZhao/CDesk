###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)
options(download.file.method = "wget")
# GO over-representation test

enrichGOtest <- function( genes, inType, species, outfile ){
	suppressMessages(library(clusterProfiler))
	# change the input gene type (convert to ENTREZID type which fastly return result)
	gene.df <- bitr(geneID = genes, fromType = inType, toType = "ENTREZID", OrgDb = paste0("org.",species,".eg.db") )
	eg <- enrichGO( gene = gene.df$ENTREZID,
					OrgDb = paste0("org.",species,".eg.db"),
					ont = "BP",
					pAdjustMethod = "BH",
					pvalueCutoff = 0.05,
					qvalueCutoff = 0.2
		)
	# return(eg)
	out_gene <- c()
	out_name <- c()
	for (each_termindex in seq(length(eg@result$ID))){
		out_name <- c(out_name,paste(eg@result$ID[each_termindex],eg@result$Description[each_termindex],sep="\t"))
		tmp_gene <- bitr(geneID = strsplit(eg@result$geneID,"/")[[each_termindex]], fromType = "ENTREZID", toType = inType, OrgDb = paste0("org.",species,".eg.db"))[,2]
		out_gene <- c(out_gene, paste(tmp_gene,collapse=","))
	}
	names(out_gene) <- out_name
	# eg@result$ID           eg@result$Description  eg@result$GeneRatio    eg@result$BgRatio      eg@result$pvalue       eg@result$p.adjust     eg@result$qvalue       eg@result$geneID       eg@result$Count
	out_matrix <- cbind(eg@result$pvalue,eg@result$p.adjust,eg@result$qvalue,out_gene)
	rownames(out_matrix) <- out_name
	colnames(out_matrix) <- c("pvalue","p.adjust","qvalue","genes")
	write.table(out_matrix,file=paste0(outfile,".GOTerms.txt"),col.names=T,row.names=T,sep="\t",quote=F)

	pdf(paste0(outfile,".GO.pdf"),width=10,height=4)
	dotplot(eg,showCategory = 12,font.size = 10)
}

enrichKEGGtest <- function( genes, inType, species, organism, outfile ){
	suppressMessages(library(clusterProfiler))
	# change the input gene type (convert to ENTREZID type which fastly return result)
	gene.df <- bitr(geneID = genes, fromType = inType, toType = "ENTREZID", OrgDb = paste0("org.",species,".eg.db") )
	eg <- enrichKEGG( gene = gene.df$ENTREZID,
					organism = organism,
					keyType = "kegg",
					pAdjustMethod = "BH",
					minGSSize = 10, 
					maxGSSize = 500,
					pvalueCutoff = 0.05,
					qvalueCutoff = 0.2
		)
	# return(eg)
	out_gene <- c()
	out_name <- c()
	for (each_termindex in seq(length(eg@result$ID))){
		out_name <- c(out_name,paste(eg@result$ID[each_termindex],eg@result$Description[each_termindex],sep="\t"))
		tmp_gene <- bitr(geneID = strsplit(eg@result$geneID,"/")[[each_termindex]], fromType = "ENTREZID", toType = inType, OrgDb = paste0("org.",species,".eg.db"))[,2]
		out_gene <- c(out_gene, paste(tmp_gene,collapse=","))
	}
	names(out_gene) <- out_name
	# eg@result$ID           eg@result$Description  eg@result$GeneRatio    eg@result$BgRatio      eg@result$pvalue       eg@result$p.adjust     eg@result$qvalue       eg@result$geneID       eg@result$Count
	out_matrix <- cbind(eg@result$pvalue,eg@result$p.adjust,eg@result$qvalue,out_gene)
	rownames(out_matrix) <- out_name
	colnames(out_matrix) <- c("pvalue","p.adjust","qvalue","genes")
	write.table(out_matrix,file=paste0(outfile,".KEGGTerms.txt"),col.names=T,row.names=T,sep="\t",quote=F)

	pdf(paste0(outfile,".KEGG.pdf"),width=10,height=4)
	dotplot(eg,showCategory = 12,font.size = 10)
}

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
gene <- as.vector(read.table(sys_argv[1], header = FALSE)[,1])
out_name <- sys_argv[2]

###################################################################################################
################################                PLOT               ################################
###################################################################################################
enrichGOtest(gene, "SYMBOL", "Hs", out_name);dev.off()

# https://www.genome.jp/kegg/catalog/org_list.html
enrichKEGGtest(gene, "SYMBOL", "Hs", "hsa" , out_name);dev.off()


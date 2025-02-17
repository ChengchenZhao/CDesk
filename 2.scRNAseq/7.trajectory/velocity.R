library(Seurat)

sys_argv <- commandArgs(T)
input_seurat = sys_argv[1]
output_directory =sys_argv[2]

if (length(sys_argv) != 2) {
  stop("Input and output arguments should be provided.")
}

data = readRDS(input_seurat)

## cell ID
write.csv(row.names(data@meta.data),file=paste0(output_directory,"/cellID_obs.csv"),row.names=FALSE)
## UMAP
write.csv(Embeddings(data,reduction="umap"),file=paste0(output_directory,"/cell_embeddings.csv"))
##seurat clustering
write.csv(Idents(data),paste0(output_directory,"/clusters.csv"))

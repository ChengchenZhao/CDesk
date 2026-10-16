###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
sys_argv <- commandArgs(T)
input_file <- sys_argv[1]
output_directory <- sys_argv[2]
group_file <- sys_argv[3]
batch <- sys_argv[4]
plot_width <- as.numeric(sys_argv[5])
plot_height <- as.numeric(sys_argv[6])

suppressMessages(library(dplyr))

library("gplots") %>% suppressMessages()
library(stringr) %>% suppressMessages()
library(data.table) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(ggrepel) %>% suppressMessages()
suppressMessages(library(edgeR))
suppressMessages(library(sva))
suppressMessages(library(pheatmap))

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
# Read data
data <- fread(input_file, data.table = FALSE)
rownames(data) = data[[1]];data = data[,-1]

meta = fread(group_file, data.table = FALSE)
if (batch!='no'){
  check_col_names = colnames(meta)
  if (any(!sort(c("sample","group","tag","batch" )) == sort(check_col_names))){
    stop('Should have 4 columns in the sample information file：sample,group,tag,batch')
  }
}else {
  check_col_names = colnames(meta)
  if (any(!sort(c("sample","group","tag")) == sort(check_col_names[1:3]))){
    stop('Should have 3 columns in the grouping file：sample,group,tag')
  }
}

meta$sample = as.character(meta$sample)
meta$group = as.character(meta$group)
meta$tag = as.character(meta$tag)

if (any(duplicated(meta$tag))) {
  stop("Duplicates in tag column")
} 

if (any(!meta$sample %in% colnames(data))){
  error_message <- paste0(
    "Can not find the sample in the sample information file\n",
    paste(meta$sample[!meta$sample %in% colnames(data)], collapse = ", "), "\n"
  )
  stop(error_message)
}

# Sort 
data = data[,colnames(data) %in% meta$sample]
data = data[,meta$sample]

data = log2(data+1)

if (batch!='no'){
  if (!"batch" %in% colnames(meta)) {
    stop("Error: 'batch' column is missing in the meta data!")
  }
  # Check NA
  if (any(is.na(meta$batch))) {
    stop("Error: 'batch' column contains missing values!")
  }
  # Check empty
  if (length(meta$batch) == 0 || all(meta$batch == "")) {
    stop("Error: 'batch' column is empty!")
  }
  if (length(unique(meta$batch)) == 1) {
    stop("Error: 'More than 1 batch")
  }
  if (batch=='removeBatchEffect'){
    experiment_batch <- as.factor(meta$batch)
    data <- removeBatchEffect(data,batch=experiment_batch)
  }else if (batch=='combat'){
    experiment_batch <- as.factor(meta$batch)
    #gene_var <- apply(data, 1, var, na.rm = TRUE)
    #data <- na.omit(data)          # 删除含 NA 的行
    #data <- data[!apply(data, 1, function(x) any(is.infinite(x))), ]
    #data <- data[gene_var > 0, ]
    data = ComBat(dat = data, batch = experiment_batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
  }
}

data = na.omit(data)

if (nrow(data)==0){
  stop('Error after batch process, please check')
}

# Process
row.names(meta) = meta$sample
lables <- meta[colnames(data),'tag']
colors <- meta[colnames(data),'group']

###################################################################################################
################################                PLOT               ################################
###################################################################################################
options(scipen = 200)
data = na.omit(data)
if (ncol(data)==0){
  stop('Error after batch process, please check')
}
# PCA
pdf(paste0(output_directory,"/pca.pdf"),width=plot_width,height=plot_height)
pca <- prcomp(t(data))

var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1_pct <- round(var_explained[1], 1)
pc2_pct <- round(var_explained[2], 1)

pca_data <- data.frame(
  PC1     = pca$x[, "PC1"],
  PC2     = pca$x[, "PC2"],
  Label  = lables,
  Color = colors
)

# Plot
ggplot(pca_data, aes(x = PC1, y = PC2, label = Label, color = Color)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf, segment.linetype = "dashed") +
  labs(
    title = "PCA",
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
)
while (!is.null(dev.list()))  dev.off()


# MDS
d <- dist(t(data))
if (ncol(data) < 3) {
  message("Not enough samples, skip MDS")
} else{
	d <- dist(t(data)) 
	fit <- cmdscale(d, eig = TRUE, k = 2)
	x <- fit$points[,1]
	y <- fit$points[,2]

	pdf(paste0(output_directory,"/MDS.pdf"), width = plot_width, height = plot_height)
	mds_data <- data.frame(
	  c1 = x,
	  c2 = y,
	  Label  = lables,
	  Color = colors
	)
	# Plot
	ggplot(mds_data, aes(x = c1, y = c2, label = Label, color = Color)) +
	  geom_point(size = 3) +  
	  geom_text_repel(size = 3, max.overlaps = Inf,segment.linetype = "dashed") + 
	  labs(title = "MDS", x = "Coordinate 1", y = "Coordinate 2") +
	  theme_bw() +
	  theme(
	    panel.grid.major = element_blank(),  
	    panel.grid.minor = element_blank())
	while (!is.null(dev.list()))  dev.off()
}

# correlation
#colnames(data) = meta$tag
#sample_cor <- cor(data)
#pdf(paste0(output_directory,"/cor_heatmap.pdf"),width=plot_width,height=plot_height)
#ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"), bias=1)(100)
#ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"))(100)

#pheatmap(
#  sample_cor,
#  color            = ColorRamp,
#  display_numbers  = TRUE,        
#  number_format    = "%.2f",     
#  fontsize_number  = 8,        
#  cluster_rows     = T,      
#  cluster_cols     = T,
#  legend            = TRUE,      
#  legend_breaks    = c(0,0.5,1),  
#  legend_labels    = c("0","0.5","1"),
#  border_color     = NA,
#  show_rownames    = TRUE,
#  show_colnames    = TRUE
#)
#while (!is.null(dev.list()))  dev.off()

# 相关性热图
colnames(data) <- meta$tag
sample_cor <- cor(data)

# --- 检查数据有效性 ---
if (any(is.na(sample_cor))) {
  message("Correlation matrix contains NA, skip heatmap.")
} else if (length(unique(as.vector(sample_cor))) < 2) {
  message("All correlation values are identical, skip heatmap.")
} else {
  # --- 尝试绘图 ---
  tryCatch(
    expr = {
      # 打开 PDF 设备
      pdf(paste0(output_directory, "/cor_heatmap.pdf"), 
          width = plot_width, height = plot_height)
      
      # 颜色调色板（100 种颜色）
      ColorRamp <- colorRampPalette(c("white", "#F2FAB0", "red"))(100)
      breaks <- seq(
 	   min(sample_cor, na.rm = TRUE),
  	   max(sample_cor, na.rm = TRUE),
  	   length.out = length(ColorRamp) + 1
	)
      # 绘制热图（不手动指定 legend_breaks，避免断点重复问题）
      pheatmap(
        sample_cor,
        color            = ColorRamp,
	breaks           = breaks,
        display_numbers  = TRUE,
        number_format    = "%.4f",
        fontsize_number  = 8,
        cluster_rows     = TRUE,
        cluster_cols     = TRUE,
        legend           = TRUE,
        #legend_breaks  = c(0, 0.5, 1),   
        #legend_labels  = c("0", "0.5", "1"),
        border_color     = 'black',
        show_rownames    = TRUE,
        show_colnames    = TRUE
      )
      
      dev.off()
      message("Correlation heatmap saved successfully.")
    },
    error = function(e) {
      message("Correlation heatmap failed: ", e$message)
      # 确保关闭所有图形设备
      while (!is.null(dev.list())) dev.off()
    }
  )
}

cat('Done, you can check the results now\n')

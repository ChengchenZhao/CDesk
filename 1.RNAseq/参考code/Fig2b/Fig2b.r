###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)
library(VennDiagram)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
# 1,985 mechanical strain induced genes. 
MSIG <- as.vector(read.table("../../Figure1/Fig1c/MechanicalStrainInducedGenes.txt")[,1])

# GAR overlapped promoters
GARG <- as.vector(read.table("../../Figure2/Fig2d/GAP.OverlappedGenes.txt")[,1])

# all genes
all_genes <- as.vector(read.table("../../Figure2/Fig2d/AllGenes.txt")[,1])
###################################################################################################
################################                PLOT               ################################
###################################################################################################

GARG_MSIG_venn <- venn.diagram(list(GARG=GARG,MSIG=MSIG),cat.pos=1, fill=cccol[c(6,7)], alpha=c(0.5,0.5), cex=1, cat.cex=1, height = 1000, width = 1000, resolution = 500, filename=NULL)
pdf(file="GARG_MSIG_venn.pdf")
grid.draw(GARG_MSIG_venn)
dev.off()

# fisher.test(matrix(c(n - union(A,B), setdiff(A,B), setdiff(B,A), intersect(A,B)), nrow=2), nrow=2)
N = all_genes; A = GARG; B = MSIG
fisher.test(matrix(c(length(N) - length(union(A,B)), length(setdiff(A,B)), length(setdiff(B,A)), length(intersect(A,B))), nrow=2))$p.value
# 0.0000000000000000000000000000320011

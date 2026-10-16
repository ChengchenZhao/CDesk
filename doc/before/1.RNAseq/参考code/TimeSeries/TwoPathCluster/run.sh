Rscript Fig.ExpClusters.FPKM.r
Rscript Fig.ExpClusters.ZScore.r

for i in {1..12}
do
	Rscript ~/bin/GO_mm.r Zscore.cluster_${i}_genes.txt Zscore.cluster_${i}_genes &
done
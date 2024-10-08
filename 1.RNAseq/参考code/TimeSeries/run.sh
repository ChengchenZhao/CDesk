Rscript Fig.MarkerGeneExpressionPattern.r
Rscript Fig.SelectedGeneHeatmap.r
Rscript Fig.ExpClusters.r

for i in {1..9}
do
	Rscript ~/bin/GO_mm.r cluster_${i}_gene.mESCData.txt cluster_${i}_gene.mESC &
	Rscript ~/bin/GO_mm.r cluster_${i}_gene.tCLCData.txt cluster_${i}_gene.tCLC &
done
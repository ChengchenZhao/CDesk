GAP=/mnt/Storage/home/zhaochengchen/Work/3.MCNT/results/1.GradientAccessibility/GainAccessPeaks.bed
LAP=/mnt/Storage/home/zhaochengchen/Work/3.MCNT/results/1.GradientAccessibility/LostAccessPeaks.bed
RAP=/mnt/Storage/home/zhaochengchen/Work/3.MCNT/results/1.GradientAccessibility/ReservedAccessPeaks.bed

Promoter=/mnt/Storage/home/zhaochengchen/Data/mm10/mm10.Promoter.bed
# awk -v "OFS=\t" '{if (($6=="+") && ($2-5000>0)) print $1,$2-5000,$2+5000,$4,$5;else if (($6=="-") && ($3-5000>0)) print $1,$3-5000,$3+5000,$4,$5}' ~/Data/mm10/mm10.refseq.bed | sort -k1,1 -k2,2n | uniq > mm10.Promoter.5k.bed
# Promoter=mm10.Promoter.5k.bed
Enhancer_tmp=/mnt/Storage/home/zhaochengchen/Data/mm10/Yanghui_mouse_ATAC.bed
intersectBed -v -a ${Enhancer_tmp} -b ${Promoter} > /mnt/Storage/home/zhaochengchen/Data/mm10/mm10.PoterntialEnhancer.bed
Enhancer=/mnt/Storage/home/zhaochengchen/Data/mm10/mm10.PoterntialEnhancer.bed

wc -l ${GAP}
# 57823
intersectBed -wa -a ${GAP} -b ${Promoter} | sort | uniq | wc -l
# 8925
intersectBed -wa -a ${GAP} -b ${Enhancer} | sort | uniq | wc -l
# 44926


awk -v "OFS=\t" '{if (($6=="+") && ($2-10000>0)) print $1,$2-10000,$2-2000,$4,$5;else if ($6=="-") print $1,$3+2000,$3+10000,$4,$5}' ~/Data/mm10/mm10.refseq.bed | sort -k1,1 -k2,2n | uniq > ~/Data/mm10/mm10.DistalRegion.bed

intersectBed -wa -b ${GAP} -a ${Promoter} | cut -f 5 | sort | uniq > GAP.OverlappedGenes.txt

intersectBed -wb -a ${GAP} -b ${Enhancer} | sort | uniq > GAPOverlappedEnhancers.bed
intersectBed -wa -b GAPOverlappedEnhancers.bed -a ~/Data/mm10/mm10.DistalRegion.bed | cut -f 5 | sort | uniq > GAP.OverlappedGenes.Enhancer.txt

# Rscript ~/bin/GO.r GAP.OverlappedGenes.txt GAP.OverlappedGenes
# Rscript ~/bin/GO.r GAP.OverlappedGenes.Enhancer.txt GAP.OverlappedGenes.Enhancer

cat GAP.OverlappedGenes.txt GAP.OverlappedGenes.Enhancer.txt | sort | uniq > GAP.OverlappedGenes.PE.txt

cat CC_14h_up_genes.txt CC_14h_dn_genes.txt > CC_14h_DE_genes.txt
wc -l CC_14h_DE_genes.txt
cat CC_14h_DE_genes.txt GAP.OverlappedGenes.PE.txt | sort | uniq -d | wc -l

wc -l CC_14h_up_genes.txt
cat CC_14h_up_genes.txt GAP.OverlappedGenes.PE.txt | sort | uniq -d | wc -l
cat CC_14h_up_genes.txt GAP.OverlappedGenes.txt | sort | uniq -d | wc -l
cat CC_14h_up_genes.txt GAP.OverlappedGenes.Enhancer.txt | sort | uniq -d | wc -l
cat CC_14h_up_genes.txt GAP.OverlappedGenes.BothPromoterEnhancer.txt | sort | uniq -d | wc -l


# Super TDU
peak_path=/mnt/Storage/home/zhaochengchen/Work/3.MCNT/data/inhouse_ATACseq/0.mapping/peak
bw_path=/mnt/Storage/home/zhaochengchen/Work/3.MCNT/data/inhouse_ATACseq/0.mapping/bw
Promoter=/mnt/Storage/home/zhaochengchen/Data/mm10/mm10.Promoter.bed
Genebody=/mnt/Storage/home/zhaochengchen/Data/mm10/mm10.Genebody.bed

wc -l ${peak_path}/CC_peaks.bed
wc -l ${peak_path}/SCC_Ctrl_peaks.bed
wc -l ${peak_path}/M14h_peaks.bed
wc -l ${peak_path}/SCC_M14_peaks.bed
# 38228 /mnt/Storage/home/zhaochengchen/Work/3.MCNT/data/inhouse_ATACseq/0.mapping/peak/CC_peaks.bed
# 39376 /mnt/Storage/home/zhaochengchen/Work/3.MCNT/data/inhouse_ATACseq/0.mapping/peak/SCC_Ctrl_peaks.bed
# 92419 /mnt/Storage/home/zhaochengchen/Work/3.MCNT/data/inhouse_ATACseq/0.mapping/peak/M14h_peaks.bed
# 68394 /mnt/Storage/home/zhaochengchen/Work/3.MCNT/data/inhouse_ATACseq/0.mapping/peak/SCC_M14_peaks.bed


GAP=/mnt/Storage/home/zhaochengchen/Work/3.MCNT/results/1.GradientAccessibility/GainAccessPeaks.bed
LAP=/mnt/Storage/home/zhaochengchen/Work/3.MCNT/results/1.GradientAccessibility/LostAccessPeaks.bed
RAP=/mnt/Storage/home/zhaochengchen/Work/3.MCNT/results/1.GradientAccessibility/ReservedAccessPeaks.bed
intersectBed -v -a ${peak_path}/M14h_peaks.bed -b ${peak_path}/SCC_M14_peaks.bed | sort | uniq > SLost_Accessibility.bed
wc -l SLost_Accessibility.bed
intersectBed -wa -a SLost_Accessibility.bed -b ${GAP} | sort | uniq | wc -l
intersectBed -wa -a SLost_Accessibility.bed -b ${RAP} | sort | uniq | wc -l
# 33605 wc -l SLost_Accessibility.bed
# 30294 intersectBed -wa -a SLost_Accessibility.bed -b ${GAP} | sort | uniq | wc -l
# 3311 intersectBed -wa -a SLost_Accessibility.bed -b ${RAP} | sort | uniq | wc -l

wc -l ${GAP}
# 57823
wc -l ${RAP}
# 34596

intersectBed -wa -a SLost_Accessibility.bed -b ${GAP} | sort | uniq > SLost_Accessibility.GAP.bed
intersectBed -v -a ${GAP} -b SLost_Accessibility.GAP.bed | sort | uniq > SRemain_Accessibility.GAP.bed

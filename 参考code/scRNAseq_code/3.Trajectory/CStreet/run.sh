conda activate CStreet

data_path=/data/zhaochengchen/Work/2.CEmbryo/results/2.Clustering/Seurat

for stage in EGKXI_1h_WE EGKXI_2h_WE EGKXII_3h_WE EGKXII_4h_WE EGKXIII_5h_WE HH2_6h_WE HH2a_10h_WE HH3_13h_WE HH3a_16h_WE HH4_19h_WE HH5_22h_WE HH7_23h_WE HH8_26h_WE HH9_29h_WE HH10_33h_WE HH14_50h_SN HH14_50h_SO HH18_66h_SO HH20_70h_SW HH20_70h_WB HH21_82h_WB HH23_94h_WB HH24_106h_WB
do
    CStreet -i ${data_path}/ChickEmbryo_EGKX_0h_WE_RawCountsMatrix.txt ${data_path}/ChickEmbryo_${stage}_RawCountsMatrix.txt \
            -s ${data_path}/EGKX_0h_WE_AllSplitClusteringResult.txt ${data_path}/${stage}_AllSplitClusteringResult.txt \
            -T --Switch_Normalize ON --Switch_LogTransform ON --FigureParam_FigureSize 7 7 --Threshold_MinProbability 0.1 --ProbParam_SamplingSize 50 \
            -n CEmbryo_${stage}
    awk -v "OFS=\t" '{if ($3>0.1) printf "%s->%s\n",$1,$2}' CEmbryo_${stage}/CEmbryo_${stage}_CellStatesConnCytoscape.txt >> CellStateTrajectory.0h.txt
done
sort CellStateTrajectory.0h.txt | uniq -c > CellStateTrajectory.0h.Stat.txt
Rscript CellStateTrajectory.0h.Stat.r

CStreet -i ${data_path}/ChickEmbryo_EGKX_0h_WE_RawCountsMatrix.txt \
           ${data_path}/ChickEmbryo_EGKXI_1h_WE_RawCountsMatrix.txt \
           ${data_path}/ChickEmbryo_EGKXI_2h_WE_RawCountsMatrix.txt \
           ${data_path}/ChickEmbryo_EGKXII_3h_WE_RawCountsMatrix.txt \
           ${data_path}/ChickEmbryo_EGKXII_4h_WE_RawCountsMatrix.txt \
           ${data_path}/ChickEmbryo_EGKXIII_5h_WE_RawCountsMatrix.txt \
           ${data_path}/ChickEmbryo_HH2_6h_WE_RawCountsMatrix.txt \
           ${data_path}/ChickEmbryo_HH2a_10h_WE_RawCountsMatrix.txt \
           ${data_path}/ChickEmbryo_HH3_13h_WE_RawCountsMatrix.txt \
        -s ${data_path}/EGKX_0h_WE_AllSplitClusteringResult.txt \
           ${data_path}/EGKXI_1h_WE_AllSplitClusteringResult.txt \
           ${data_path}/EGKXI_2h_WE_AllSplitClusteringResult.txt \
           ${data_path}/EGKXII_3h_WE_AllSplitClusteringResult.txt \
           ${data_path}/EGKXII_4h_WE_AllSplitClusteringResult.txt \
           ${data_path}/EGKXIII_5h_WE_AllSplitClusteringResult.txt \
           ${data_path}/HH2_6h_WE_AllSplitClusteringResult.txt \
           ${data_path}/HH2a_10h_WE_AllSplitClusteringResult.txt \
           ${data_path}/HH3_13h_WE_AllSplitClusteringResult.txt \
        -T --Switch_Normalize ON --Switch_LogTransform ON --FigureParam_FigureSize 24 7 --Threshold_MinProbability 0.1 --ProbParam_SamplingSize 50 -n CEmbryo_All &

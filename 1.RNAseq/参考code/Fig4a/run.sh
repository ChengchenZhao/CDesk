####################################################################################################
#######################            SCNT Falsely Activated/Repressed            #####################
####################################################################################################
Rscript ExpressionVolcanoPlot.r

cat Ctrl_1C_dn_genes.txt Ctrl_4C_dn_genes.txt Ctrl_E2C_dn_genes.txt Ctrl_L2C_dn_genes.txt | sort | uniq -c > Ctrl_All_dn_genes.txt
cat Ctrl_1C_up_genes.txt Ctrl_4C_up_genes.txt Ctrl_E2C_up_genes.txt Ctrl_L2C_up_genes.txt | sort | uniq -c > Ctrl_All_up_genes.txt
cat M14h_1C_dn_genes.txt M14h_4C_dn_genes.txt M14h_E2C_dn_genes.txt M14h_L2C_dn_genes.txt | sort | uniq -c > M14h_All_dn_genes.txt
cat M14h_1C_up_genes.txt M14h_4C_up_genes.txt M14h_E2C_up_genes.txt M14h_L2C_up_genes.txt | sort | uniq -c > M14h_All_up_genes.txt

wc -l Ctrl_All_dn_genes.txt >> SCNT_FalaselyExpressedGeneNumber.txt
wc -l Ctrl_All_up_genes.txt >> SCNT_FalaselyExpressedGeneNumber.txt
wc -l M14h_All_dn_genes.txt >> SCNT_FalaselyExpressedGeneNumber.txt
wc -l M14h_All_up_genes.txt >> SCNT_FalaselyExpressedGeneNumber.txt
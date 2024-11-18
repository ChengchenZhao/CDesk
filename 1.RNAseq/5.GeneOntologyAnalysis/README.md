# packages required

* ggplot2
* clusterProfiler
* dplyr



# overview of the functions and usage

## enrichment_analysis_customer.R

> Discription: A R script for performing the enrichment analysis based on the customization gene sets

```bash
#! /usr/bin/bash

/usr/bin/ Rscript enrichment_analysis_customer.R <gene_file(SYMBOL)> <customer_term(term    symbol)>
```

| arguments                       | discription                                                  |
| ------------------------------- | ------------------------------------------------------------ |
| <gene_file(SYMBOL)>             | path  of  a single column file containing genes of  interests(SYMBOL) |
| <customer_term(term    symbol)> | A customer file containing two columns, every row is consisted of two factors: the customer term name and gene of interests separated by tab |

### output example

```bash
ID      Description     GeneRatio       BgRatio pvalue  p.adjust        qvalue  geneID  Count
term_172        hsa04310        40/638  174/8850        2.53440393878609e-11    8.1607806828912e-09     6.34934881506410e-09    9350/2239/6422/22943/27121/4920/85407/55366/4919/7474/4163/83439/7090/59343/23411/144165/7855/81839/8549/64321/6932/23401/25805/6934/8325/7473/51701/817/818/84870/4088/147495/5579/51176/5176/147111/8321/7478/81029/11197     40
term_231        hsa04820        47/638  232/8850        5.10757594904185e-11    8.22319727795738e-09    6.3979109256419e-09     1462/7431/7414/1291/2201/2192/7168/1605/375790/1289/1832/3655/3678/1277/1282/3685/1290/10398/29766/6383/481/287/6385/7057/288/1299/3679/3680/2273/7058/58/70/478/2274/89/5318/1298/7273/3676/7111/6641/10529/1297/72/22989/4620/7139    47
term_188        hsa04550        33/638  144/8850        1.58546336506169e-09    1.39250161555254e-07    1.08341086793561e-07    2260/3720/2247/7547/207/6657/657/5604/7474/3399/3977/5291/4838/7855/2263/6932/3572/659/91/3398/8325/7473/6926/5295/652/5293/4088/83729/8321/7478/9421/81029/3624        33
term_183        hsa04512        25/638  89/8850 1.81028955252498e-09    1.39250161555254e-07    1.08341086793561e-07    1291/1605/375790/3912/3655/3911/3678/1277/1282/3685/3913/6385/7057/6696/1299/3679/3680/7058/5649/1298/3676/3918/158326/3371/1297        25
term_176        hsa04360        38/638  184/8850        2.16226958936729e-09    1.39250161555254e-07    1.08341086793561e-07    57556/5362/5364/1969/1808/7474/10501/6608/9423/2043/137970/10398/5291/4756/2770/6405/659/219699/9037/2049/10507/10371/5295/817/818/80031/5293/4233/7223/2051/54437/59277/7852/64218/8633/84628/57689/81029      38
term_333        hsa05217        19/638  63/8850 4.61205980570961e-08    2.47513876239749e-06    1.92573725220857e-06    7474/83439/6608/7855/650/6932/2736/6934/8325/7473/652/51684/2735/51176/8321/7478/10912/81029/64399      19
term_160        hsa04151        55/638  362/8850        7.40781986675929e-08    3.40759713870927e-06    2.65121974178753e-06    1291/2260/2064/2247/3912/8822/207/3655/5156/54541/3911/3678/1277/5604/1282/1969/3685/5106/9223/5291/2263/3913/2321/2790/598/23239/7057/6696/1299/2066/2065/3679/3680/7058/2249/3643/7010/5295/5649/55970/5521/5293/1298/6446/3676/7039/23533/3918/3371/4233/285/1297/10161/5154/627     55
term_182        hsa04510        37/638  203/8850        1.14286636239406e-07    4.60003710863607e-06    3.57897624012875e-06    7414/1291/2317/2064/3912/207/3655/5156/3911/3678/1277/5604/1282/3685/10398/5291/3913/2321/7057/6696/824/1299/3679/3680/7058/5295/5649/5293/1298/3676/3918/857/3371/4233/5579/1297/5154  37
term_322        hsa05205        36/638  204/8850        3.89752652841402e-07    1.3944483801659e-05     1.08492551317256e-05    2260/2064/2247/207/3678/1277/5604/3685/7474/6608/6383/4313/5291/7855/287/3709/6385/2719/7057/288/2066/2065/8325/7473/5295/817/818/5293/857/4233/5579/8321/7478/355/81029/286    36
term_270        hsa04974        23/638  105/8850        1.13982027273613e-06    3.67022127821034e-05    2.85554973590736e-05    54407/1291/6510/6520/80781/1289/23428/1277/1282/1290/481/4311/1305/1299/7373/136227/8645/478/1298/1294/1297/1303/81578  23
```



## enrichment_analysis.R

> Discription: A R script for performing the GO/KEGG enrichment analysis based on the pubsihed GO/KEGG gene sets

```bash
#! /usr/bin/bash

# for GO analysis
/usr/bin/Rscript enrichment_analysis.R <gene_file(SYMBOL)> GO <species> <similarity cutoff level[0:1]>

# for KEGG analysis
/usr/bin/Rscript enrichment_analysis.R <gene_file(SYMBOL)> KEGG <species>  <organism>
```

| arguments                      | discription                                                  |
| :----------------------------- | ------------------------------------------------------------ |
| <gene_file(SYMBOL)>            | path  of  a single column file containing genes of  interests(SYMBOL) |
| < species >                    | Chose the functional annotation database of specific species: human, rat, mouse, chicken, or pig |
| <similarity cutoff level[0:1]> | Define a threshold between 0~1, the GO terms will be merged into one once the gene sets matched to them exceed the threshold |
| < organism>                    | Chose the organism for KEGG analysis                         |

### Out put example

```bash
ONTOLOGY        ID      Description     GeneRatio       BgRatio pvalue  p.adjust        qvalue  geneID  Count
BP      GO:0050673      epithelial cell proliferation   71/1114 481/18903       6.63909688535626e-13    7.8284840294988e-11     5.75884036662593e-11    APLNR/HAS2/FGFR1/SFRP1/RBPJ/MMP14/ERBB2/TRIM24/FGF2/LAMB1/AKT1/SOX2/BMPR1A/CDH3/CNMD/LAMA5/LGR4/MAP2K1/EPHA2/WNT5A/MTA3/MCC/SMO/APELA/SIRT1/PIK3CB/NODAL/FGFR2/FLT1/LIPA/VASH2/LGR5/BCL2L1/BMP2/GPC3/THBS1/BMPR2/TCF7L2/ERRFI1/SEMA5A/HMOX1/FST/EDNRB/ID2/ZNF703/CFLAR/ESRP2/TEK/BMP4/BCL2L2/MARVELD3/PTPRM/PIK3CD/SMAD3/ITGA4/TGFA/BID/CAV1/IGFBP3/PTPRN/VASH1/HTRA1/GLI1/HES1/BMPER/SERPINF1/SNAI2/ST8SIA1/ZEB1/TNFSF11/CDH13 71
BP      GO:0097106      postsynaptic density organization       8/1114  31/18903        0.00033108613801222     0.0035286141862764      0.00259574213060998     CDH2/SYNGAP1/LRFN4/SHANK2/RELN/NPTX1/CBLN1/NRXN1        8
BP      GO:0060068      vagina development      4/1114  10/18903        0.00189096059228965     0.0144947491043835      0.0106627216625724      LHX1/WNT5A/TYRO3/LRP2   4
BP      GO:0030856      regulation of epithelial cell differentiation   21/1114 164/18903       0.000660480424665218    0.00621457133021161     0.00457160340403608     CD24/VCL/LHX1/FGF2/ESRP1/PROM1/SMO/NODAL/FOXC1/BMP2/S1PR3/SPRED2/ERRFI1/FST/TBX3/BMP4/MESP1/CAV1/HES1/MSX2/ZEB1 21
CC      GO:0005902      microvillus     15/1153 93/19869        0.000285135591355247    0.00282080567162155     0.00215995570143166     FSCN1/PODXL/SLC7A5/SLC6A6/SLC7A8/PDGFRA/ITGAV/PROM1/PDPN/IQGAP2/TEK/VIL1/PDGFA/CUBN/MYO7A       15
BP      GO:0001890      placenta development    30/1114 150/18903       3.0249367133784e-09     1.67985015670836e-07    1.23574230411432e-07    EOMES/RBPJ/BPTF/HS6ST1/ADAM19/AKT1/MAP2K1/STC2/APELA/STC1/NODAL/FZD5/FGFR2/VASH2/MME/SPP1/BMPR2/SNAI1/MDFI/EGLN1/RSPO3/ST14/VASH1/HTRA1/HES1/LEF1/LHX4/HAND1/SOX15/PRDM1        30

```

## **plot_enrichment.R** 

> Discription: A R script for drawing the bubble plot and bar plot for the enrichment information generated by the enrichment_analysis_customer.R or enrichment_analysis.R

```bash
#! /usr/bin/bash

/usr/bin/Rscript plot_enrichment.R <enrichment_file> <color for high value> <color for low value>
```

| arguments               | discription                                                 |
| :---------------------- | ----------------------------------------------------------- |
| <enrichment_file>       | Enrichment information of interests which can be            |
| < color for high value> | Assign a color representing the maximum –lg(p.adjust) value |
| < color for low value>  | Assign a color representing the minim –lg(p.adjust) value   |

### output example

![image-20241117220922583](C:\Users\ztzhang\AppData\Roaming\Typora\typora-user-images\image-20241117220922583.png)
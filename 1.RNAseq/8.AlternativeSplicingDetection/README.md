# Softwares installation

```bash
# install the rMATs

conda install -c bioconda rmats

# install the rmats2sashimiplot

conda install -c bioconda rmats2sashimiplot
```

# overview of the functions and usage

## enrichment_analysis_customer.R

> Discription: A R script for detecting the differentail splicing events by rMATs based on the bam files.

```bash
#! /usr/bin/bash

/usr/bin/bash alternative_splicing.sh --b1 <b1 file path>  --b2 <b2 file path> --gtf <GTF file path> --od <output directory>

```

| arguments           | discription                                                  |
| ------------------- | ------------------------------------------------------------ |
| < b1 file path>     | Path  to the b1 file containing BAM file paths separated by comma |
| < b2  file path>    | Path  to the b2 file containing BAM file paths separated by comma |
| < GTF file  path >  | Path  to the GTF file                                        |
| < output directory> | Output directory for the results                             |

**the explaination of output: rmats-turbo/README.md at v4.1.2 · Xinglab/rmats-turbo (github.com)**

## **draw_sashimiplot.sh** 

> Discription: A R script for drawing the sashimi plot to visualize the differential splicing events based on the bam files..

```bash
#! /usr/bin/bash

/usr/bin/bash draw_sashimiplot.sh --b1 <b1 file path> --b2 <b2 file path> --gff <GFF file path> --od <output directory> --r <interval> --g <group file path> 

```

| arguments           | discription                                                  |
| ------------------- | ------------------------------------------------------------ |
| < b1 file path>     | Path  to the b1 file containing BAM file paths separated by comma |
| < b2  file path>    | Path  to the b2 file containing BAM file paths separated by comma |
| < GTF file  path >  | Path  to the GTF file                                        |
| < output directory> | Output directory for the results                             |
| < interval>         | Interval  of interests,chromosome:orientation:start:end.  example: chr11:-:4752513:4854957 |
| <  group file path> | Path to the file containing group information of samples. example: WT: 1-3 |

### output example

![image-20241117224916907](C:\Users\ztzhang\AppData\Roaming\Typora\typora-user-images\image-20241117224916907.png)
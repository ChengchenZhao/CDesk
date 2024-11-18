### Language:
```
#!python
```

### Enviroment:
The DEGD.py script uses a few python packages. You can create a conda environment by following commands below:
```
conda create --name degd
conda activate degd
conda install numpy, pandas, matplotlib
pip install seaborn, omicverse
```
Now you can easily use this script by following the commands below in the terminal.Parameters includes:

| 参数         | 描述                                                 | 默认值    |
|--------------|-----------------------------------------------------|-----------|
| ``           | deseq2 gfold or adjust_t                            | /         |
| `--input`    | Input count matrix for DESeq2/adjusted t-test<br> or Input bam file directory for GFOLD                      | /         |
| `--output_dir` | Directory for saving output                       | /         |
| `--species`  | Species for choosing .gtf for mapping gene id to gene name | `mouse`       |
| `--meta_file`  | Group informations for analysis                   | /         |
| `--plot`    | Whether to draw plot or not                          | /         |
| `--heatmap_types`   | Data types for drawing heatmap plots         | /         |
| `--vocano_types`  | Data types for drawing vocano plots            | `qvalue`  |
| `--gene_of_interest`  | Interest gene for plot heatmap             | /         |
| `--top_genes`    | Number of top differential expression genes for plot heatmap     | 30        |
| `--fc_threshold`   | fc_threshold for vocano plot                  | 1.5       |

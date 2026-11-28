<div align="center">

# CDesk Multi-Omics Analysis Pipeline

📚 **Documentation**: [**User Guide**](https://zsifu.github.io/CDesk.html)

</div>


CDesk is an integrated multi-omics analysis sand-box designed for processing data from various sequencing-based assays, including RNA-seq, scRNA-seq, ATAC-seq, CUT&Tag, ChIP-seq, and Hi-C. It comprises multiple modules that cover a comprehensive range of analysis tasks, from raw sequencing data process to downstream various advanced functions. Dedicated conda environment YAML files are supplied. To install, simply download the pre-built Docker image and prepare the necessary species-specific data. Once configured, users can perform the desired analyses by entering the corresponding command on the command line.

<div align="center">

![](./doc/images/CDesk.png)
</div>

## Installation

### Download the Docker Image

CDesk and its required tools, Conda environments, R packages, and Python packages are bundled in the Docker image. Download and load the image **【待填写：最新镜像下载链接】** to use CDesk.

### Load and Start the Image

#### Linux

Install Docker Engine, then load the downloaded archive(this step will take a few minutes):

~~~bash
docker load -i "[path to the archive]"
~~~

Start the container with a working directory:

~~~bash
docker run -it --name [container name] \
  -v "[working directory of computer(host)]:/[working directory name of container]" \
  -w /[working directory name of container] \
  "cdesk:v8" \
  /bin/bash
~~~

Keep the image's /CDesk directory intact; do not mount your working directory over it.

#### Windows

Install WSL 2 with a Linux distribution such as Ubuntu, and install Docker Desktop. In Docker Desktop, enable the WSL 2 engine and WSL integration for that distribution. Run the following commands in the **WSL terminal**.

Load the downloaded image. 

~~~bash
docker load -i "[path to the archive]"
~~~

Start the container with a working directory:

~~~bash
docker run -it --name [container name] \
  -v "[working directory of computer(host)]:/[working directory name of container]" \
  -w /[working directory name of container] \
  "cdesk:v8" \
  /bin/bash
~~~


After exiting the container, Linux and Windows/WSL users can return to the same container with:

~~~bash
docker start -ai [container name]
~~~

### After Installation, You Should Be Able to Run CDesk Now

Inside the container, run:

~~~bash
/CDesk/CDesk/CDesk -h
~~~

If the CDesk help is displayed, the installation is successful. The available modules include:

~~~text
usage: CDesk [-h] {tools,bulkRNA,scRNA,ATAC,ChIPseqCUTTag,HiC} ...

CDesk multiomics pipeline

positional arguments:
  {tools,bulkRNA,scRNA,ATAC,ChIPseqCUTTag,HiC}
                        Choose an omic to run the script
    tools               Some other tools
    bulkRNA             bulkRNA-seq pipeline
    scRNA               scRNA pipeline
    ATAC                ATAC pipeline
    ChIPseqCUTTag       ChIPseq&CUTTag pipeline
    HiC                 HiC pipeline

options:
  -h, --help            show this help message and exit
  --config CONFIG       The configuration files
~~~

### Prepare the Data and Write the Configuration File

Prepare reference data for the species of interest and specify the corresponding file paths in the config.json configuration file. This file stores paths to genomic data, annotation files, and additional software installations.

The example below is from the current config.json recording the reference data of the Arabidopsis thaliana. You can customize it to support additional species or data types. The reference data and tools needed depend on the task you intend to run.

<details>
<summary>An example of the configuration file</summary>

~~~json
{
  "software": {
    "cellranger": "/root/cellranger-10.0.0/bin/cellranger",
    "DrSeq": "/root/Dr.seq.1.2.0/bin/Drseq.py",
    "dnbc4tools": "/root/dnbc4tools3.1/dnbc4tools",
    "juicer_tools_jar": "/root/juicer_tools_1.22.01.jar",
    "celescope_path": "/opt/conda/envs/celescope/bin"
  },
  "conda_env": {
    "CDesk": "/opt/conda/envs/CDesk",
    "CDesk_R": "/opt/conda/envs/CDesk_R",
    "CDesk_py3.7": "/opt/conda/envs/CDesk_py3.7",
    "CDesk_py2.7": "/opt/conda/envs/CDesk_py2.7"
  },
  "data": {
    "TAIR10": {
      "hisat2_index": "/CDesk/TAIR/hisat2_index/TAIR10",
      "refseq_gtf": "/CDesk/TAIR/Araport11.gtf",
      "refseq_bed": "/CDesk/TAIR/Araport11.bed",
      "chromInfo": "/CDesk/TAIR/TAIR10.chrom.sizes",
      "fasta": "/CDesk/TAIR/TAIR10.fasta",
      "tf_file": "/CDesk/TAIR/Ath_TF_list.txt",
      "promoter_file": "/CDesk/TAIR/TAIR10_promoter_transcipt_updown_2kb_promoter.bed",
      "TE_idx": "/CDesk/TAIR/TAIR10.exclusive.idx",
      "scRef10x": "/CDesk/TAIR/scRef10x/TAIR",
      "effective_genome_size": "1.19e8",
      "gff3": "/CDesk/TAIR/Araport11.gff",
      "refgenes": "/CDesk/TAIR/TAIR.refgenes.txt",
      "bowtie2_mapindex": "/CDesk/TAIR/bowtie2_index/TAIR10",
      "singleron_mapindex": "/CDesk/TAIR/singleron_mapindex",
      "dnbc_mapindex": "/CDesk/TAIR/dnbcIndex/TAIR"
    }
  }
}

~~~

</details>

- **software**: Paths to tools installed in the image. Keep the entries required by your workflows.
- **conda_env**: Paths to the CDesk environments inside the image.
- **data**: Reference data configurations for different species. The TAIR10 key corresponds to --species TAIR10.

The reference data fields are described below.

| Field | Meaning and Purpose | What to Specify |
|---|---|---|
| `hisat2_index` | HISAT2 alignment index used by workflows such as bulk RNA-seq to align sequencing reads to the reference genome | The shared prefix of the index files, not a directory or a FASTA file |
| `refseq_gtf` | GTF gene annotation containing the coordinates of genes, transcripts, exons, and other features; used for RNA analysis and for locating genes in `accessbility --mode gene` | `.gtf` file |
| `refseq_bed` | BED annotation of transcript structures, read as an RSeQC annotation during bulk RNA preprocessing for read distribution quality control | An annotation BED file that meets the tool's requirements, typically BED12; an arbitrary three-column peak BED file cannot be used as a substitute |
| `chromInfo` | Chromosome size table containing chromosome names and lengths; used for signal statistics, BigWig conversion, and related tasks | A text file, typically with two columns |
| `fasta` | Reference genome sequences containing the DNA nucleotide sequences of chromosomes | `.fa`, `.fasta`, or `.fna` file |
| `tf_file` | Table of transcription factor genes, used by ATAC dynamic analysis to identify transcription factors among associated genes | The current code reads a tab-delimited table and requires a column named `Symbol` |
| `promoter_file` | Table of promoter regions, used to calculate promoter signals and associate ATAC regions with genes | The relevant current code requires at least five columns: chromosome, start, end, transcript, and gene |
| `TE_idx` | Index for transposable element analysis; the bulk RNA TE analysis workflow passes it to `scTE -x` | An index file recognized by scTE, not a HISAT2 or Bowtie2 index |
| `scRef10x` | Cell Ranger reference database used for preprocessing 10x single-cell RNA-seq data | A prebuilt Cell Ranger reference directory, not a single FASTA file |
| `effective_genome_size` | Effective genome size used in MACS peak-calling calculations for ATAC/ChIPseqCUTTag | A numeric value or another value supported by the tool, not a file path |
| `gff3` | GFF3 gene annotation read by the current alternative splicing visualization script | `.gff3` file; changing a GTF file's extension does not convert it to GFF3 |
| `refgenes` | Reference gene annotation used by DrSeq; the source code passes it to DrSeq's `-g` argument | A reference annotation file that meets DrSeq's requirements, not a simple single-column list of gene names |
| `bowtie2_mapindex` | Bowtie2 alignment index used by workflows such as ATAC, ChIPseqCUTTag, and DrSeq | The shared prefix of the index files |
| `singleron_mapindex` | Singleron/CeleScope reference database; the source code passes it to `--genomeDir` | A reference directory prepared for the corresponding tool |
| `dnbc_mapindex` | DNBelab C Series/dnbc4tools reference database; the source code passes it to `--genomeDir` | A reference directory prepared for the corresponding tool |

Use consistent chromosome names and a matching genome assembly across reference sequences, annotations, indices, and analysis inputs.

#### Use a Configuration File Stored Outside the Image

You can create and edit customized config.json in the host working directory without modifying the image. Copy the complete configuration example above, retain the software and conda_env paths for the image you are using, and edit the reference data entries as needed.

The host configuration file is available inside the container . CDesk uses its bundled config.json by default; select your external file explicitly with **--config before the module name**:

~~~bash
CDesk \
  --config [path to customized config.json] \
  ATAC preprocess \
  --input [path to samples.csv] \
  --output [output directory] \
  --species mm10
~~~


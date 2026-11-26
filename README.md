<div align="center">

# CDesk Multi-Omics Analysis Pipeline

📚 **Documentation**: [**User Guide**](https://zsifu.github.io/CDesk.html)

</div>


CDesk is an integrated multi-omics analysis sand-box designed for processing data from various sequencing-based assays, including RNA-seq, scRNA-seq, ATAC-seq, CUT&Tag, ChIP-seq, and Hi-C. It comprises multiple modules that cover a comprehensive range of analysis tasks, from raw sequencing data process to downstream various advanced functions. Dedicated conda environment YAML files are supplied. To install, simply download the pre-built Docker image and prepare the necessary species-specific data. Once configured, users can perform the desired analyses by entering the corresponding command on the command line.

<div align="center">

![](./doc/images/CDesk.png)
</div>

## Installation

### 1. Choose an installation method

Use a pre-built CDesk Docker image for routine analysis. The image provides the software environment; input data and analysis results are stored in directories mounted from the host. Linux users can use Docker Engine. Windows and macOS users can use Docker Desktop with Linux containers, provided that the image supports their CPU architecture.

For the code under development in this branch, use **the development source inside a compatible CDesk image** as described in section 3. Loading an image does not update its bundled CDesk code to this branch.

**Image release status:** the previous README linked to [Zenodo](https://doi.org/10.5281/zenodo.19709171), but compatibility of that archive with the current development branch has not been established. A current release needs to identify the download URL, archive filename and SHA256 checksum, loaded image name and tag (or registry digest), supported architecture, bundled reference data, and corresponding code commit. Do not assume that the old archive name or image tag identifies the current branch. The image names in the commands below are placeholders.

Installing directly from Conda environment files is an advanced alternative (section 6). The environment definitions alone do not include every separately distributed tool or reference dataset.

### 2. Run a pre-built image

Install Docker and prepare host directories for the input data and output results. Replace every example path below with an existing absolute path on your machine.

After obtaining the appropriate image archive and its published checksum, verify the archive and load it:

~~~bash
# Linux; compare the output with the checksum published for this exact archive.
sha256sum /absolute/path/to/cdesk-image.tar.gz

# On macOS, use: shasum -a 256 /absolute/path/to/cdesk-image.tar.gz
# Docker can load a gzip-compressed image archive directly.
docker load -i /absolute/path/to/cdesk-image.tar.gz
~~~

Use the image name and tag printed by the load command. The value below is an example placeholder, not a published image name.

~~~bash
CDESK_IMAGE='IMAGE_NAME:TAG'

docker run --rm -it \
  --mount type=bind,src=/absolute/path/to/input,dst=/data,readonly \
  --mount type=bind,src=/absolute/path/to/results,dst=/results \
  --workdir /results \
  "$CDESK_IMAGE" \
  /bin/bash
~~~

Inside the container:

~~~bash
CDesk -h
~~~

This command uses the CDesk code bundled in the image. To use this development branch, follow section 3 instead.

- Input files are available under **/data** in the container. Use these container-visible paths in sample information CSV files, for example, /data/sample_R1.fastq.gz.
- Write analysis outputs under **/results** so that they remain available on the host after the container exits.
- Do not mount a host source directory over **/CDesk**. It would hide the image's existing files at that location, including any reference data under /CDesk/TAIR.
- If additional reference data are needed, mount their host directory separately, for example at /reference with a read-only bind mount, and use /reference paths in the configuration.
- Host directories used with --mount must already exist. On a cluster, run analysis on an allocated compute node using the container runtime permitted by the cluster.

### 3. Run this development branch inside an image

Get the source on the host:

~~~bash
git clone --branch Zsifu_edit2026 --single-branch \
  https://github.com/ChengchenZhao/CDesk.git
cd CDesk

# Several shell scripts are invoked directly by the dispatcher.
# Set execution permissions in this local checkout before mounting it read-only.
chmod +x CDesk
find src -type f -name '*.sh' -exec chmod +x {} \;

# Record the exact source revision used for an analysis.
git rev-parse HEAD
~~~

Mount the checkout at a location separate from the image's bundled code and references:

~~~bash
CDESK_IMAGE='IMAGE_NAME:TAG'

docker run --rm -it \
  --mount type=bind,src=/absolute/path/to/CDesk,dst=/opt/cdesk-dev,readonly \
  --mount type=bind,src=/absolute/path/to/input,dst=/data,readonly \
  --mount type=bind,src=/absolute/path/to/results,dst=/results \
  --workdir /results \
  "$CDESK_IMAGE" \
  /bin/bash
~~~

Inside the container, copy the branch configuration to the results directory. Edit that copy if the image's tool, environment, or reference paths differ:

~~~bash
cp /opt/cdesk-dev/config.json /results/config.json

python3 /opt/cdesk-dev/CDesk \
  --config /results/config.json \
  -h

python3 /opt/cdesk-dev/CDesk \
  --config /results/config.json \
  ATAC preprocess -h
~~~

For an analysis, keep the same explicit entry point and place --config **before** the module name:

~~~text
python3 /opt/cdesk-dev/CDesk --config /results/config.json MODULE COMMAND [arguments]
~~~

Replace MODULE, COMMAND, and [arguments] with a command documented in the user guide. Running the image's bare CDesk command would normally select its bundled copy instead of /opt/cdesk-dev/CDesk.

The checked-out source and the image's dependencies must be compatible. This procedure selects the branch code; it does not install or repair software dependencies. Use an image release validated for the intended source revision.

### 4. Configure software and reference data

CDesk reads config.json next to its resolved entry script by default. Use --config to select a different file. The configuration has three sections:

| Section | Purpose |
| --- | --- |
| software | Paths to separately installed tools and the Celescope executable directory. |
| conda_env | Root directories of the four main CDesk environments; do not append /bin. |
| data | Species/reference keys and paths to their genome, annotation, and index resources. |

The [config.json in this branch](./config.json) uses container paths under /opt/conda, /root, and /CDesk/TAIR. Those paths must exist **inside the running container**. They are not host installation paths.

The default reference key is **TAIR10**, which points to resources under **/CDesk/TAIR/**. Use -s TAIR10 for commands that accept the -s species option, or the species option shown by the command's help. This configuration does not define a TAIR10_chr key. For another reference, add its own data entry and use that entry's exact key.

<details>
<summary>Default configuration for this branch</summary>

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

Keep the software entries required by the workflows you intend to run. For example, the Cell Ranger preprocessing implementation invokes software.cellranger directly, even when Cell Ranger is installed inside the container. A celescope_path value is a directory containing the Celescope executables; juicer_tools_jar is the JAR file path.

For native installations, replace all container-specific paths with the actual local installation paths. CDesk searches for environments named CDesk, CDesk_py3.7, CDesk_py2.7, and CDesk_R using Conda, then Mamba, before falling back to paths in conda_env. Thus a detected environment with a matching name takes precedence over its configured path.

Reference fields depend on the selected workflow:

| Field | Resource |
| --- | --- |
| hisat2_index | HISAT2 index prefix, not just its parent directory. |
| refseq_gtf | Gene annotation in GTF format. |
| refseq_bed | Gene annotation in BED format. |
| chromInfo | Chromosome names and lengths. |
| fasta | Reference genome FASTA. |
| tf_file | Transcription factor list. |
| promoter_file | Promoter regions in BED format. |
| TE_idx | scTE index. |
| scRef10x | Cell Ranger reference directory. |
| effective_genome_size | Effective genome size accepted by the relevant workflow. |
| gff3 | Gene annotation in GFF/GFF3 format. |
| refgenes | Gene annotation resource used by the DrSeq workflow. |
| bowtie2_mapindex | Bowtie2 index prefix. |
| singleron_mapindex | Singleron/Celescope reference index directory. |
| dnbc_mapindex | DNBC4Tools reference index directory. |

Use a consistent genome assembly and chromosome naming scheme across FASTA, annotations, indices, BAM, BigWig, BED, and Hi-C resources. A reference key such as TAIR10 is a configuration label; it does not convert chromosome names or make resources from different annotations interchangeable. Cloning this repository does not download the reference datasets or build their indices.

### 5. Verify the installation

Run the checks in the environment where the analysis will execute. Use the explicit branch entry point from section 3 when testing development code.

~~~bash
CDesk -h
CDesk ATAC preprocess -h
CDesk HiC preprocess -h
~~~

Help output confirms that the command-line interface can be parsed; it does **not** validate all dependencies or reference files. In the current entry script, help exits before the runtime environment initialization, although configuration loading and module imports occur before help is displayed.

Before starting a workflow, confirm its configured executables, environment directories, and reference files exist, then run a documented example with suitable input data. Save the selected image identifier, source commit (for development runs), configuration, sample/group tables, command, and log with the analysis results. A successful help command is not evidence that all five analysis modules have passed.

### 6. Advanced: install from source with Conda/Mamba

Use this route when you need to maintain the software environment yourself. The commands below describe the repository layout; they are not a guarantee that the current environment files resolve on every platform or reproduce a released image. Use a compatible Linux system and consult the individual tools' platform requirements.

Clone the branch as in section 3, then run the following from the repository root:

~~~bash
# Run from env/ so the pip_requirements.txt reference in CDesk_py3.7.yml resolves.
cd env
mamba env create -f CDesk.yml
mamba env create -f CDesk_py3.7.yml
mamba env create -f CDesk_py2.7.yml
mamba env create -f CDesk_R.yml
cd ..

chmod +x CDesk
find src -type f -name '*.sh' -exec chmod +x {} \;
~~~

Install tools not supplied by those environments as required by your workflow, including Cell Ranger, DrSeq, DNBC4Tools, Celescope, or the Juicer tools JAR. Obtain and build the relevant reference resources, copy config.json to a local configuration file, and replace its container paths with your actual paths.

You can invoke the entry script without a system-wide installation:

~~~bash
python3 /absolute/path/to/CDesk/CDesk \
  --config /absolute/path/to/config.local.json \
  -h
~~~

Alternatively, expose the executable through your user-local bin directory:

~~~bash
mkdir -p "$HOME/.local/bin"
ln -s /absolute/path/to/CDesk/CDesk "$HOME/.local/bin/CDesk"
export PATH="$HOME/.local/bin:$PATH"

CDesk --config /absolute/path/to/config.local.json -h
~~~

Add the PATH export to your shell configuration if you want it to persist. Complete the workflow-specific verification in section 5 before relying on a native installation.

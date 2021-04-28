<p align="center">
    <img width="50%" src="https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210401150558.png">
</p>

[![Anaconda-Server Badge](https://anaconda.org/bioconda/tritimap/badges/version.svg)](https://anaconda.org/bioconda/tritimap) [![Anaconda-Server Badge](https://anaconda.org/bioconda/tritimap/badges/latest_release_date.svg)](https://anaconda.org/bioconda/tritimap) [![Anaconda-Server Badge](https://anaconda.org/bioconda/tritimap/badges/downloads.svg)](https://anaconda.org/bioconda/tritimap) ![Docker Pulls](https://img.shields.io/docker/pulls/fei0810/tritimap) ![Anaconda-Server Badge](https://anaconda.org/bioconda/tritimap/badges/license.svg)

**Triti-Map** is a Snakemake-based pipeline for gene mapping in Triticeae, which contains a suite of user-friendly computational packages and [web-interface](http://bioinfo.cemps.ac.cn/tritimap/) integrating multi-omics data from Triticeae species including genomic, epigenomic, evolutionary and homologous information.

**Triti-Map** could efficiently explore trait-related genes or functional elements not present in the reference genome and reduce the time and labor required for gene mapping in large genome species.

**More thorough information and explanations are provided in the [Triti-Map Wiki](https://github.com/fei0810/Triti-Map/wiki).**

<h3 align="center">Triti-Map workflow overview</h3>

<p align="center">
    <img width="80%" src="https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210330151508.png">
</p>

Triti-Map ptimization steps to address specific challenges of Triticeae gene-mapping

<p align="center">
    <img width="80%" src="https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210330151509.png">
</p>

## Getting Started with Triti-Map

### Installation

#### Installing from Bioconda

First, to install Triti-Map you need a UNIX environment contains [Bioconda](http://bioconda.github.io/). See [how to install Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

```sh
# create new environment and install Triti-Map
conda create -c conda-forge -c bioconda -n tritimap tritimap
# activate Triti-Map environment
conda activate tritimap
# test Triti-Map
tritimap --help
```

#### Installing from Docker

You can also use [Triti-Map](https://hub.docker.com/r/fei0810/tritimap) via Docker. See [how to install Docker](https://docs.docker.com/engine/install/), then download and run this image using the following commands:

```sh
# docker pull command
docker pull fei0810/tritimap:v0.9.5
# run docker
docker run -i -t fei0810/tritimap:v0.9.5 /bin/bash
```

#### Installing from GitHub

```sh
# download Triti-Map
git clone https://github.com/fei0810/Triti-Map.git
cd Triti-Map
# install Triti-Map
python setup.py install
```

When using source code for installation, you need to install other dependencies of Triti-Map by yourself. You can view Triti-Map dependent software via `tritimap_env.yaml`

### Preparing relevant files

#### Genome and annotation files

Downloading the genome and annotation files you need. Here are some links to download the genome of the Triticeae species.

- [_Triticum aestivum_](http://plants.ensembl.org/Triticum_aestivum/Info/Index)
- [_Triticum turgidum_](http://plants.ensembl.org/Triticum_turgidum/Info/Index)
- [_Triticum dicoccoides_](https://plants.ensembl.org/Triticum_dicoccoides/Info/Index)
- [_Hordeum vulgare_](http://plants.ensembl.org/Hordeum_vulgare/Info/Index)
- [_Aegilops tauschii_](https://plants.ensembl.org/Aegilops_tauschii/Info/Index)
- [_Triticum urartu_](http://www.mbkbase.org/Tu/)

Building GATK and samtools index file

```sh
# for example
# GATK index
gatk CreateSequenceDictionary -R /genome/path/genome.fasta
# samtools index
samtools faidx /genome/path/genome.fasta
```

DNA-seq data use [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)。

```sh
# for example
bwa-mem2 index /genome/path/genome.fasta
```

RNA-seq data use [STAR](https://github.com/alexdobin/STAR)。

```sh
# for example
STAR --runThreadN 30 --runMode genomeGenerate \
--genomeDir /star/index/path/genome_star \
--genomeFastaFiles /genome/path/genome.fasta \
--sjdbOverhang 100 \
--sjdbGTFfile /anntotaion/path/genome.gtf \
--genomeChrBinNbits 18 \
--limitGenomeGenerateRAM 50805727274
```

#### Configuration files

```sh
# generate configuration file in current directory
tritimap init

#Or generate configuration file in running directory
tritimap init -d /your/work/path
```

When `tritimap init` is run successfully, the working directory will generate three configuration files that you need to modify.

- `config.yaml`: Triti-Map configuration file
- `sample.csv`: Sample information file
- `region.csv`: Chromosome region file used to filter the raw results (required only when running `Assembly Module` alone)

**The [Triti-Map wiki](https://github.com/fei0810/Triti-Map/wiki/preparing_config) contains detailed information about the meaning and usage of each parameter in the configuration file.**

### Running Triti-Map

```sh
conda activate tritimap
# running directory
cd /your/work/path

# three types of analysis method

# running both Interval Mapping Module and Assembly Module
tritimap run -j 30 all
# only running Interval Mapping Module
tritimap run -j 30 only_mapping
# only running Assembly Module
tritimap run -j 30 only_assembly

```

Triti-Map supports three types of analysis method.

- **`tritimap run -j 30 only_mapping`**: If you only need to identify trait association intervals and mutations, then run the **Interval Mapping Module**.
- **`tritimap run -j 30 only_assembly`**: If you only need to identify trait association new genes, then run the **Assembly Module**.
- **`tritimap run -j 30 all`**: run the **Interval Mapping Module** and the **Assembly Module** together.

Note: Triti-Map pipeline may take a long time to run(1 to 2 days). The `screen` command is useful for the cases when you need to start a long-running process. Learn more about [GNU Screen](https://www.gnu.org/software/screen/).

## Exploring Triti-Map's results

A complete catalog of Triti-Map results is shown below.

```
├── results
│   ├── 01_cleandata
│   ├── 02_mergedata
│   ├── 03_mappingout
│   ├── 04_GATKout
│   ├── 05_vcfout
│   ├── 06_regionout
│   ├── 07_assembleout
│   ├── logs

```

**You can learn about the results generated by Triti-Map in the [Triti-Map wiki](https://github.com/fei0810/Triti-Map/wiki/results_tritimap).**

## Triti-Map Annotation Platform

[Triti-Map Annotation Platform](http://bioinfo.cemps.ac.cn/tritimap/) is an online analysis module of Triti-Map. To locate causal variants and candidate genes or regulatory elements, Triti-Map integrated multi-omics data and various information from Triticeae species to provide a functional and evolutionary characterization of SNPs, genes, genomic regions, and new sequences related to the target trait.

The platform can perform various analyses, including SNP annotation and visual display, homologous gene analysis, collinearity analysis, and new sequence function annotation, providing richer reference information for Triticeae gene mapping.

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210411112445.png)

## More

### Frequently Asked Questions

**You can also check out some of the [FAQ](https://github.com/fei0810/Triti-Map/wiki/faq_tritimap) you may encounter during use Triti-Map**

### Citing

### Author/Support

Fei Zhao zhaofei920810@gmail.com

Lab Home Page：http://bioinfo.cemps.ac.cn/zhanglab/

Issues can be raised at: https://github.com/fei0810/Triti-Map/issues

We also encourage you to contribute to TriTi-Map! To fix bugs or add new features you need to create a Pull Request.

### Acknowledgements

Thanks to [@xuzhougeng](https://github.com/xuzhougeng) who provided installation support for Bioconda, and thanks to [@zwbao](https://github.com/zwbao) who offered Docker installation support.

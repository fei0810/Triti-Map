## Introduction

<p align="center">
    <img width="50%" src="https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210401150558.png">
</p>

**[Triti-Map](https://github.com/fei0810/Triti-Map)** is a Snakemake-based pipeline for gene mapping in Triticeae, which contains a suite of user-friendly computational packages and web-interface integrating multi-omics data from Triticeae species including genomic, epigenomic, evolutionary and homologous information.

Triti-Map could efficiently explore trait-related genes or functional elements not present in the reference genome and reduce the time and labor required for gene mapping in large genome species.

Triti-Map including three analysis modules:

- Interval Mapping Module
- Assembly Module
- Web-based Annotation Module

**Interval Mapping Module** and **Assembly Module** are command-line software. The input is mixed pool sequencing DNA-Seq(ChIP-Seq/WGS) or RNA-Seq data, which can be generated in one step including trait-related genomic interval, mutation sites and new genes.

Web-based Annotation Module is an [online analysis platform](http://bioinfo.cemps.ac.cn/tritimap/). We collected genomic information from each species of the Triticeae and hexaploid wheat varieties that have been sequenced. Gene function annotation and transcription factor binding site prediction are uniformly performed, while representative histone modification data of each species are integrated. The platform can perform various analyses, including SNP annotation and visualization, homologous gene analysis, Triticeae collinearity analysis, and new sequence function annotation, providing richer reference information for wheat family gene cloning.

<h3 align="center">Triti-Map workflow overview</h3>

<p align="center">
    <img width="80%" src="https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210330151508.png">
</p>

Triti-Map ptimization steps to address specific challenges of Triticeae gene-mapping

<p align="center">
    <img width="80%" src="https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210330151509.png">
</p>

## Author information

Developer: zhaofei920810@gmail.com

Lab web page: [http://bioinfo.cemps.ac.cn/zhanglab/](http://bioinfo.cemps.ac.cn/zhanglab/)

## Installation

### Installing from Bioconda

First, to install Triti-Map you need a UNIX environment contains [Bioconda](http://bioconda.github.io/). See [how to install Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

```sh
# create new environment and install Triti-Map
conda create -c conda-forge -c bioconda -n tritimap tritimap
# activate Triti-Map environment
conda activate tritimap
# test Triti-Map
tritimap --help
```

### Installing from Docker

You can also use [Triti-Map](https://hub.docker.com/r/fei0810/tritimap) via Docker. See [how to install Docker](https://docs.docker.com/engine/install/), then download and run this image using the following commands:

```sh
# docker pull command
docker pull fei0810/tritimap:v0.9.4
# run docker
docker run -i -t fei0810/tritimap:v0.9.4 /bin/bash
```

### Installing from GitHub

```sh
# download Triti-Map
git clone https://github.com/fei0810/Triti-Map.git
cd Triti-Map
# install Triti-Map
python setup.py install
```

When using source code for installation, you need to install other dependencies of Triti-Map by yourself. You can view Triti-Map dependent software via `tritimap_env.yaml`

After correct installation you can see a help message similar to the following.

```sh
$ tritimap -h
Usage: tritimap [OPTIONS] COMMAND [ARGS]...

  Triti-Map: A Snakemake-based pipeline for gene mapping in Triticeae. For
  more information, see: https://github.com/fei0810/Triti-Map

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  init  Generate snakemake configuration file and other needed file.
  run   Run Triti-Map pipeline
```

Note: running Triti-Map in an independent environment can avoid the influence of other software.

## Prepare relevant documents

### Genome and annotation files

Triti-Map mapping DNA-seq(ChIP-Seq and WGS data) to genome file with bwa-mem2, and mapping RNA-Seq data to genome file with STAR, Before formal analysis, you need to prepare relevant genomes and annotation files, and build index files.

The genome and annotation files of some sequenced species of the wheat race can be downloaded through the following link.

- [_Triticum aestivum_](http://plants.ensembl.org/Triticum_aestivum/Info/Index)
- [_Triticum turgidum_](http://plants.ensembl.org/Triticum_turgidum/Info/Index)
- [_Triticum dicoccoides_](https://plants.ensembl.org/Triticum_dicoccoides/Info/Index)
- [_Hordeum vulgare_](http://plants.ensembl.org/Hordeum_vulgare/Info/Index)
- [_Aegilops tauschii_](https://plants.ensembl.org/Aegilops_tauschii/Info/Index)
- [_Triticum urartu_](http://www.mbkbase.org/Tu/)

**Building a genome index file**

```
conda activate tritimap
# bulid gatk genome index
gatk CreateSequenceDictionary -R /genome/path/genome.fasta
# build samtools genome index
samtools faidx /genome/path/genome.fasta
```

**Building an index file for the alignment software**

DNA-seq data use [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) .

```
bwa-mem2 index /genome/path/genome.fasta
```

RNA-seq data use [STAR](https://github.com/alexdobin/STAR) .

For example:

```
STAR --runThreadN 30 --runMode genomeGenerate \
--genomeDir /star/index/path/genome_star \
--genomeFastaFiles /genome/path/genome.fasta \
--sjdbOverhang 100 \
--sjdbGTFfile /anntotaion/path/genome.gtf \
--genomeChrBinNbits 18 \
--limitGenomeGenerateRAM 50805727274
```

When using STAR to construct a hexaploid wheat genome index, parameters are required `--limitGenomeGenerateRAM` to increase available memory. The value varies slightly from genome to genome.

When running Triti-Map, Triti-Map automatically detects if there is a correct matching index file based on the genome file path.

### Configuration files

Use `tritimap init` generate configuration files required for running the main program.

```
conda activate tritimap
cd /your/work/path
tritimap init
```

Three configuration files will be generated in the working directory after the command is run successfully.

- **config.yaml**:Triti-Map main program operating parameters.
- **sample.csv**: Raw analysis data sample information.
- **region.csv**: Sequence assembly filter interval file. (Required only when using `Assembly Module` alone)

### Modifying the sample information file `cample.csv`

The sample information is a comma-separated CSV file that includes the sample name, sequencing type, bulk genotype, and file location of the data.

Triti-Map can accept 2 sets of input data with only pools or 4 sets of input data containing both pools and parents. To improve mutation identification and assembly, it is recommended to use ChIP-seq data with 2-3 different histone modifications, such as H3K27me3, H3K4me3 and H3K36me3.

**Description**

- `sample`: the name of the sample, each sample name should be unique.
- `bulk`: the name of the bulk, the names of different samples belonging to the same bulk should be the same; if the parents are included, the names of the parents and the corresponding bulk should also be the same.
- `type`: sample data type, this column can be filled with either **pool** (corresponding to bulk data) or **parent** (corresponding to parent data).
- `genotype`: genotype of the trait-related genes, either **recessive** (corresponding to recessive) or **dominant** (corresponding to dominant), or **unknown** if you not sure.
- `marker`: sequencing type, such as h3k36me3, h3k27me3, RNAseq and WGS, etc. This column is used to distinguish different samples from the same bulk.
- `fq1` and `fq2`: the specific paths of sequencing files R1 and R2, **absolute paths need to be used**.

**Example 1**

A typical multiple-modification ChIP-seq sample information file sample.csv is shown below. The input file is 2 sets of bulk data without parents, where each bulk contains ChIP-seq data with 3 different histone modifications.

```
sample,bulk,type,genotype,marker,fq1,fq2
res1,Res,pool,dominant,h3k27me3,/data/path/Res_H3K27me3_R1.fq.gz,/data/path/Res_H3K27me3_R2.fq.gz
res2,Res,pool,dominant,h3k36me3,/data/path/Res_H3K36me3_R1.fq.gz,/data/path/Res_H3K36me3_R2.fq.gz
res3,Res,pool,dominant,h3k4me3,/data/path/Res_H3K4me3_R1.fq.gz,/data/path/Res_H3K4me3_R2.fq.gz
sus1,Sus,pool,recessive,h3k27me3,/data/path/Sus_H3K27me3_R1.fq.gz,/data/path/Sus_H3K27me3_R2.fq.gz
sus2,Sus,pool,recessive,h3k36me3,/data/path/Sus_H3K36me3_R1.fq.gz,/data/path/Sus_H3K36me3_R2.fq.gz
sus3,Sus,pool,recessive,h3k4me3,/data/path/Sus_H3K4me3_R1.fq.gz,/data/path/Sus_H3K4me3_R2.fq.gz
```

**Example 2**

WGS bulk data with parents:

```
sample,bulk,type,genotype,marker,fq1,fq2
res1,Res,pool,dominant,wgs,/data/path/Res_wgs_R1.fq.gz,/data/path/Res_wgs_R2.fq.gz
sus1,Sus,pool,recessive,wgs,/data/path/Sus_wgs_R1.fq.gz,/data/path/Sus_wgs_R2.fq.gz
resp,Res,parent,dominant,wgs,/data/path/data/p_Res_wgs_R1.fq.gz,/data/path/data/p_Res-wgs_R2.fq.gz
susp,Sus,parent,recessive,wgs,/data/path/data/p_Sus_wgs_R1.fq.gz,/data/path/data/p_Sus_wgs_R2.fq.gz
```

**Example 3**

RNA-seq pool data without parent:

```
sample,bulk,type,genotype,marker,fq1,fq2
res1,Res,pool,dominant,rnaseq,/data/path/Res_rnaseq_R1.fq.gz,/data/path/Res_rnaseq_R2.fq.gz
sus1,Sus,pool,recessive,rnaseq,/data/path/Sus_rnaseq_R1.fq.gz,/data/path/Sus_rnaseq_R2.fq.gz
```

### Modifying the `region.csv` file

**Note: This file needs to be modified only when using Assembly Module alone.**

When using Triti-Map for de novo sequence assembly, you need to use the trait association interval to filter the assembled sequences, which will be generated automatically if you run the full analysis modules; if you only run the Assembly Module, you need to fill in the filter interval in the region.csv configuration file in advance.

A typical filter interval information file region.csv is shown below.

```
chrom,start,end
chr7A,724000000,730000000
```

- `chrom`: chromosome where the interval is positioned
- `start`: position of the start of the interval
- `end`: position of the end of the interval

### Modifying the configuration file `config.yaml`

To run the Triti-Map main program, you need to modify the parameters in the configuration file.

**Description of the main parameters of the configuration file**

- `email`: **Important** , you need to provide your personal email when using EMBL-EBI API for related analysis.
- `samples`: No modification is needed, the path of the sample information file. The default sample information file is `sample.csv` in the Triti-Map running directory.
- `datatype`: **Important** , the sequencing type of the sample data, **dna**: ChIP-seq or WGS; **rna**: RNA-seq.
- `maxthreads`: the maximum number of threads that can be used. Default value: `30` .
- `ref`: **Important** , reference genome related file paths, contains 3 sub-parameters.

  - `genome`: path to the reference genome file, use absolute path.
  - `annotation`: path to the reference genome gene annotation file, use the absolute path.
  - `STARdir`: STAR reference genome index directory (required for RNA-seq analysis).

- `gatk`: GATK-related analysis parameters, including one sub-parameter.

  - `min_SNP_DP` No modification is needed, the minimum depth that needs to be met for each pool of valid SNPs. default is `10`.

- `snpindex`: Parameters required for BSA localization analysis, including 6 sub-parameters.

  - `pop_struc`: **Important** , the population structure of the pool samples. If the data of the pool is F2 generation, fill in **F2**; if the data of the pool is RIL population or all individuals are homozygous, fill in **RIL**.
  - `bulksize`: The number of samples in the pool, e.g., if each pool consists of 30 samples, then fill in `30`.
  - `winsize`: No modification is needed, the length of the sliding window for data correction, default is `1000000`(1Mb).
  - `filter_percentage`: **Important**, the percentage of the original results to filter based on Delta SNPindex and SNPconut/Mb. If the value is 0.75, it means that the average Delta SNPindex and the average number of SNPs per 1Mb of the candidate interval should be both greater than 75% of the corresponding values of all original results.
  - `fisher_p`: No modification is needed, filter the SNP loci of the trait association interval and calculate the pvalue of each locus using fisher test, default is `0.0001`.
  - `min_length`: No modification is needed, the minimum length of the candidate trait association interval. For large genome species such as wheat, the default length is `1000000`(1Mb).

- `merge_lib`: How to handle multiple sets of different ChIP-seq data of the same pool. **The default is `merge`**, i.e. samples are merged first and then assembled to get better results; if dealing with large genomic data such as hexaploid wheat and your server memory is less than 300G, you can modify it to `split`, i.e. each group of data is assembled separately and then merged for subsequent analysis.
- `memory`: the maximum memory available when assembling transcriptome sequences using SPAdes, `300` means 300G.
- `denovo_filter_method`: **Important** , the pool-specific sequence filtering method, valid when running `only_assembly` module. Set to `external_region` means user customize filter interval in `region.csv` file; set to `external_fasta` means user use own prepared external fasta sequence as filter database, please refer to Q&A.
- `filter_region_file`: If you use `external_region`, you need to fill in the location of the filter interval file here, the default is `region.csv`.
- `filter_fasta_file`: If you use `external_fasta`, you need to fill in the location of the fasta file, e.g. `/your/path/region.fasta`.
- `blast_database`: No modification is needed, the database used for BLAST of assembly sequences. `em_cds_pln` means using EBI ENA plant coding sequence database, `em_std_pln` means using EBI ENA plant standard sequence database. Default value: `em_cds_pln`.

## Running Triti-Mapping

After successfully running `tritimap init` and modifying each configuration file, you can run the main Triti-Map program.

```sh
conda activate tritimap
cd /your/work/path
tritimap run -j 30 all
```

`tritimap run` supports a total of three analysis modes, corresponding to the following parameters and commands.

- **`only_mapping`**
  Running only the **Interval Mapping Module** for BSA analysis, which is used to identify trait association intervals and mutations.
  `tritimap run -j 30 only_mapping`
- **`only_assembly`**
  Running **Assembly Module** for sequence assembly only, used to identify new genes for trait associations.
  `tritimap run -j 30 only_assembly`
- **`all`**
  Running both the **Interval Mapping Module** and the **Assembly Module**.
  `tritimap run -j 30 all`

### Testing run and creating a DAG (directed acyclic graph) of jobs

Before the actual run, you can view all the proposed jobs via flag `-n`.

```
tritimap run -j 30 -n all
```

A message similar to the following can be seen.

```
Job counts:
	count	jobs
	2	AssembleMerge
	2	BlastAnnoEBI
	1	GATK4MergeVcf
	1	MergeCleanPoolFastq
	2	PfamAnnoEBI
	1	QTLseqr_plot
	2	Rawfastq2Cleanfastq
	2	STAR1Mapping
	2	STAR2Mapping
	1	all
	22	gatk4HC
	2	rnaFilter2Uniqmap
	2	rnaGATK4FixMateInfo
	2	rnaGATK4Markdup
	2	rnaGATK4ReplaceRG
	2	rnaGATK4SplitNCigar
	1	uniqScaffold
	2	unmap_BlastAnnoEBI
	2	unmap_PfamAnnoEBI
	1	vcf2tab
	1	vcfGenotypeFiltration
	1	vcfHardFiltration
	56
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

You can create a visualization of the DAG using flag `--dag` and the `dot` command provided by [Graphviz](https://www.graphviz.org/).

```sh
tritimap run -j 30 -n all --dag | dot -Tpdf > tritimap.dag.pdf
```

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210401162907.png)

`tritimap run --help` show help message:

```sh
$ tritimap run -h
Usage: tritimap run [OPTIONS] [[all|only_mapping|only_assembly]]
                    [SNAKEMAKE_ARGS]...

  Triti-Map main command. The pipeline supports three execute modules: all,
  only_mapping and only_assembly. First, you need to fill in the
  configuration file correctly.

Options:
  -d, --working-dir PATH  Triti-Map running directory.
  -c, --config-file PATH  Triti-Map config file, generated with tritimap init.
  -j, --jobs INTEGER      Use at most N CPU cores/jobs in parallel.
  -n, --dryrun            Do not execute anything, and display what would be
                          done.  [default: False]

  --profile TEXT          Name of profile to use for configuring Snakemake.
  -h, --help              Show this message and exit.
```

Because a complete analysis can take a long time, we recommend that you can run Triti-Map to the background using `screen`（https://www.gnu.org/software/screen/) program.

### Viewing and handling error messages

If an error occurs during the analysis process, you can locate the error step through the error message, or you can find the cause of the error through the corresponding log file.

If an error occurs in the middle of the analysis process, you can directly run `tritimap run` again after modifying the correct parameters, and the steps that have already generated results will be automatically skipped and run directly from the unfinished steps.

## Triti-Map results

### Directory structure

The complete Triti-Map result directory is as follows:

```sh
├── results
│   ├── 01_cleandata
│   ├── 02_mergedata
│   ├── 03_mappingout
│   ├── 04_GATKout
│   ├── 05_vcfout
│   ├── 06_regionout
│   ├── 07_assembleout
│   ├── logs
```

The contents within each result directory are as follows.

- `01_cleandata` Clean data after pre-processing of the raw data
- `02_mergedata` The merged pool of data with different histone modification types for subsequent analysis
- `03_mappingout` BAM files after reference genome alignment
- `04_GATKout` BAM file for mutation identification after preprocessing
- `05_vcfout` Original VCF file with mutation information and filtered VCF file
- `06_regionout` SNP input and result files for BSA localization analysis
- `07_assembleout` Sequence assembly generated by pool-specific sequences and corresponding functional annotation files
- `logs` Log files generated at different steps of the analysis process

Temporary files that may be of value during analysis are saved in the `temp_output` sub-directory of each directory.

Temporary files do not need to be noticed under normal circumstances, you can check them in combination with log files when the program has errors.

### Interval Mapping Module Results

The Interval Mapping Module results are located in the directory `06_regionout`, the `xxx` and `yyy` in the following file names indicate the pool names in the `sample.csv` file.

```
06_regionout
├── xxx_yyy_snpindex_input.txt
├── xxx_yyy_qtlseqr_output.txt
├── xxx_yyy_qtlseqr_raw_region.txt
├── xxx_yyy_qtlseqr_filter_region.txt
├── xxx_yyy_qtlseqr_filter_snpinfo.txt
├── xxx_yyy_qtlseqr_filter_indelinfo.txt
├── xxx_yyy_qtlseqr_filter_snp.bed
├── xxx_yyy_qtlseqr_filter_indel.bed
├── xxx_yyy_qtlseqr_SNPcounts_line.pdf
├── xxx_yyy_qtlseqr_SNPcounts_point.pdf
├── xxx_yyy_qtlseqr_SNPindex_line.pdf
├── xxx_yyy_qtlseqr_SNPindex_point.pdf
└── xxx_pool_vs_yyy_pool_candidateregion_1.pdf
```

- `xxx_yyy_snpindex_input.txt` SNP input file for Interval Mapping analysis
- `xxx_yyy_qtlseqr_output.txt` QTLseqr SNP result information obtained by Delta SNPindex analysis
- `xxx_yyy_qtlseqr_raw_region.txt` QTLseqr Raw interval result information from Delta SNPindex analysis
- `xxx_yyy_qtlseqr_filter_region.txt` High confidence trait-related interval information by filtering
- `xxx_yyy_qtlseqr_filter_snpinfo.txt` SNP information in the high confidence interval
- `xxx_yyy_qtlseqr_filter_indelinfo.txt` INDEL information in the high confidence interval
- `xxx_yyy_qtlseqr_filter_snp.bed` SNP bed format file for the high confidence interval, **which can be annotated in Triti-Map online analysis platform**
- ` xxx_yyy_qtlseqr_filter_indel.bed` INDEL bed format file in the high confidence interval
- `*.pdf` pdf file for visualization of the results, including all chromosome Delta SNPindex, snp counts plots, and Delta SNPindex plots for each high confidence interval

Generally speaking, the results that are most worthy of further analysis are `xxx_yyy_qtlseqr_filter_region.txt`, `xxx_yyy_qtlseqr_filter_snpinfo.txt` and `xxx_yyy_qtlseqr_filter_indelinfo.txt`.

#### Results Format Description

**xx_yyy_snpindex_input.txt**

```
#CHROM  POSm1   POS     REF     ALT     Res_pool_ref    Res_pool_alt    Res_pool_depth  Res_pool_ratio  Sus_pool_ref    Sus_pool_alt    Sus_pool_depth  Sus_pool_ratiosnpindex
chr1A   1144541 1144542 G       A       13      147     160     0.91875 1       144     145     0.993103        -0.0743534
```

Left to right: chromosome, SNP postion-1, SNP postion, REF , ALT, pool A Ref number, pool A ALT number, SNPindex of pool A, pool B Ref number, pool B ALT number, SNPindex of pool B, Delta SNPindex.

**xxx_yyy_qtlseqr_filter_region.txt**

```
CHROM   qtl     start   end     length  nSNPs   avgSNPs_Mb      peakDeltaSNP    posPeakDeltaSNP avgDeltaSNP
chr7A   7       724111912       730119678       6007766 3262    543     0.954374003095134       729506093       0.812296971085005
```

Left to right: interval chromosome, qtlseqr raw results interval number, interval start position, interval end position, interval length; number of SNPs in the interval, average number of SNPs per 1Mbp; Delta SNPindex peak; Delta SNPindex peak location; average Delta SNPindex. For more information, please refer to [QTLseqr](https://github.com/bmansfeld/QTLseqR/).

**xxx_yyy_qtlseqr_filter_snpinfo.txt**

```
CHROM   POS     REF     ALT     AD_REF.Sus_pool AD_ALT.Sus_pool DP.Sus_pool     SNPindex.Sus_pool       AD_REF.Res_pool AD_ALT.Res_pool DP.Res_pool     SNPindex.Res_pool      REF_FRQ DeltaSNP        nSNPs   tricubeDeltaSNP minDP   tricubeDP       CI_99   fisher_p
chr7A   724111912       T       G       0       12      12      1       21      0       21      0       0.636363636363636       -1      445     0.505705636608414     12       20      -0.5    2.8183517084228e-09`
```

Left to right: chromosome, SNP postion, REF , ALT, pool A Ref number, pool A ALT number, pool A depth, SNPindex of pool A, pool B Ref number, pool B ALT number, SNPindex of pool B, qtlseqr total reference allele frequnce, Delta SNPindex, snp count, tricube, Delta SNPindex, min depth, trichbe depth, 99% confidence interval Delta SNPindex, fisher test pvalue

#### figure

**All chromosome Delta SNPindex plot**

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210310150108.png)

X: Chromosome, Y:Delta SNPindex. The blue line indicates the 99% confidence interval calculated by QTLseqR corresponding to the Delta SNPindex value. The gray dots indicate that the corresponding SNP is below this value, and the black dots indicate that the corresponding SNP is above this value.

**All chromosome SNP count plot**

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210310150210.png)

X: Chromosome, Y:SNP count.

**Candidate interval Delta SNPindex plot**

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210401170136.png)

The blue line indicates the 99% confidence interval calculated by QTLseqR corresponding to the Delta SNPindex value. The shaded area indicates the trait-related interval.

### Assembly Module Results

The Assembly Module results are located in the directory `07_assembleout`, the `xxx` in the following file names indicate the pool names in the `sample.csv` file.

```sh
07_assembleout
├── xxx_merge_denovo_scaffolds.fasta
├── xxx_candiadate_denovo.fasta
├── xxx_candidate_denovo2ref.info.txt
├── xxx_candiadate_denovo_pfam_anno.fasta
├── xxx_candiadate_denovo_pfam_anno.txt
├── xxx_candiadate_denovo_blast_anno.txt
├── xxx_unmap_denovo.fasta
├── xxx_unmap_denovo_pfam_anno.fasta
├── xxx_unmap_denovo_pfam_anno.txt
└── xxx_unmap_denovo_blast_anno.txt
```

- `*_merge_denovo_scaffolds.fasta`: pool-specific sequences that do not match the reference genome properly
- `*_candiadate_denovo.fasta`: pool-specific sequences that cannot be matched normally with the reference genome and partially matched with the trait-related region
- `*_candidate_denovo2ref.info.txt`: Information on the position of `*_candiadate_denovo.fasta` sequences relative to the reference genome
- `*_candiadate_denovo_pfam_anno.txt`: `*_candiadate_denovo.fasta` sequence function annotation results
- `*_candiadate_denovo_pfam_anno.fasta`: the functional sequence in `*_candiadate_denovo.fasta`.
- `*_candiadate_denovo_blast_anno.txt`: similarity annotation information of the functional sequences
- `*_unmap_*`: indicates the sequence information that does not match with the reference genome at all, corresponding to `*_candidater_*`

Generally speaking, the results that are most worthy of further analysis are `*_candiadate_denovo_pfam_anno.txt`, `*_candiadate_denovo_blast_anno.txt` and `xxx_candiadate_denovo_pfam_anno.fasta`.

**\*\_candiadate_denovo.fasta**

```
>17655292:chr1A:2501389:378S202M:580:Res
TTTAGTTCTAGCCACATAATTATTTCCAGGATTTTATTAAATGATTTAGTATTTTTACTA
AACCAAAACAAAACAGAACAATAAAAAAAAACTGAAATAGAAAGCAAATAGAAATAAGAG
AGAGAGAGAGAGAGAGAGAGTAGGCCACTAGGCCACTTACCTGGACCTACCTGGCCCAGC
GCCAGGCCACCAACGGCCCAGCCCACCAGCACCTCGCCCGTCGTCTTCCTCCCCTCGCCC
CGAAGCAGCTGCGTGCCAGCAGCAGCGCACACCGGAGGTCTCCACCTCCTGCTTCGCGGT
GGCAGCCTCCCTGGCCTCCTAACCGCGCCACGGAGACACCCGCGGCCCCCTCGACCCTCT
CTGGCTCTCCCTGGACCTCTCTCCCTCCTCTGCTCTCTCCCTACGACGCCACCGAACACG
CCCATCGCCGCCGTTGGCCGTAGCCGCGACCACCGTACCTTCCTAGCCCATCCAACATGC
TCGAGAATTCCGCCTCGACCCCCTCTTTCTTCCACCCCAAGCCACGCCACCGCGGGAGCC
CCACATCGCCGCCACGAGTCGTCTTCCCCGAGCACGGCCG
```

**examle ID**:`17655292:chr1A:2501389:378S202M:580:Res`
**explanation** `original assembly sequence ID:located chromosome:chromosome position:CIGAR information:sequence length:pool name`

**\*\_candiadate_denovo_pfam_anno.txt**

```
seqid   Description    E-value       Model
23558405:chr7B:733811080:14S493M9I1522M:2038:res        NB-ARC NB-ARC domain    2.3e-37       PF00931.22
```

Left to right: sequence ID, function description, E-value, pfam domain ID

**\*\_candiadate_denovo_blast_anno.txt**

```
seqid	hit_db	hit_id	hit_desc	hit_url	hsp_bit_score	hsp_align_len	hsp_identity	hsp_query_from	hsp_query_to	hsp_hit_from	hsp_hit_to
23558405:chr7B:733811080:14S493M9I1522M:2038:res	EM_CDS	AUO29721.1	Triticum urartu powdery mildew resistance protein	https://www.ebi.ac.uk/ena/data/view/AUO29721.1	3514.25	1958	99.8	81	2038	1	1958
```

Left to right: sequence ID, blast database, hit sequence id, hit sequence description, hit sequence ENA url, blast bit score, align length, blast identity, query start position, query end position, hit start position, hit end position

## Triti-Map Annotation Platform

[Triti-Map Annotation Platform](http://bioinfo.cemps.ac.cn/tritimap/) is an online analysis module of Triti-Map.

To locate causal variants and candidate genes or regulatory elements, Triti-Map integrated multi-omics data and various information from Triticeae species to provide a functional and evolutionary characterization of SNPs, genes, genomic regions, and new sequences related to the target trait.

The platform can perform various analyses, including SNP annotation and visual display, homologous gene analysis, collinearity analysis, and new sequence function annotation, providing richer reference information for Triticeae gene mapping. Detailed instructions can be found in the manual.

Functional diagram of Triti-Map Web-based Annotation Platform

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210411112445.png)

## Support

### FAQ

**Q: Which species can Triti-Map be used to analyze?**

A: The Triti-Map analysis pipeline is designed for Triticeae species represented by hexaploid wheat, with several steps optimized for large genomic crops and related adjustments to the way results are visualized. However, in principle, Triti-Map can perform BSA analysis and sequence assembly analysis of trait-related genes on pool sequencing data of any species.

It should be noted that Triti-Map is not suitable for processing genomic data with a large number of scaffolds. If the genome file has a large number of scaffolds in addition to the regular chromosomes, it is recommended to use `samtools faidx` command or `seqkit grep` command of seqkit to first extract the chromosome-level genome as the reference genome file before analysis.

**Q: What kind of data can Triti-Map analyze?**

A: In terms of high-throughput sequencing data types, Triti-Map is designed and developed for ChIP-seq data, but in principle, the Interval Mapping Module can analyze most high-throughput sequencing data including whole genome resequencing (WGS), exome sequencing and transcriptome sequencing.

It is important to note that exome sequencing data is not applicable to Assembly Module; using WGS data of large genomic species for Assembly Module analysis requires sufficient computational resources and is time-consuming.

As far as the type of sequenced material is concerned, the Interval Mapping Module of Triti-Map currently only supports using the QTLseq Delta snpindex method and does not set corresponding filters for EMS mutagenesis data, so it is currently more suitable for analyzing materials with specific traits obtained through hybrid recombination. The Assembly Module is particularly suitable for cases such as crosses of wheat with different ploidy to obtain superior genes, which often have large segments of variation compared to the reference genome or do not exist in the reference genome.

**Q: Which histone modifications should be used by Triti-Map when using ChIP-seq data?**

A: It has been tested that the combination of H3K27me3 and H3K4me3 is more suitable for performing the Interval Mapping Module analysis, and H3K36me3 is suitable for the Assembly Module analysis. If both modules are performed, it is recommended to use the above three modifications simultaneously.

**Q: When using only the Assembly Module of Triti-Map, how can I filter the assembled sequences using the specified intervals?**

A: Some plant materials may have known approximate locus intervals but the trait-related genes cannot be found in the genome, so you can use the Assembly Module of Triti-Map only for new gene identification.

If the known trait locus interval length is too small (less than 1Mb), it is recommended to modify the start and stop sites appropriately to increase the interval length, which is recommended to be no less than 5Mb, and then fill the interval in the `region.csv` file generated by `tritimap init`.

In the configuration file `config.yaml`, set the `denovo_filter_method` parameter to `external_region` and `filter_region_file` to the file path of `region.csv`.

**Q: When using only the Assembly Module of Triti-Map, how can I filter the assembly sequences using the specified database?**

A: In view of the complex homology and variability of subgenomes among Triticeae, Triti-Map specially introduces external sequence database parameters for assembled sequence filtering when only the Assembly Module is used for analysis.

After high confidence trait-related intervals have been obtained by Triti-Map's Interval Mapping Module or other methods, you can download the collinearity sequences of this interval in all Triticeae subgenomes by Triti-Map [Web-based Annotation Module](http://bioinfo.cemps.ac.cn/tritimap).

In the configuration file `config.yaml`, set the parameter `denovo_filter_method` to `external_fasta`, and set `filter_region_file` to the path of the downloaded fasta file, such as `region.fasta`.

In this mode, the assembled pool-specific sequences will be compared with all collinearity regions, thus screening for new genes that may be associated with the trait in a larger scale.

**Q: How should the parameter `filter_percentage` be set when Triti-Map filters the original trait-related intervals?**

A: There is no fixed setting standard for `filter_percentage`, but it is related to the accuracy of sample sampling in the pool, the sequencing depth of the pool and the purity of the sample material.

For two mixed pools with pure material and sufficient sequencing depth, the `filter_percentage` can be set to `0.75` or above; for non-pure material or low sequencing depth, the `filter_percentage` can be set to `0.5`.

If there is no high confidence interval due to unreasonable setting of `filter_percentage`, Triti-Map will report an error automatically. Then you can reset the appropriate `filter_percentage` according to the original result in `xxx_yyy_qtlseqr_raw_region.txt`, and just re-run the main program.

**Q: What performance servers are required to run Triti-Map?**

Past analysis experience shows that the most memory-consuming part of the Triti-Map analysis process is the assembly step in the Assembly Module.

Taking the 400M reads pool ChIP-seq data of hexaploid wheat as an example, it is recommended to have at least 400G of server memory for the Assembly Module analysis.

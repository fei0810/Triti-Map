## Triti-Map 介绍

Triti-Map 是借助 Snakemake 流程语言开发的一套基于多组学 BSA 定位分析同时结合 _de novo_ 序列组装的性状关联基因定位流程。Triti-Map 包括 Interval Mapping Module
、Assembly Module 和 Web-based Annotation Module 三个分析模块。该分析流程针对小麦族大基因组物种开发，可用来进行常规 BSA 定位分析并挖掘参考基因组中不存在的新功能基因。

其中 Interval Mapping Module 和 Assembly Module 为命令行软件，输入数据为混池测序的 DNA-Seq 数据（ChIP-Seq/WGS）或 RNA-Seq 数据，可以一步生成包括性状关联区间、性状关联突变位点和性状关联新基因在内的多种结果。

Web-based Annotation Module 是在线分析平台。我们收集了小麦族各物种和六倍体小麦已有测序品种的基因组信息，对小麦族各物种统一进行了基因功能注释和转录组因子结合位点预测，同时整合了各物种有代表性的表观修饰数据，搭建了配套 Triti-Map 分析结果使用的在线平台。该平台可以进行包括 SNP 注释与可视化分析、同源基因分析、小麦族共线性区间分析和新序列功能注释等在内的多种分析，方便用户对已有定位结果进一步深入分析，为基因克隆提供更加丰富的信息。

**Triti-Map 分析流程示意图**

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210330151508.png)

**流程模块与主要使用软件示意图**

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210330151509.png)

## 作者信息

开发与维护 zhaofei920810@gmail.com
实验室主页 http://bioinfo.sibs.ac.cn/zhanglab/

## 软件安装

安装 conda，推荐安装[miniconda](https://conda.io/en/latest/miniconda.html)即可。

```sh
wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

conda 安装过程中如果有问题，可以参考[miniconda](https://conda.io/en/latest/miniconda.html)和[bioconda](https://bioconda.github.io/user/index.html)的官方详细说明。

conda 安装完成后，新建 Triti-Map 运行环境并通过 bioconda 安装 Triti-Map。

```sh
# 新建运行环境
conda env create -n tritimap
# 启动运行环境
conda activate tritimap
# 安装 Triti-Map
conda install -c bioconda tritimap
```

此外，你也可以通过 github 下载 Triti-Map，自行安装 Triti-Map 和其它依赖软件。

```sh
# 下载Triti-Map分析流程
git clone xxxxxxx
# 进入Triti-Map运行目录
cd Triti-Map
# 安装Triti-Map
python setup.py install
# 安装依赖软件
conda env create -f Triti-Map_env.yml
```

## 准备基因组和注释文件

Triti-Map 针对 DNA-seq（ChIP-Seq 和 WGS 数据）测序数据使用 bwa-mem2 进行比对，针对 RNA-Seq 测序数据使用 STAR 进行比对，运行 Triti-Map 之前首先需要准备对应的基因组文件和 gff3 注释文件。

### 下载基因组和注释文件

小麦族部分以测序物种基因组和注释文件可通过以下链接下载

- [中国春 Chinese Spring](http://plants.ensembl.org/Triticum_aestivum/Info/Index)
- [硬粒小麦 Triticum turgidum durum](http://plants.ensembl.org/Triticum_turgidum/Info/Index)
- [野生二粒小麦 wild emmer wheat Zavitan](https://plants.ensembl.org/Triticum_dicoccoides/Info/Index)
- [大麦 Hordeum vulgare (barley) ](http://plants.ensembl.org/Hordeum_vulgare/Info/Index)
- [粗山羊草 Aegilops tauschii](https://plants.ensembl.org/Aegilops_tauschii/Info/Index)
- [乌拉尔图 Triticum urartu](http://www.mbkbase.org/Tu/)

### 生成基因组 index 文件

```sh
# 启动Triti-Map conda 运行环境
conda activate tritimap
# 建立参考基因组 gatk 索引
gatk CreateSequenceDictionary -R /genome/path/161010_Chinese_Spring_v1.0_pseudomolecules.fasta
# 建立参考基因组 samtools 索引
samtools faidx /genome/path/161010_Chinese_Spring_v1.0_pseudomolecules.fasta

```

### 生成基因组比对 index 文件

DNA-seq 数据使用 bwa-mem2：

注：详细信息可以参考 [bwa-mem2 官方介绍](https://github.com/bwa-mem2/bwa-mem2)。

```sh
# 启动Triti-Map conda 运行环境
conda activate tritimap

# 使用bwa-mem2 建立参考基因组比对index文件
bwa-mem2 index /genome/path/161010_Chinese_Spring_v1.0_pseudomolecules.fasta
```

RNA-seq 数据使用 STAR：

注：详细信息可以参考[STAR 官方文档](https://github.com/alexdobin/STAR)。

```sh
# 启动Triti-Map conda 运行环境
conda activate tritimap
# 使用 STAR 建立六倍体小麦基因组索引需要通过参数 --limitGenomeGenerateRAM 提高默认可用内存

STAR --runThreadN 30 --runMode genomeGenerate \
--genomeDir /star/index/path/wheat_STAR_wholegenome1.1 \
--genomeFastaFiles /genome/path/161010_Chinese_Spring_v1.0_pseudomolecules.fasta \
--sjdbOverhang 100 \
--sjdbGTFfile /anntotaion/path/IWGSC_v1.1_HC_20170706.gtf \
--genomeChrBinNbits 18 \
--limitGenomeGenerateRAM 50805727274
```

说明：Triti-Map 分析流程会根据填写的基因组文件路径自动检测是否有争取匹配的索引文件，如果检测失败会显示对应的生成命令。

## 生成配置文件

```sh
# 启动Triti-Map conda 运行环境
conda activate tritimap

# 进入分析目录
cd /your/work/path

# 生成配置文件
tritimap init
```

`tritimap init` 运行成功后所在目录会生成三个文件：

- sample.csv：原始分析数据样本信息文件
- config.yaml：Triti-Map 运行参数配置文件
- region.csv：过滤区间文件，仅使用`Assembly Module`时需要

#### 填写样本信息文件 cample.csv

样本信息文件为`,`分割的 CSV 格式文件，该表格包括数据的样本名、测序类型、混池基因型和文件位置等信息。

Triti-Map 可以接受仅有混池的 2 组输入数据，或者同时包含混池和亲本的 4 组输入数据。为了提高突变鉴定效果和序列组装效果，使用 ChIP-seq 数据分析时推荐同时使用 2-3 种不同组蛋白修饰数据，如 H3K27me3、H3K4me3 和 H3K36me3。

一个典型的多种修饰 ChIP-seq 样本信息文件`sample.csv`如下所示

```
sample,bulk,type,genotype,marker,fq1,fq2
res1,Res,pool,dominant,h3k27me3,/data/path/Res_H3K27me3_R1.fq.gz,/data/path/Res_H3K27me3_R2.fq.gz
res2,Res,pool,dominant,h3k36me3,/data/path/Res_H3K36me3_R1.fq.gz,/data/path/Res_H3K36me3_R2.fq.gz
res3,Res,pool,dominant,h3k4me3,/data/path/Res_H3K4me3_R1.fq.gz,/data/path/Res_H3K4me3_R2.fq.gz
sus1,Sus,pool,recessive,h3k27me3,/data/path/Sus_H3K27me3_R1.fq.gz,/data/path/Sus_H3K27me3_R2.fq.gz
sus2,Sus,pool,recessive,h3k36me3,/data/path/Sus_H3K36me3_R1.fq.gz,/data/path/Sus_H3K36me3_R2.fq.gz
sus3,Sus,pool,recessive,h3k4me3,/data/path/Sus_H3K4me3_R1.fq.gz,/data/path/Sus_H3K4me3_R2.fq.gz
```

- `sample`：样本名称，需要保证每个样本名唯一；
- `bulk`：混池名称，属于一个混池的不同样本该列名称需要一致，如果包含亲本数据，则亲本和对应混池的名称也需一致；
- `type`：样本数据类型，该列可以填写 **pool**（对应混池数据）或者 **parent**（对应亲本数据）；
- `genotpe`：性状相关基因的基因型，可以填写 **recessive**（对应隐性）或者 **dominant**（对应显性），如果不确定则全部填写 **unknown**；
- `marker`：测序数据修饰类型，如 h3k36me3、h3k27me3、rnaseq 和 wgs 等，该列用于区分同一混池的不同样本；
- `fq1` 和`fq2`：测序文件 R1 和 R2 的具体路径，需要使用绝对路径。

通过上述示例样本信息文件可以看出该组输入文件为不包含亲本的 2 组混池数据，其中每组混池包含 3 种不同组蛋白修饰的 ChIP-seq 数据。

一个包含亲本的 WGS 混池数据样本信息文件示例如下所示：

```
sample,bulk,type,genotype,marker,fq1,fq2
res1,Res,pool,dominant,wgs,/data/path/Res_wgs_R1.fq.gz,/data/path/Res_wgs_R2.fq.gz
sus1,Sus,pool,recessive,wgs,/data/path/Sus_wgs_R1.fq.gz,/data/path/Sus_wgs_R2.fq.gz
resp,Res,parent,dominant,wgs,/data/path/data/testp_Res_wgs_R1.fq.gz,/data/path/data/testp_Res-wgs_R2.fq.gz
susp,Sus,parent,recessive,wgs,/data/path/data/testp_Sus_wgs_R1.fq.gz,/data/path/data/testp_Sus_wgs_R2.fq.gz
```

一个不包含亲本的 RNA-seq 混池数据样本信息文件示例如下所示：

```
sample,bulk,type,genotype,marker,fq1,fq2
res1,Res,pool,dominant,rnaseq,/data/path/Res_rnaseq_R1.fq.gz,/data/path/Res_rnaseq_R2.fq.gz
sus1,Sus,pool,recessive,rnaseq,/data/path/Sus_rnaseq_R1.fq.gz,/data/path/Sus_rnaseq_R2.fq.gz
```

### 填写过滤区间文件 region.csv

**注意：该文件在仅使用 Assembly Module 时需要**

使用 Triti-Map 进行 _de novo_ 序列组装时，需要使用已知的性状关联区间对序列进行过滤，如果运完整分析模块该数据会自动生成，如果仅运行 Assembly Module 模块时需要用户在`region.csv`配置文件中提前填写过滤区间。

一个典型的过滤区间信息文件 region.csv 如下所示

```
chrom,start,end
chr7A,724000000,730000000
```

- `chrom`：定位区间所在染色体；
- `start`：定位区间所在染色体的起始位置；
- `end`：定位区间所在染色体的终止位置

### 填写配置文件 config.yaml

完成上述准备工作后，可以在配置文件 config.yaml 中填写 Triti-Map 运行时需要的分析参数。

**配置文件关键参数说明**

- `email`：必填，在使用 EMBL-EBI API 进行相关分析时需要提供个人邮箱。
- `samples`：无需修改，指定样本信息文件路径。默认样本信息文件为 Triti-Map 运行目录内的 `sample.csv`。
- `datatype`：必填，指定样本数据的测序类型，`dna`：ChIP-seq 或 WGS；`rna`： RNA-seq。
- `maxthreads`: 指定可以使用的最大线程数，默认为`30`。
- `ref`：该父参数下的三个子参数分别为参考基因组文件路径、STAR 参考基因组 index 目录（进行 RNA-seq 分析时需要）、基因组注释文件路径
- `gatk`: 该父参数下的`min_SNP_DP` 表示过滤 SNP 时每个混池需要满足的最小深度。该深度数据和样本测序数据相关，通常建议设置为 10 即可。
- `snpindex`：该父参数下的子参数为进行 BSA 定位分析时 QTLseqr 软件需要指定的参数以及对分析结果进行过滤时所需要的参数。
  - `pop_struc`：该参数为混池样本的群体结构。QTLseqr 支持 `F2` 和 `RIL` 两个参数，如果混池数据为 F2 代则填写 F2 如果混池数据为 RIL 群体或所有个体均为纯合则填写 RIL。
  - `bulksize`：该参数为 QTLseqr 分析是需要填写的混池样本数量，例如两个混池各自有 30 个样本，则填写 30。（无需绝对准确）
  - `winsize`：该参数为 QTLseqr 进行数据矫正时的滑窗长度，建议使用 1e6(1Mb) 即可，无需修改。
  - `filter_probs`：根据 delta snpindex 和 snpconut/1Mb 过滤的百分比，如果该值为 0.75，则表示定位区间的平均 snpindex 和定位区间平均每 Mb 的 SNP 数量需要同时大于 qtlseqr 所有原始结果的 75%。
  - `fisher_p`: 过滤性状关联区间的 snp 信息，使用 fisher test 计算每个突变位点的可信度。
- `merge_lib`：该参数用于指定在进行 de novo 组装时如何处理一个混池内的多组不同组蛋白修饰 ChIP-seq 数据，建议设置为`merge`，即先进行样本合并再进行一次组装，可以得到较好的组装结果；如果在处理六倍体小麦等大基因组数据且服务器内存小于 300G 时，可以填写`split`，将会对每一组数据单独进行组装然后再将各自组装结果进行合并进行后续分析。
- `memory`：该参数指定进行转录组序列组装时可用的最大内存，`400` 表示 400G。
- `abyss`：该参数下的自参数`kmer`表示进行 ChIP-seq 或 WGS 数据组装时，abyss 软件使用的 kmer，已有分析经验显示在分析六倍体小麦等大基因组数据时设置为 `90` 效果较好同时不会占用大量服务其内存。
- `blast_database`：该参数表示在进行 de novo 组装序列进行相似序列比对注释时使用的数据库，`em_cds_pln` 使用 EBI 的 ENA 植物编码序列数据库，`em_std_pln` 使用 EBI 的 ENA 植物标准序列数据库。这里建议只需要使用`em_cds_pln`皆可，`em_std_pln`会造成分析时间过长和 API 返回数据连接不稳定等问题。

配置文件 config.yaml 如下所示，`#`后为填写说明

```yaml
# The contents after the '#' symbol are the option description.

###############################
###### Basic Information ######
###############################

# use EMBL-EBI hmmerscan and blast API need user email, required by EMBL-EBI API.
email: zhaofei920810@gamil.com

# sample information file path. No change needed
samples: sample.csv

# datatype option: dna/rna
# dna:ChIP-seq or WGS data
# rna:RNA-seq
datatype: dna

# the max threads that Triti-Map can use
maxthreads: 30

##########################################
###### Reference Genome Information ######
##########################################

# genome reference and annotation file path
ref:
  genome: /genome/path/161010_Chinese_Spring_v1.0_pseudomolecules.fasta
  annotation: /annotation/path/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gtf
  # STARdir is required if datatype is RNA-seq(rna)
  STARdir: /star/path/wheat_STAR_wholegenome1.1

#####################################
###### Raw Fastq Data Pretreat ######
#####################################

# fastp software parameters. No change needed!
trim:
  quality: 20
  minlength: 50

#########################################
###### BSA-based Module Parameters ######
#########################################

# min SNP depth in each bulk to filter gatk raw vcf file, related to sequencing depth，10 is recommended. No change needed!
gatk:
  min_SNP_DP: 10

# Important!
snpindex:
  # population structure
  # pop_struc option: RIL/F2, required by QTLseqr
  pop_struc: RIL
  # bulk size, required by QTLseqr
  bulksize: 30
  # smooth window size, required by QTLseqr. No change needed!
  winsize: 1000000
  # delta snpindex and snp-conut/1Mb percentage of filtration.
  # 0.75: candidate region with average delta snpindex and snp-conut/1Mb both greater than 75% of all raw regions will be retained.
  filter_probs: 0.75
  # SNP fesher test filter pvalue, used to filter candidate SNP. No change needed!
  fisher_p: 0.0001

########################################
###### Assemble Module Parameters ######
########################################

# If your chip-seq data have multiple histone modification types, option 'split' is recommended when the server memory is less than 300G.
# merge_lib option: merge/split
merge_lib: merge
# the memory(G) spades software can use when assemble RNA-seq data
memory: 400
# the kmer abyss software can use when assemble ChIP-seq data, 90 is recommended. No change needed!
abyss:
  kmer: 90
# the de novo scaffolds min length. No change needed!
scaffold:
  minlength: 500

# de novo scaffolds filter method: Triti-Map/external_region/external_fasta
# Triti-Map: auto filter by Triti-Map method, don't need external file
# external_region: filter by user-defined region file
# external_fasta: filter by user-defined fasta file
denovo_filter_method: Triti-Map
# if denovo_filter_method: external_region, set region file path
filter_region_file: ""
# if denovo_filter_method: external_fasta, set fasta file path
filter_fasta_file: ""

# Database to be used for de novo scaffolds blast
# blast_database option: em_cds_pln/em_std_pln
# em_cds_pln (plant cds database) is highly recommended, em_std_pln (plant standard database) will consume a lot of time
blast_database: em_cds_pln
```

## 运行 Triti-Map

Triti-Map conda 运行环境配置完成，相关物种文件下载并建立 index 文件后，样本信息文件和配置文件填写完成即可以运行 Triti-Map 主程序。

- `module`: 必填，Triti-Map 分析模块，`only_bas` 只进行 BSA 定位分析；`only_assemble` 只进行 de novo 序列组装；`all` 同时进行上述两种分析。

首先测试配置文件是否正确且流程是否可以正常运行。

```sh
# 首先进入 Triti-Map 目录内
cd Triti-Map
# 确保启用 Triti-Map conda 运行环境
conda activate Triti-Map
# 使用 snakemake -n 参数进行伪运行，使用 -j 参数指定最大运行线程数
snakemake -j 30 -n
# 生成DAG图

tritimap run -d test -j 30 -n all --dag | dot -Tpdf > dag_rna.pdf
```

正常情况会输出 Triti-Map 将运行的任务名称和任务数量，如下所示：

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
	1	all
	2	dnaBWAmem2Mapping
	2	dnaGATK4Markdup
	22	gatk4HC
	1	uniqScaffold
	1	vcf2tab
	1	vcfGenotypeFiltration
	1	vcfHardFiltration
	42
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

没有其它错误信息即可进行正式分析，因为完成的流程分析时间较长建议将任务使用 `screnn` 或其它类似程序置于服务器后台运行。

```sh

# 进入 Triti-Map 目录内
cd Triti-Map
# 确保启用 Triti-Map conda 运行环境
conda activate Triti-Map
# 使用 snakemake -j 参数指定最大运行线程数
snakemake -j 30

```

在分析过程中如果出现错误可以通过报错信息定位出现错误的分析步骤，通过报错信息和对应步骤的 log 文件排除错误原因。

## Triti-Map 结果说明

一个完整的 Triti-Map 流程结果目录如下图所示：

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

每个目录内的文件信息如下，产生的中间过程文件可以直接用于其它分析：

- `01_cleandata` 原始数据经过预处理的清洁数据
- `02_mergedata` 不同修饰类型数据合并后用于后续分析的混池数据
- `03_mappingout` 与参考基因组比对后的 bam 文件
- `04_GATKout` 经过 GATK 预处理后用于鉴定变异的 bam 文件
- `05_vcfout` GATK 处理生成的变异信息原始 VCF 文件和经过过滤筛选的 VCF 文件
- `06_regionout` 用于进行 BSA 分析的 snp 输入文件和 BSA 分析结果文件
- `07_assembleout` de novo 组装生成的混池序列和功能注释文件
- `logs` 分析过程中产生的日志文件

如果分析过程中有中间文件生成，每一个子目录中的中间临时文件都会保存在子目录的`temp_output`目录中，中间临时文件正常情况下无需使用，如果流程出现报错可以结合日志文件进行检查。

### 基于 BSA 定位模块的主要结果文件

基于 BSA 定位模块的分析结果位于目录 06_regionout, 文件名中的`xxx`和`yyy`表示输入文件中的混池名称。

```
06_regionout
├── xxx_yyy_snpindex_input.txt
├── xxx_yyy_qtlseqr_output.txt
├── xxx_yyy_qtlseqr_raw_region.txt
├── xxx_yyy_qtlseqr_filter_region.txt
├── xxx_yyy_qtlseqr_filter_snpinfo.txt
├── xxx_yyy_qtlseqr_SNPcounts_line.pdf
├── xxx_yyy_qtlseqr_SNPcounts_point.pdf
├── xxx_yyy_qtlseqr_SNPindex_line.pdf
├── xxx_yyy_qtlseqr_SNPindex_point.pdf
└── xxx_pool_vs_yyy_pool_candidateregion_1.pdf
```

- `xxx_yyy_snpindex_input.txt` 所有用于 BSA 分析的 SNP 信息
- `xxx_yyy_qtlseqr_output.txt` qtlseqr 通过 delta SNPindex 计算后得到的所有 SNP 结果信息
- `xxx_yyy_qtlseqr_raw_region.txt` qtlseqr 通过 delta SNPindex 计算后得到的所有定位区间结果信息
- `xxx_yyy_qtlseqr_filter_region.txt` 通过筛选过后的高可信度定位区间信息
- `xxx_yyy_qtlseqr_filter_snpinfo.txt` 高可信度定位区间内的 SNP 结果信息
- `*.pdf` pdf 文件为分析结果可视化图片，包括全基因 delta SNPindex、全基因组 snp counts 和定位区间 delta SNPindex

### 基于 de novo 序列组装模块的主要结果文件

基于 de novo 序列组装的结果位于 07_assembleout 中，文件名中的`xxx`和`yyy`表示输入文件中的混池名称。

```
07_assembleout
├── xxx_merge_denovo_scaffolds.fasta
├── xxx_candiadate_denovo.fasta
├── xxx_candidate_denovo2ref.info.txt
├── xxx_candiadate_denovo_pfam_anno.fasta
├── xxx_candiadate_denovo_pfam_anno.txt
├── xxx_candiadate_denovo_blast_anno.txt
├── yyy_merge_denovo_scaffolds.fasta
├── yyy_candiadate_denovo.fasta
├── yyy_candidate_denovo2ref.info.txt
├── yyy_candiadate_denovo_pfam_anno.fasta
├── yyy_candiadate_denovo_pfam_anno.txt
└── yyy_candiadate_denovo_blast_anno.txt
```

- `*_merge_denovo_scaffolds.fasta` 每个混池组装后各自的特异性序列
- `*_candiadate_denovo.fasta` 每个混池组装后各自可以和性状关联区域部分匹配的特异性序列
- `*_candidate_denovo2ref.info.txt` 可以和性状关联区域部分匹配的特异性序列与参考基因组对应的位置信息
- `*_candiadate_denovo_pfam_anno.txt` 可以和性状关联区域部分匹配特异性序列的功能注释结果
- `*_candiadate_denovo_pfam_anno.fasta` 可以和性状关联区域部分匹配且具有功能的特异性序列
- `*_candiadate_denovo_blast_anno.txt` 可以和性状关联区域部分匹配且具有功能的特异性序列在 EBI 植物数据库中的相似性比对注释信息

### 具体结果文件格式说明

`xxx_yyy_snpindex_input.txt` 文件示例如下：

```
#CHROM  POSm1   POS     REF     ALT     Res_pool_ref    Res_pool_alt    Res_pool_depth  Res_pool_ratio  Sus_pool_ref    Sus_pool_alt    Sus_pool_depth  Sus_pool_ratiosnpindex
chr1A   1144541 1144542 G       A       13      147     160     0.91875 1       144     145     0.993103        -0.0743534
chr1A   1144550 1144551 A       G       13      154     167     0.922156        1       153     154     0.993506        -0.0713508
chr1A   1669712 1669713 G       A       13      4       17      0.235294        21      1       22      0.0454545       0.18984
```

从左到右每一列数据依次表示：突变所在染色体、所在位点-1、所在位点、参考基因组碱基、突变碱基、混池 A 的参考基因组碱基数、混池 A 的 突变碱基数、混池 A 的 SNPindex、混池 B 的参考基因组碱基数、混池 B 的突变碱基数、混池 B 的 SNPindex、该位点的 delta SNPindex。

`xxx_yyy_qtlseqr_filter_region.txt` 文件示例如下：

```
CHROM   qtl     start   end     length  nSNPs   avgSNPs_Mb      peakDeltaSNP    posPeakDeltaSNP avgDeltaSNP
chr7A   7       724111912       730119678       6007766 3262    543     0.954374003095134       729506093       0.812296971085005

```

从左到右每一列数据依次表示：定位区间所在染色体、qtlseqr 原始结果对应区间编号、区间起始位置、区间终止位置、区间长度、区间内 SNP 个数、每 1M 内的平均 SNP 个数、区间内 Delta SNPindex 峰值、区间内 Delta SNPindex 峰值所在位置、区间内平均 Delta SNPindex。具体信息可参考 QTLseqr。

`xxx_yyy_qtlseqr_filter_snpinfo.txt` 文件示例如下：

```
CHROM   POS     REF     ALT     AD_REF.Sus_pool AD_ALT.Sus_pool DP.Sus_pool     SNPindex.Sus_pool       AD_REF.Res_pool AD_ALT.Res_pool DP.Res_pool     SNPindex.Res_pool      REF_FRQ deltaSNP        nSNPs   tricubeDeltaSNP minDP   tricubeDP       CI_99   fisher_p
chr7A   724111912       T       G       0       12      12      1       21      0       21      0       0.636363636363636       -1      445     0.505705636608414     12       20      -0.5    2.8183517084228e-09
chr7A   724111983       G       C       0       14      14      1       16      0       16      0       0.533333333333333       -1      445     0.505753155444727     14       20      -0.5    6.87650670708678e-09
```

从左到右每一列数据依次表示：突变所在染色体、突变所在位点、参考基因组碱基、突变碱基、混池 A 的参考基因组碱基数、混池 A 的 突变碱基数、混池 A 的深度、混池 A 的 SNPindex、混池 B 的参考基因组碱基数、混池 B 的 突变碱基数、混池 B 的深度、混池 B 的 SNPindex、、qtlseqr 计算时使用的参考等位基因频率、delta SNPindex 、位点深度、矫正后的 delta SNPindex、最小深度、矫正后的位点深度、99%置信区间的 SNPindex 值、fisher test 的 p value

`*_candiadate_denovo.fasta` 文件示例如下：

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

de novo 组装序列示例 ID：`17655292:chr1A:2501389:378S202M:580:res`
序列 ID 命名规则说哦：`原始组装序列ID:比对基因组序列所在染色体:比对基因组序列对应染色体位置:比对相关的CIGAR信息:序列长度:对应混池`

`*_candiadate_denovo_pfam_anno.txt` 文件示例如下：

```
seqid   Description    E-value       Model
23558405:chr7B:733811080:14S493M9I1522M:2038:res        NB-ARC NB-ARC domain    2.3e-37       PF00931.22
```

从左到右每一列数据依次表示：序列 ID、注释功能描述、domain 比对 E-value、pfam domain ID。

`candiadate_denovo_blast_anno.txt` 文件示例如下：

```
seqid	hit_db	hit_id	hit_desc	hit_url	hsp_bit_score	hsp_align_len	hsp_identity	hsp_query_from	hsp_query_to	hsp_hit_from	hsp_hit_to
23558405:chr7B:733811080:14S493M9I1522M:2038:res	EM_CDS	AUO29721.1	Triticum urartu powdery mildew resistance protein	https://www.ebi.ac.uk/ena/data/view/AUO29721.1	3514.25	1958	99.8	81	2038	1	1958
```

从左到右每一列数据依次表示：序列 ID、比对数据库、数据库序列 ID、数据库序列描述、数据库序列信息 URL、比对分数、比对长度、比对 identity、组装序列比对起点、组装序列比对终点、数据库序列比对起点、数据库序列比对终点。

### 可视化结果说明

**全基因组 delta SNPindex**

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210310150108.png)

纵坐标表示`delta SNPindex`，横坐标染色体位置。蓝色线表示 QTLseqR 计算得出的 99% 置信区间对应 delta SNPindex 数值，灰颜色点表示对应的 SNP 低于此值，黑颜色的点表示对应的 SNP 高于此值。

**全基因组 SNP count**

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210310150210.png)

纵坐标表示 SNPcount,横坐标表示染色体位置。

**性状关联区间局部 delta SNPindex** （每个性状关联区间都会单独生成对应的可视化结果）

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210310150307.png)

该图为全基因组 delta SNPindex 的性状关联区间局部放大图。其中纵坐标表示`delta SNPindex`，横坐标染色体位置，蓝色线表示 QTLseqR 计算得出的 99% 置信区间对应 delta SNPindex 数值。红色部分表示该位置为计算出的性状关联区间。

## 在线辅助平台后续分析

Triti-Map 的分析结果可以配套的[在线辅助分析平台](http://bioinfo.sibs.ac.cn/wheatsnpanno)进行更多后续分析和可视化展示。如对 SNP 信息进行功能和表观修饰注释，对基因进行小麦族同源基因分析，对性状关联区间进行共线性分析，对特异性片段进行功能注释等。可访问[在线平台](http://bioinfo.sibs.ac.cn/wheatsnpanno)阅读详细说明。

在线分析平台主要功能

![](https://kaopubear-1254299507.cos.ap-shanghai.myqcloud.com/picgo/20210330151510.png)

## 支持

### 更新日志

### 常见问题

Q：可以分析哪些物种？

A：Triti-Map 分析流程最初针对以六倍体小麦为代表的小麦族等大基因组物种设计，其中比对和突变检测步骤针对大基因组作物进行了优化，部分结果可视化方式也根据大基因组作物进行了调整。

但是原则上 Triti-Map 可以对任何物种的混池测序数据进行 BSA 分析和性状相关 de novo 序列分析。需要注意的是，Triti-Map 不适合处理含有大量 scaffold 的基因组数据，如果下载的基因组文件除了染色体以外还有大量 scaffold 序列，建议使用 `samtools` 的 `samtools faidx` 命令或者 `seqkit` 的`seqkit grep`命令首先提取出染色体级别的基因组作为参考基因组文件然后再进行分析。

Q：可以分析哪些数据？

A：针对高通量测序类型。Triti-Map 分析流程最初针对 ChIP-seq 数据设计，但是原则上 BSA 模块可以分析包括全基因组重测序(WGS)、外显子组测序和转录组测序在内的大部分高通量测序数据。但是需要注意的是，一方面外显子组测序数据使用 de novo 组装模块意义不大，另一方面使用大基因组物种的 WGS 数据进行 de novo 组装模块分析时，需要足够的计算资源且耗时较长。

针对测序材料类型，因为 Triti-Map 的 BSA 模块目前仅支持使用 QTLseq delta snpindex 方法进行定位区间筛选，并未针对 EMS 诱变数据设置特殊的过滤条件，所以目前更适合通过杂交重组获得特殊性状的材料。其中 de novo 组装模块尤其适合类似于通过不同倍性小麦杂交获得优良基因的情况，这类基因往往和参考基因组相比存在着大片段的变异或者本身在参考基因组中不存在。

Q：使用 ChIP-seq 数据时应该使用哪些修饰。

A：经测试，H3K27me3、H3K4me3 组合更适用与进行 BSA 定位分析，进行 de novo 序列组装以挖掘新基因时建议增加 H3K36me3。

Q：如何使用指定区间对*de novo*组装序列进行过滤

A：如果已知的定位区间长度过小（小于 1Mb），建议适当延伸起始位点和终止位点，增加区间长度。同时，如果是多倍体小麦，考虑到小麦亚基因组之间存在的变异情况，可以首先通过[在线辅助分析平台](http://bioinfo.sibs.ac.cn/wheatsnpanno)获得对应区间在其它同一物种不同亚基因组的共线性区域，一并填写到该定位区间信息文件中。

Q：如何使用指定数据库对*de novo*组装序列进行过滤

A：

Q: filter_probs 应该如何设置

A：

Q：需要怎样的服务器性能？

通过已有材料的分析经验表明，Triti-Map 分析流程中最消耗内存的是 de novo 模块中的序列组装步骤。以六倍体小麦 400M reads 的混池 ChIP-seq 数据为例，如果使用 de novo 组装模块，建议服务器内存不少于 400G。

Q：如果只进行 de novo 组装分析需要注意哪些问题？
A：如果已经通过已有的 BSAseq 相关分析或传统的图位克隆获得了性状相关基因所在区段但是在已有的参考基因组中并未鉴定到和性状相关的已知基因时，可以仅使用 de novo 组装模块进行分析。这时需要在`region.csv`中设置已知的一个或多个性状相关区间，如果定位区间少于 1Mb，可以适当增加区间长度不少于 1Mb。

## 引用

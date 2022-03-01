# Triti-Map 中文使用说明

## Triti-Map 介绍

<p align="center">
    <img src="https://kaopubear-1254299507.file.myqcloud.com/picgo/20210401150558.png">
</p>

[Triti-Map](https://github.com/fei0810/Triti-Map) 是一套结合多组学 BSA 定位分析和*de novo*序列组装的小麦族功能基因定位分析流程。通过 [Snakemake](https://snakemake.github.io/) 开发，共包含三个分析模块：

- Interval Mapping Module
- Assembly Module
- Web-based Annotation Module

该分析流程针对小麦族大基因组物种开发和优化，进行常规 BSA 定位分析的同时可以挖掘参考基因组中不存在的新功能基因，还可以结合在线注释平台进行更多深入的下游分析。

其中，Interval Mapping Module 和 Assembly Module 为命令行软件，输入为混池测序 DNA-Seq（ChIP-Seq/WGS）或 RNA-Seq 数据，可以一步生成包括性状关联区间、突变位点和新基因在内的多种结果。

Web-based Annotation Module 为[在线分析平台](http://bioinfo.cemps.ac.cn/tritimap/)。Triti-Map 收集了小麦族各物种和六倍体小麦已有测序品种的基因组信息。统一进行了基因功能注释和转录因子结合位点预测，同时整合了各物种有代表性的表观修饰数据。该平台可以进行包括 SNP 注释与可视化展示、同源基因分析、小麦族共线性区间分析和新序列功能注释在内的多种分析，可以为小麦族基因克隆提供更加丰富的参考信息。

### Triti-Map 分析流程图

![](https://kaopubear-1254299507.file.myqcloud.com/picgo/20210330151508.png)

### Triti-Map 具体分析方法和针对性优化内容

![](https://kaopubear-1254299507.file.myqcloud.com/picgo/20210330151509.png)

## 作者信息

开发者：zhaofei920810@gmail.com

实验室主页： http://bioinfo.cemps.ac.cn/zhanglab/

## 安装 Triti-Map

### 安装 Conda

运行 Triti-Map 首先需要在 UNIX 环境中安装 Conda（建议使用 [Miniconda](https://conda.io/en/latest/miniconda.html)），同时保证可以正常使用[Bioconda](http://bioconda.github.io/)。

```sh
# 安装 Miniconda
wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### 使用 Conda 安装 Triti-Map

新建 Triti-Map 运行环境并通过 bioconda 安装 Triti-Map（推荐）。

```sh
# 新建运行环境同时安装Triti-Map
conda create -c conda-forge -c bioconda -n tritimap tritimap
# 启动运行环境
conda activate tritimap
# 测试是否安装成功
tritimap --help
```

正常安装后可以看到如下类似的帮助信息：

```
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

说明：使用独立环境运行 Triti-Map，可以避免其它软件的影响。

### 使用源码安装 Triti-Map

你也可以通过 github 下载 Triti-Map 源码后自行安装 Triti-Map 和其它依赖软件。

```sh
# 下载Triti-Map源码
git clone https://github.com/fei0810/Triti-Map.git
# 进入Triti-Map目录
cd Triti-Map
# 安装Triti-Map
python setup.py install
# 安装依赖软件
```

说明：使用源码进行安装时，需要自行安装 Triti-Map 其它依赖软件，可以通过目录中的`tritimap_env.yaml`配合`conda`进行安装。

### 使用 Docker 安装使用 Triti-Map

你也可以通过 Docker 安装使用 Triti-Map, 查看 [如何安装 Docker](https://docs.docker.com/engine/install/)，然后通过下面的命令安装 Triti-Map:

```sh
# 下载 Triti-Map
docker pull fei0810/tritimap:v0.9.7
# 运行 Docker 镜像
docker run -i -t fei0810/tritimap:v0.9.7 /bin/bash
```

## 准备相关文件

### 基因组和注释文件

Triti-Map 针对 DNA-seq（ChIP-Seq 和 WGS 数据）数据使用 bwa-mem2 进行比对，针对 RNA-Seq 数据使用 STAR 进行比对，运行 Triti-Map 前需要准备对应的基因组文件和 gff3 注释文件。进行正式分析前需要准备相关基因组和注释文件，并构建索引文件。

小麦族部分已测序物种基因组和注释文件可通过以下链接下载。

- [中国春 _Triticum aestivum_](http://plants.ensembl.org/Triticum_aestivum/Info/Index)
- [硬粒小麦 _Triticum turgidum_](http://plants.ensembl.org/Triticum_turgidum/Info/Index)
- [野生二粒小麦 _Triticum dicoccoides_](https://plants.ensembl.org/Triticum_dicoccoides/Info/Index)
- [大麦 _Hordeum vulgare_](http://plants.ensembl.org/Hordeum_vulgare/Info/Index)
- [粗山羊草 _Aegilops tauschii_](https://plants.ensembl.org/Aegilops_tauschii/Info/Index)
- [乌拉尔图 _Triticum urartu_](http://www.mbkbase.org/Tu/)

**生成基因组 index 文件**

```sh
# 启动 Triti-Map 运行环境
conda activate tritimap
# 建立参考基因组的 gatk 索引
gatk CreateSequenceDictionary -R /genome/path/genome.fasta
# 建立参考基因组的 samtools 索引
samtools faidx /genome/path/genome.fasta
```

**生成基因组比对软件的 index 文件**

DNA-seq 数据使用 [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)。

```sh
# 启动 Triti-Map 运行环境
conda activate tritimap
# 使用 bwa-mem2 建立参考基因组比对index文件
bwa-mem2 index /genome/path/genome.fasta
```

RNA-seq 数据使用 [STAR](https://github.com/alexdobin/STAR)。

```sh
# 启动Triti-Map 运行环境
conda activate tritimap
# 使用 STAR 建立
STAR --runThreadN 30 --runMode genomeGenerate \
--genomeDir /star/index/path/genome_star \
--genomeFastaFiles /genome/path/genome.fasta \
--sjdbOverhang 100 \
--sjdbGTFfile /anntotaion/path/genome.gtf \
--genomeChrBinNbits 18 \
--limitGenomeGenerateRAM 50805727274
```

使用 STAR 构建六倍体小麦基因组索引时需要通过参 `--limitGenomeGenerateRAM` 提高可用内存。具体数值不同基因组略有不同。

在运行过程中 Triti-Map 会根据填写的基因组文件路径自动检测是否有正确匹配的索引文件。

### 配置文件

使用 `tritimap init` 生成主程序运行所需要的配置文件。

```sh
# 启动 Triti-Map 运行环境
conda activate tritimap
# 进入运行目录
cd /your/work/path
# 生成配置文件
tritimap init
```

`tritimap init` 运行成功后所在运行目录会生成三个配置文件：

- config.yaml：Triti-Map 主程序运行参数配置文件
- sample.csv：原始分析数据样本信息文件
- region.csv：序列组装过滤区间文件（仅在单独使用`Assembly Module`时需要）

#### 修改样本信息文件 sample.csv

样本信息文件为`,`分割的 CSV 文件，包括数据的样本名、测序类型、混池基因型和文件位置等信息。

Triti-Map 可以接受仅有混池的 2 组输入数据，或同时包含混池和亲本的 4 组输入数据。为了提高突变鉴定和组装效果，使用 ChIP-seq 数据时建议同时使用 2-3 种不同组蛋白修饰数据，如 H3K27me3、H3K4me3 和 H3K36me3。

**填写说明**

- `sample`：样本名称，需要保证每个样本名唯一；
- `bulk`：混池名称，属于同一混池的不同样本该列名称需要一致，如包含亲本，则亲本和对应混池的名称也需一致；
- `type`：样本数据类型，该列可以填写 **pool**（对应混池数据）或者 **parent**（对应亲本数据）；
- `genotpe`：性状相关基因的基因型，可以填写 **recessive**（对应隐性）或者 **dominant**（对应显性），如果不确定则全部填写 **unknown**；
- `marker`：测序类型，如 h3k36me3、h3k27me3、rnaseq 和 wgs 等，该列用于区分同一混池的不同样本；
- `fq1` 和`fq2`：测序文件 R1 和 R2 的具体路径，**需要使用绝对路径**。

**示例 1**

一个典型的多种修饰 ChIP-seq 样本信息文件`sample.csv`如下所示。输入文件为不包含亲本的 2 组混池数据，其中每组混池包含 3 种不同组蛋白修饰的 ChIP-seq 数据。

```
sample,bulk,type,genotype,marker,fq1,fq2
res1,Res,pool,dominant,h3k27me3,/data/path/Res_H3K27me3_R1.fq.gz,/data/path/Res_H3K27me3_R2.fq.gz
res2,Res,pool,dominant,h3k36me3,/data/path/Res_H3K36me3_R1.fq.gz,/data/path/Res_H3K36me3_R2.fq.gz
res3,Res,pool,dominant,h3k4me3,/data/path/Res_H3K4me3_R1.fq.gz,/data/path/Res_H3K4me3_R2.fq.gz
sus1,Sus,pool,recessive,h3k27me3,/data/path/Sus_H3K27me3_R1.fq.gz,/data/path/Sus_H3K27me3_R2.fq.gz
sus2,Sus,pool,recessive,h3k36me3,/data/path/Sus_H3K36me3_R1.fq.gz,/data/path/Sus_H3K36me3_R2.fq.gz
sus3,Sus,pool,recessive,h3k4me3,/data/path/Sus_H3K4me3_R1.fq.gz,/data/path/Sus_H3K4me3_R2.fq.gz
```

**示例 2**

包含亲本的 WGS 混池数据：

```
sample,bulk,type,genotype,marker,fq1,fq2
res1,Res,pool,dominant,wgs,/data/path/Res_wgs_R1.fq.gz,/data/path/Res_wgs_R2.fq.gz
sus1,Sus,pool,recessive,wgs,/data/path/Sus_wgs_R1.fq.gz,/data/path/Sus_wgs_R2.fq.gz
resp,Res,parent,dominant,wgs,/data/path/data/testp_Res_wgs_R1.fq.gz,/data/path/data/testp_Res-wgs_R2.fq.gz
susp,Sus,parent,recessive,wgs,/data/path/data/testp_Sus_wgs_R1.fq.gz,/data/path/data/testp_Sus_wgs_R2.fq.gz
```

不包含亲本的 RNA-seq 混池数据：

```
sample,bulk,type,genotype,marker,fq1,fq2
res1,Res,pool,dominant,rnaseq,/data/path/Res_rnaseq_R1.fq.gz,/data/path/Res_rnaseq_R2.fq.gz
sus1,Sus,pool,recessive,rnaseq,/data/path/Sus_rnaseq_R1.fq.gz,/data/path/Sus_rnaseq_R2.fq.gz
```

#### 修改过滤区间文件 region.csv

**注：该文件仅在单独使用 Assembly Module 时才需要修改。**

使用 Triti-Map 进行 _de novo_ 序列组装时，需要使用性状关联区间对组装序列进行过滤，如果运行完整分析模块该区间会自动生成；如果仅运行 Assembly Module 模块时需要用户在`region.csv`配置文件中提前填写过滤区间。

一个典型的过滤区间信息文件 region.csv 如下所示：

```
chrom,start,end
chr7A,724000000,730000000
```

- `chrom`：定位区间所在染色体
- `start`：定位区间起始位置
- `end`：定位区间终止位置

#### 修改配置文件 config.yaml

运行 Triti-Map 主程序需要在配置文件 config.yaml 中修改各项参数。

**配置文件重要参数说明**

- `email`：**重要**，使用 EMBL-EBI API 进行相关分析时需要提供个人邮箱。
- `samples`：无需修改，样本信息文件路径。默认样本信息文件为 Triti-Map 运行目录的 `sample.csv`。
- `datatype`：**重要**，样本数据的测序类型，`dna`：ChIP-seq 或 WGS；`rna`： RNA-seq。
- `maxthreads`: 可以使用的最大线程数，默认为`30`。
- `ref`：**重要**，参考基因组相关文件路径，包含 3 个子参数。
  - `genome`：参考基因组文件路径，使用绝对路径。
  - `annotation`：参考基因组基因注释文件路径，使用绝对路径。
  - `TARdir`：STAR 参考基因组 index 目录（进行 RNA-seq 分析时需要填写）。
- `gatk`: GATK 相关分析参数，包含 1 个子参数。
  - `min_SNP_DP` 无需修改，有效 SNP 每个混池需要满足的最小深度。该深度和样本测序数据量相关，默认为 10。
- `snpindex`：进行 BSA 定位分析时需要的各项参数，包含 6 个子参数。
  - `pop_struc`：**重要**，混池样本的群体结构。如果混池数据为 F2 代则填 **F2** ；如果混池数据为 RIL 群体或所有个体均为纯合则填写 **RIL**。
  - `bulksize`：混池样本数量，如每个混池由 30 个样本混合，则填写 **30**。
  - `winsize`：无需修改，进行数据矫正时的滑窗长度，默认为 1000000(1Mb)。
  - `filter_percentage`：**重要**，根据 Delta snpindex 和 snpconut/Mb 对原始定位结果进行过滤的百分比。如果该值为 0.75，则表示候选定位区间的平均 Delta snpindex 和每 1Mb 的平均 SNP 数量要同时大于所有原始结果对应数值的 75%。
  - `fisher_p`: 无需修改，过滤性状关联区间的 SNP 位点，使用 fisher test 计算每个位点的 pvalue，默认为 0.0001。
  - `min_length`：候选性状关联区间的最小长度，针对小麦族等大基因组物种，默认为 1000000(1Mb)，无需修改。
- `bulk_specific`: **重要**，如何定义混池特异性序列
  - `identical_percentage`: 两组混池组装序列 blast 比对结果的 identical 过滤值，默认为 0.85，表示 identical 大于 85%的序列将会被过滤掉。
  - `length_percentage`：两组混池组装序列 blast 比对结果中比对长度占序列总长的比例，默认为 0.85，表示比对长度超过总长 85% 的序列将会被过滤掉。
- `merge_lib`：如何处理同一混池的多组不同修饰数据。**默认为`merge`**，即先进行样本合并再进行组装，可以得到更好的结果；如果处理六倍体小麦等大基因组数据且服务器内存小于 300G 时，可以修改为`split`，即对每一组数据单独组装再合并进行后续分析。
- `memory`：转录组序列使用 SPAdes 进行组装时可用最大内存，`300` 表示 300G。
- `denovo_filter_method`: **重要**，混池特异性序列过滤方式，在仅运行组装模块(only_assembly)时有效。设置为`external_region`表示用户在`region.csv`文件中自定义过滤区间；设置为`external_fasta`表示用户使用自己准备的外部 fasta 序列作为过滤数据库，请参考 Q&A。
- `filter_region_file`: 如果使用`external_region`，此处需填写过滤区间文件位置，默认为`region.csv`。
- `filter_fasta_file`: 如果使用`external_fasta`，此处需填写 fasta 文件位置，例如`/your/path/region.fasta`。
- `blast_database`：无需修改，对组装序列进行相似序列比对使用的数据库。`em_cds_pln` 使用 EBI ENA 植物编码序列数据库，`em_std_pln` 使用 EBI ENA 植物标准序列数据库。默认为`em_cds_pln`。

## 运行 Triti-Map

成功运行 `tritimap init` 并修改各配置文件后可以运行 Triti-Map 主程序。

```sh
# 启用 Triti-Map 运行环境
conda activate tritimap
# 进入分析目录
cd /your/work/path
# 实际运行程序
tritimap run -j 30 all
```

tritimap run 共支持三种分析模式，分别对应如下参数和命令：

- **only_mapping**：
  仅运行 Interval Mapping Module 进行定位分析，用于鉴定性状关联区间和突变。
  `tritimap run -j 30 only_mapping`
- **only_assembly**
  仅运行 Assembly Module 进行序列组装，用于鉴定性状关联新基因。
  `tritimap run -j 30 only_assembly`
- **all**
  同时运行 Interval Mapping Module 和 Assembly Module。
  `tritimap run -j 30 all`

### 测试运行和生成分析流程图

在实际运行之前，可以通过`-n`查看所有拟执行命令。

```
tritimap run -j 30 -n all
```

可以看到类似如下信息：

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

可以通过`--dag`生成分析流程示意图，查看命令执行顺序。

```
tritimap run -j 30 -n all --dag | dot -Tpdf > tritimap.dag.pdf
```

可以生成类似的下面的流程图。

![](https://kaopubear-1254299507.file.myqcloud.com/picgo/20210401162907.png)

通过`tritimap run --help`可以查看帮助信息。

```sh
Usage: tritimap run [OPTIONS] [[all|only_mapping|only_assembly]]
                    [SNAKEMAKE_ARGS]...

Options:
  -d, --working-dir PATH  Triti-Map running directory.
  -c, --config-file PATH  Triti-Map config file, generated with tritimap init.
  -j, --jobs INTEGER      Use at most N CPU cores/jobs in parallel.
  --profile TEXT          Name of profile to use for configuring Snakemake.
  -n, --dryrun            Do not execute anything, and display what would be done.  [default: False]

  -h, --help              Show this message and exit.
```

因为完整的分析时间可能较长，建议使用 `screen`（https://www.gnu.org/software/screen/）程序在服务器后台运行 Triti-Map。

### 查看和处理错误

在分析过程中如果出现错误，可以通过报错信息定位出现错误的分析步骤，也可以通过相应的 log 文件查找错误原因。如果分析流程在中途出现错误，修改正确的配置文件参数后可以直接再次运行`tritimap run`，已经生成结果的步骤会自动跳过，直接从尚未完成的步骤运行。

## Triti-Map 结果说明

### 目录结构

完整的 Triti-Map 结果目录如下所示：

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

每个结果目录内的内容如下：

- `01_cleandata` 原始数据经过预处理后的清洁数据
- `02_mergedata` 不同修饰类型数据合并用于后续分析的混池数据
- `03_mappingout` 参考基因组比对后的 BAM 文件
- `04_GATKout` 经过预处理后用于鉴定突变的 BAM 文件
- `05_vcfout` 变异信息原始 VCF 文件和经过过滤筛选的 VCF 文件
- `06_regionout` 用于进行 BSA 定位分析的 SNP 输入文件和结果文件
- `07_assembleout` 序列组装生成的混池特异性序列和对应功能注释文件
- `logs` 分析过程中不同步骤产生的日志文件

分析过程中可能有价值的临时文件会保存在每个目录的`temp_output`子目录中。临时文件正常情况下无需使用，程序出现错误可以结合日志文件进行检查。

### 定位模块的主要结果文件

基于 BSA 定位模块的分析结果位于目录 **06_regionout**, 以下文件名中的`xxx`和`yyy`表示 sample.csv 文件中的混池名称。

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
├── xxx_pool_vs_yyy_pool_candidateregion_1.pdf
└── xxx_pool_vs_yyy_pool_candidateregion_1.pdf
```

- `xxx_yyy_snpindex_input.txt` 用于 BSA 定位分析的 SNP 输入文件
- `xxx_yyy_qtlseqr_output.txt` QTLseqr 通过 Delta SNPindex 分析得到的 SNP 结果信息
- `xxx_yyy_qtlseqr_raw_region.txt` QTLseqr 通过 Delta SNPindex 分析得到的原始定位区间结果信息
- `xxx_yyy_qtlseqr_filter_region.txt` 通过筛选的高可信度定位区间信息
- `xxx_yyy_qtlseqr_filter_snpinfo.txt` 高可信度定位区间内的 SNP 信息
- `xxx_yyy_qtlseqr_filter_indelinfo.txt` 高可信度定位区间内的 INDEL 信息
- `xxx_yyy_qtlseqr_filter_snp.bed` 高可信度定位区间内的 SNP bed 格式文件，可以使用该文件在 Triti-Map 在线分析平台进行注释
- `xxx_yyy_qtlseqr_filter_indel.bed` 高可信度定位区间内的 INDEL bed 格式文件
- `*.pdf` pdf 文件为分析结果可视化展示，包括所有染色体 Delta SNPindex 和 snp counts 图，每个高可信度定位区间的 Delta SNPindex 图

### 序列组装模块的主要结果文件

基于序列组装的结果位于 **07_assembleout** 中，文件名中的`xxx`表示表示 sample.csv 文件中的混池名称。

```
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

- `*_merge_denovo_scaffolds.fasta`：不能和参考基因组正常匹配的混池特异性序列
- `*_candiadate_denovo.fasta`：不能和参考基因组正常匹配且和性状关联区域部分匹配的混池特异性序列
- `*_candidate_denovo2ref.info.txt`：`*_candiadate_denovo.fasta`中的序列在参考基因组中对应的位置信息
- `*_candiadate_denovo_pfam_anno.txt`：`*_candiadate_denovo.fasta`中的序列功能注释结果
- `*_candiadate_denovo_pfam_anno.fasta`：`*_candiadate_denovo.fasta`中的功能序列。
- `*_candiadate_denovo_blast_anno.txt`： `*_candiadate_denovo.fasta`中的功能序列在 EBI 数据库中的相似性注释信息
- `*_unmap_*`：表示和参考基因组完全不能匹配的序列信息，对应`*_candidater_*`

### 具体结果文件说明

**xx_yyy_snpindex_input.txt**

```
#CHROM  POSm1   POS     REF     ALT     Res_pool_ref    Res_pool_alt    Res_pool_depth  Res_pool_ratio  Sus_pool_ref    Sus_pool_alt    Sus_pool_depth  Sus_pool_ratiosnpindex
chr1A   1144541 1144542 G       A       13      147     160     0.91875 1       144     145     0.993103        -0.0743534
```

从左到右每一列数据依次表示：突变所在染色体、所在位点-1、所在位点、参考基因组碱基、突变碱基、混池 A 的参考基因组碱基数、混池 A 的 突变碱基数、混池 A 的 SNPindex、混池 B 的参考基因组碱基数、混池 B 的突变碱基数、混池 B 的 SNPindex、该位点的 Delta SNPindex。

**xxx_yyy_qtlseqr_filter_region.txt**

```
CHROM   qtl     start   end     length  nSNPs   avgSNPs_Mb      peakDeltaSNP    posPeakDeltaSNP avgDeltaSNP
chr7A   7       724111912       730119678       6007766 3262    543     0.954374003095134       729506093       0.812296971085005
```

从左到右每一列数据依次表示：定位区间所在染色体、qtlseqr 原始结果对应区间编号、区间起始位置、区间终止位置、区间长度、区间内 SNP 个数、每 1M 内的平均 SNP 个数、区间内 Delta SNPindex 峰值、区间内 Delta SNPindex 峰值所在位置、区间内平均 Delta SNPindex。具体信息可参考 QTLseqr。

**xxx_yyy_qtlseqr_filter_snpinfo.txt**

```
CHROM   POS     REF     ALT     AD_REF.Sus_pool AD_ALT.Sus_pool DP.Sus_pool     SNPindex.Sus_pool       AD_REF.Res_pool AD_ALT.Res_pool DP.Res_pool     SNPindex.Res_pool      REF_FRQ DeltaSNP        nSNPs   tricubeDeltaSNP minDP   tricubeDP       CI_99   fisher_p
chr7A   724111912       T       G       0       12      12      1       21      0       21      0       0.636363636363636       -1      445     0.505705636608414     12       20      -0.5    2.8183517084228e-09
```

从左到右每一列数据依次表示：突变所在染色体、突变所在位点、参考基因组碱基、突变碱基、混池 A 的参考基因组碱基数、混池 A 的 突变碱基数、混池 A 的深度、混池 A 的 SNPindex、混池 B 的参考基因组碱基数、混池 B 的 突变碱基数、混池 B 的深度、混池 B 的 SNPindex、、qtlseqr 计算时使用的参考等位基因频率、Delta SNPindex 、位点深度、矫正后的 Delta SNPindex、最小深度、矫正后的位点深度、99%置信区间的 SNPindex 值、fisher test 的 pvalue

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

de novo 组装序列示例 ID：`17655292:chr1A:2501389:378S202M:580:res`

序列 ID 命名规则说明：`原始组装序列ID:比对基因组序列所在染色体:比对基因组序列对应染色体位置:比对相关的CIGAR信息:序列长度:对应混池`

`*_candiadate_denovo_pfam_anno.txt`：

```
seqid   Description    E-value       Model
23558405:chr7B:733811080:14S493M9I1522M:2038:res        NB-ARC NB-ARC domain    2.3e-37       PF00931.22
```

从左到右每一列数据依次表示：序列 ID、注释功能描述、domain 比对 E-value、pfam domain ID。

**\*\_candiadate_denovo_blast_anno.txt**

```
seqid	hit_db	hit_id	hit_desc	hit_url	hsp_bit_score	hsp_align_len	hsp_identity	hsp_query_from	hsp_query_to	hsp_hit_from	hsp_hit_to
23558405:chr7B:733811080:14S493M9I1522M:2038:res	EM_CDS	AUO29721.1	Triticum urartu powdery mildew resistance protein	https://www.ebi.ac.uk/ena/data/view/AUO29721.1	3514.25	1958	99.8	81	2038	1	1958
```

从左到右每一列数据依次表示：序列 ID、比对数据库、数据库序列 ID、数据库序列描述、数据库序列信息 URL、比对分数、比对长度、比对 identity、组装序列比对起点、组装序列比对终点、数据库序列比对起点、数据库序列比对终点。

### 可视化结果说明

**全基因组 Delta SNPindex**

![](https://kaopubear-1254299507.file.myqcloud.com/picgo/20210310150108.png)

纵坐标表示`Delta SNPindex`，横坐标染色体位置。蓝色线表示 QTLseqR 计算得出的 99% 置信区间对应 Delta SNPindex 数值，灰颜色点表示对应的 SNP 低于此值，黑颜色的点表示对应的 SNP 高于此值。

**全基因组 SNP count**

![](https://kaopubear-1254299507.file.myqcloud.com/picgo/20210310150210.png)

纵坐标表示 SNPcount,横坐标表示染色体位置。

**性状关联区间局部 Delta SNPindex** （每个性状关联区间都会单独生成对应的可视化结果）

![](https://kaopubear-1254299507.file.myqcloud.com/picgo/20210401170136.png)

该图为全基因组 Delta SNPindex 的性状关联区间局部放大图。其中纵坐标表示`Delta SNPindex`，横坐标染色体位置，蓝色线表示 QTLseqR 计算得出的 99% 置信区间对应 Delta SNPindex 数值。红色部分表示该位置为计算出的性状关联区间。

## 在线注释模块

Triti-Map 的上游分析结果可以配套[在线分析平台](http://bioinfo.cemps.ac.cn/tritimap/)进行后续分析和可视化展示。包括对 SNP 进行功能和表观修饰注释，对基因进行同源基因分析，对性状关联区间进行共线性分析，对特异性片段进行功能注释等。

**在线分析平台主要功能示意图**

![](https://kaopubear-1254299507.file.myqcloud.com/picgo/20210401170258.png)

## 支持

### 常见问题

**Q：Triti-Map 可以用来分析哪些物种？**

A：Triti-Map 分析流程针对以六倍体小麦为代表的小麦族物种设计，其中多个步骤针对大基因组作物进行了优化，结果可视化方式也进行了相关调整。但是原则上 Triti-Map 可以对任何物种的混池测序数据进行 BSA 定位分析和性状相关基因的序列组装分析。

需要注意的是，Triti-Map 不适合处理含有大量 scaffolds 的基因组数据，如果基因组文件除了常规染色体外还有大量 scaffold，建议使用 `samtools` 的 `samtools faidx` 命令或 `seqkit` 的`seqkit grep`命令首先提取出染色体级别的基因组作为参考基因组文件再进行分析。

**Q：Triti-Map 可以分析哪些数据？**

A：就高通量测序数据类型而言，Triti-Map 针对 ChIP-seq 数据设计和开发，但原则上 BSA 定位模块可以分析包括全基因组重测序(WGS)、外显子组测序和转录组测序在内的大部分高通量测序数据。需要注意的是，外显子组测序数据使用序列组装模块意义不大；使用大基因组物种 WGS 数据进行序列组装模块分析需要足够的计算资源且耗时较长。

就测序材料类型而言，Triti-Map 的 BSA 定位模块目前仅支持使用 QTLseq Delta snpindex 方法进行定位区间，并未针对 EMS 诱变数据设置相应的过滤条件，所以目前更适合分析通过杂交重组获得特殊性状的材料。其中序列组装模块尤其适合类似通过不同倍性小麦杂交获得优良基因的情况，这类基因往往和参考基因组相比存在着大片段的变异或并不存在于参考基因组中。

**Q：Triti-Map 使用 ChIP-seq 数据时应该使用哪些修饰？**

A：经若干组测试，H3K27me3、H3K4me3 组合更适用于进行 BSA 定位模块分析，H3K36me3 适合进行序列组装模块分析。如果同时进行定位和新基因挖掘，建议同时使用上述三种修饰数据。

**Q：仅使用 Triti-Map 的序列组装模块时，如何使用指定区间对组装序列进行过滤？**

A：部分材料可能已知大致的定位区间但是无法在基因组中找到性状相关基因，这时可以仅使用 Triti-Map 的序列组装模块进行新基因鉴定。

如果已知的定位区间长度过小（小于 1Mb），建议适当修改起始位点和终止位点，增加区间长度，建议不小于 5Mb，然后将该区间填写在通过`tritimap init` 生成的 `region.csv` 文件中。然后，在配置文件 `config.yaml` 中将 `denovo_filter_method` 参数设置为 `external_region`，将 `filter_region_file` 设置为 `region.csv`的文件路径。

**Q：仅使用 Triti-Map 的序列组装模块时，如何使用指定数据库对组装序列进行过滤？**

A：针对多倍体小麦之间各亚基因组复杂的同源和变异特性，在仅进行序列组装模块分析时，Triti-Map 特别的引入了外部序列数据库参数来进行组装序列过滤。

通过 Triti-Map 的定位模块或其它方法已经得到高可信度性状关联区间后，可以首先通过[Triti-Map 在线分析平台](https://github.com/fei0810/Triti-Map)下载该关联区间在所有小麦族物种不同亚基因组的共线性区域序列。在配置文件 `config.yaml` 中将 `denovo_filter_method` 参数设置为 `external_fasta`，将 `filter_region_file` 设置为已经下载好的 fasta 文件路径，如 `region.fasta`。

在该模式下，组装出的混池特异性序列将会和性状关联区间的所有共线性序列进行比对，从而在更大范围内筛选可能和性状相关的新基因。

**Q: Triti-Map 对原始性状关联区间进行过滤时，参数 `filter_percentage` 应该如何设置？**

A：`filter_percentage` 参数并没有一个固定的设置标准，但是和混池样本采样的准确性、混池测序深度以及样本材料是否纯合有关。对于两个混池全部为纯合材料且测序深度足够的情况下，可以将`filter_percentage`设置为 0.75 或以上；对于非纯和材料或者测序深度较低的情况，可以将`filter_percentage`设置为 0.5。

如果因为`filter_percentage`设置不合理而导致没有定位区间，Triti-Map 会自动报错。然后可以根据`xxx_yyy_qtlseqr_raw_region.txt`中的原始结果重新设置合适的`filter_percentage`，重新运行主程序即可。

**Q：运行 Triti-Map 需要什么性能的服务器？**

过往分析经验表明，Triti-Map 分析流程中最消耗内存的是序列组装模块中的组装步骤。以六倍体小麦的 400M reads 混池 ChIP-seq 数据为例，进行序列组装模块分析时，建议服务器内存不少于 400G。

### 错误和需求提交

错误信息和需求可以通过 [GitHub issue](https://github.com/fei0810/Triti-Map/issues) 进行提交。

开发者联系方式 zhaofei920810@gmail.com

### 更新日志

### 引用

def get_merged_fastq(wildcards):
    return expand(dir_path+"/02_mergedata/{bulk}_{bulktype}_merge_fastp_R{group}.fq.gz", bulk=wildcards.bulk, bulktype=wildcards.bulktype, group=[1, 2])

rule STAR1Mapping:
    input:
        get_merged_fastq
    params:
        sharemem = "NoSharedMemory",
        stardir = config['ref']['STARdir'],
        bulk  = lambda wildcards: wildcards.bulk,
        bulktype = lambda wildcards: wildcards.bulktype
    output:
        dir = directory(dir_path+"/03_mappingout/{bulk}_{bulktype}_step1/"),
        sj = dir_path+"/03_mappingout/{bulk}_{bulktype}_step1/{bulk}_{bulktype}SJ.out.tab",
    threads: thread
    message: "\nMapping {input} with STAR (step 1)\n"
    log:
        dir_path+"/logs/{bulk}_{bulktype}_STARstep1.log"
    shell: """
    STAR --genomeLoad {params.sharemem} --runThreadN {threads} \
    --genomeDir {params.stardir} --readFilesCommand zcat \
    --readFilesIn {input} \
    --outFileNamePrefix {output.dir}/{params.bulk}_{params.bulktype} \
    --outSAMtype BAM Unsorted --limitBAMsortRAM 10000000000 > {log} 2>&1
    """

rule STAR2Mapping:
    input:
        fq = get_merged_fastq,
        sj = dir_path+"/03_mappingout/{bulk}_{bulktype}_step1/{bulk}_{bulktype}SJ.out.tab"
    params:
        stardir = config['ref']['STARdir'],
        sharemem = "NoSharedMemory",
        bulk  = lambda wildcards: wildcards.bulk,
        bulktype = lambda wildcards: wildcards.bulktype
    output:
        dir = directory(dir_path+"/03_mappingout/{bulk}_{bulktype}_step2"),
        step2bam = dir_path+"/03_mappingout/{bulk}_{bulktype}_step2/{bulk}_{bulktype}Aligned.sortedByCoord.out.bam"
    threads: thread
    message: "\nMapping {input.fq} with STAR (Step 2)\n"
    log:
        dir_path+"/logs/{bulk}_{bulktype}_STARstep2.log"
    shell: """
    STAR --genomeLoad {params.sharemem} --runThreadN {threads} \
    --genomeDir {params.stardir} --readFilesCommand zcat \
    --readFilesIn {input.fq} \
    --outFileNamePrefix {output.dir}/{params.bulk}_{params.bulktype} \
    --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 \
    --sjdbFileChrStartEnd {input.sj} \
    --outReadsUnmapped Fastx > {log} 2>&1 ;
    samtools index -c {output.step2bam}
    """

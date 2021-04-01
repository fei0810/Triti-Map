def get_clean_fastq(wildcards):
    return expand(dir_path+"/01_cleandata/{s.sample}_{s.marker}_{s.bulk}_{s.type}_fastp_R{group}.fq.gz",s=samples.itertuples(),group=[1,2])

rule MergeCleanPoolFastq:
    input:
       get_clean_fastq
    params:
        bulk1 = bulkname[0],
        bulk2 = bulkname[1],
        dir_path = dir_path
    output:
        bulk1_fq1 = join(dir_path+"/02_mergedata", bulkname[0] + "_pool_merge_fastp_R1.fq.gz"),
        bulk1_fq2 = join(dir_path+"/02_mergedata", bulkname[0] + "_pool_merge_fastp_R2.fq.gz"),
        bulk2_fq1 = join(dir_path+"/02_mergedata", bulkname[1] + "_pool_merge_fastp_R1.fq.gz"),
        bulk2_fq2 = join(dir_path+"/02_mergedata", bulkname[1] + "_pool_merge_fastp_R2.fq.gz")
    message: "\nMerge different type of ChIP-seq data to one file to calling snp. Input fils: {input}\n"
    run:
        shell("""
        ls {params.dir_path}/01_cleandata/*{params.bulk1}_pool_fastp_R1.fq.gz | xargs cat > {output.bulk1_fq1} ; \
        ls {params.dir_path}/01_cleandata/*{params.bulk1}_pool_fastp_R2.fq.gz | xargs cat > {output.bulk1_fq2} ; \
        ls {params.dir_path}/01_cleandata/*{params.bulk2}_pool_fastp_R1.fq.gz | xargs cat > {output.bulk2_fq1} ; \
        ls {params.dir_path}/01_cleandata/*{params.bulk2}_pool_fastp_R2.fq.gz | xargs cat > {output.bulk2_fq2}
        """)

rule MergeCleanParentFastq:
    input:
        get_clean_fastq
    params:
        bulk1 = bulkname[0],
        bulk2 = bulkname[1],
        dir_path = dir_path
        # prefix  = lambda wildcards: wildcards.sample
    output:
        bulk1_fq1 = join(dir_path+"/02_mergedata", bulkname[0] + "_parent_merge_fastp_R1.fq.gz"),
        bulk1_fq2 = join(dir_path+"/02_mergedata", bulkname[0] + "_parent_merge_fastp_R2.fq.gz"),
        bulk2_fq1 = join(dir_path+"/02_mergedata", bulkname[1] + "_parent_merge_fastp_R1.fq.gz"),
        bulk2_fq2 = join(dir_path+"/02_mergedata", bulkname[1] + "_parent_merge_fastp_R2.fq.gz")
    message: "\nMerge different type of ChIP-seq data to one file to calling snp. Input file: {input}\n"
    run:
        shell("""
        ls {params.dir_path}/01_cleandata/*{params.bulk1}_parent_fastp_R1.fq.gz | xargs cat > {output.bulk1_fq1} ; \
        ls {params.dir_path}/01_cleandata/*{params.bulk1}_parent_fastp_R2.fq.gz | xargs cat > {output.bulk1_fq2} ; \
        ls {params.dir_path}/01_cleandata/*{params.bulk2}_parent_fastp_R1.fq.gz | xargs cat > {output.bulk2_fq1} ; \
        ls {params.dir_path}/01_cleandata/*{params.bulk2}_parent_fastp_R2.fq.gz | xargs cat > {output.bulk2_fq2}
        """)
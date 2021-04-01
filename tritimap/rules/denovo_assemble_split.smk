def get_pool_R1_trimmed(wildcards):
    return expand(dir_path+"/01_cleandata/{poolsample}_{poolmarker}_{poolbulk}_{pooltype}_fastp_R1.fq.gz", poolsample=wildcards.poolsample, poolmarker=samples.loc[wildcards.poolsample].marker,poolbulk=samples.loc[wildcards.poolsample].bulk,pooltype=samples.loc[wildcards.poolsample].type)
def get_pool_R2_trimmed(wildcards):
    return expand(dir_path+"/01_cleandata/{poolsample}_{poolmarker}_{poolbulk}_{pooltype}_fastp_R2.fq.gz", poolsample=wildcards.poolsample, poolmarker=samples.loc[wildcards.poolsample].marker,poolbulk=samples.loc[wildcards.poolsample].bulk,pooltype=samples.loc[wildcards.poolsample].type)

rule AssembleSplit:
    input:
        fq1 = get_pool_R1_trimmed,
        fq2 = get_pool_R2_trimmed
    output:
        dir_path+"/07_assembleout/{poolsample}_{poolmarker}_{poolbulk}_{pooltype}_split_denovo_scaffolds.fasta"
    params:
        memory = config['memory'],
        datatype = config['datatype'],
        dir=dir_path,
        scriptdir = script_dir
    threads: thread
    message: "\nAssembly {input} with trinity\n"
    log:
        dir_path+'/logs/{poolsample}_{poolmarker}_{poolbulk}_{pooltype}_split_denovo_scaffolds.log'
    shell:"""
    set +e
    bash {params.scriptdir}/denovo_assemble.sh {input.fq1} {input.fq2} {params.datatype} {threads} {params.memory} {params.dir} {output} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """

rule MergeSplitAssemble:
    input:
        expand(dir_path+"/07_assembleout/{s.sample}_{s.marker}_{s.bulk}_{s.type}_split_denovo_scaffolds.fasta",s=samples[samples.type=='pool'].itertuples())
    params:
        bulk1 = bulkname[0],
        bulk2 = bulkname[1],
        dir_path = dir_path
    output:
        bulk1 = join(dir_path+"/07_assembleout", bulkname[0] + "_merge_denovo_scaffolds.fasta"),
        bulk2 = join(dir_path+"/07_assembleout", bulkname[1] + "_merge_denovo_scaffolds.fasta")
    message: "\nMerge different type of ChIP-seq data to one file. Input files: {input}\n"
    run:
        shell("""
        ls {params.dir_path}/07_assembleout/*_{params.bulk1}_pool_split_denovo_scaffolds.fasta | xargs cat > {output.bulk1} ; \
        ls {params.dir_path}/07_assembleout/*_{params.bulk2}_pool_split_denovo_scaffolds.fasta | xargs cat > {output.bulk2}
    """)
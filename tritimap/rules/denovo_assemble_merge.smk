def get_R1_trimmed(wildcards):
    return expand(dir_path+"/02_mergedata/{bulk}_pool_merge_fastp_R1.fq.gz", bulk=wildcards.bulk)
def get_R2_trimmed(wildcards):
    return expand(dir_path+"/02_mergedata/{bulk}_pool_merge_fastp_R2.fq.gz", bulk=wildcards.bulk)

rule AssembleMerge:
    input:
        bulk1_fq1 = get_R1_trimmed,
        bulk1_fq2 = get_R2_trimmed
    output:
        dir_path+"/07_assembleout/{bulk}_merge_denovo_scaffolds.fasta"
    params:
        memory = config['memory'],
        datatype = config['datatype'],
        dir=dir_path,
        scriptdir = script_dir
    threads: thread
    message: "\nAssembly {input} \n"
    log:
        dir_path+'/logs/{bulk}_pool_assemble.log'
    shell:"""
        set +e
        bash {params.scriptdir}/denovo_assemble.sh {input.bulk1_fq1} {input.bulk1_fq2} {params.datatype} {threads} {params.memory} {params.dir} {output} > {log} 2>&1
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
    """
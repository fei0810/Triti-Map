def get_all_scaffolds(wildcards):
    return expand(dir_path+"/07_assembleout/{bulk}_merge_denovo_scaffolds.fasta", **wildcards)
rule uniqScaffoldByFasta:
    input:
        fa=expand(dir_path+"/07_assembleout/{bulk}_merge_denovo_scaffolds.fasta", bulk=bulkname),
    params:
        length = config['scaffold']['minlength'],
        memory = config['memory'],
        datatype = config['datatype'],
        genome = config['ref']['genome'],
        runmodule = config['module'],
        database = config['filter_fasta_file'],
        scriptdir = script_dir
    output:
        scaffold1 = join(dir_path+"/07_assembleout", bulkname[0] + "_candidate_denovo.fasta"),
        scaffold2 = join(dir_path+"/07_assembleout", bulkname[1] + "_candidate_denovo.fasta"),
        unmap1 = join(dir_path+"/07_assembleout", bulkname[0] + "_unmap_denovo.fasta"),
        unmap2 = join(dir_path+"/07_assembleout", bulkname[1] + "_unmap_denovo.fasta")
    message: "\nGet bulk uniq scaffolds by fasta file. Input file: {input.fa}\n"
    threads: thread
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates()) + "_uniqscaffold.log")
    shell:"""
    set +e
    bash {params.scriptdir}/getuniqscaffold_by_fasta.sh {input.fa} {params.database} {params.length} {params.genome} {threads} {params.datatype} {params.memory} {output.scaffold1} {output.scaffold2} {output.unmap1} {output.unmap2} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """

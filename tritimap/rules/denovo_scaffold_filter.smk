def get_all_scaffolds(wildcards):
    return expand(dir_path+"/07_assembleout/{bulk}_merge_denovo_scaffolds.fasta", **wildcards)
rule uniqScaffold:
    input:
        fa=expand(dir_path+"/07_assembleout/{bulk}_merge_denovo_scaffolds.fasta", bulk=bulkname),
        region=join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_filter_region.txt")
    params:
        length = config['scaffold']['minlength'],
        memory = config['memory'],
        datatype = config['datatype'],
        genome = config['ref']['genome'],
        runmodule = config['module'],
        scriptdir = script_dir,
        blastfilter = config['scaffold_uniq_percentage']
    output:
        scaffold1 = join(dir_path+"/07_assembleout", bulkname[0] + "_candidate_denovo.fasta"),
        scaffold2 = join(dir_path+"/07_assembleout", bulkname[1] + "_candidate_denovo.fasta"),
        unmap1 = join(dir_path+"/07_assembleout", bulkname[0] + "_unmap_denovo.fasta"),
        unmap2 = join(dir_path+"/07_assembleout", bulkname[1] + "_unmap_denovo.fasta"),
        scaffold1summary = join(dir_path+"/07_assembleout", bulkname[0] + "_candidate_denovo2ref.info.txt"),
        scaffold2summary = join(dir_path+"/07_assembleout", bulkname[1] + "_candidate_denovo2ref.info.txt")
    message: "\nGet uniq scaffolds by Triti-Map. Input file: {input.fa}\n"
    threads: thread
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates()) + "_uniqscaffold.log")
    shell:"""
    set +e
    bash {params.scriptdir}/getuniqscaffold.sh {input.fa} {input.region} {params.length} {params.genome} {threads} {params.datatype} {params.memory} {output.scaffold1} {output.scaffold2} {output.unmap1} {output.unmap2} {output.scaffold1summary} {output.scaffold2summary} {params.blastfilter} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """
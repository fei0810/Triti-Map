def get_all_scaffolds(wildcards):
    return expand(dir_path+"/07_assembleout/{bulk}_merge_denovo_scaffolds.fasta", **wildcards)
rule uniqScaffoldByRegion:
    input:
        fa = expand(dir_path+"/07_assembleout/{bulk}_merge_denovo_scaffolds.fasta", bulk=bulkname),
        region = config['filter_region_file']
    params:
        length = config['scaffold']['minlength'],
        memory = config['memory'],
        datatype = config['datatype'],
        genome = config['ref']['genome'],
        scriptdir = script_dir
    output:
        region = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_filter_region.txt"),
        scaffold1 = join(dir_path+"/07_assembleout", bulkname[0] + "_candidate_denovo.fasta"),
        scaffold2 = join(dir_path+"/07_assembleout", bulkname[1] + "_candidate_denovo.fasta"),
        unmap1 = join(dir_path+"/07_assembleout", bulkname[0] + "_unmap_denovo.fasta"),
        unmap2 = join(dir_path+"/07_assembleout", bulkname[1] + "_unmap_denovo.fasta")
    message: "\nget uniq scaffolds\n"
    threads: thread
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates()) + "_uniqscaffold.log")
    shell:"""
    set +e
    cat {input.region} | tr ',' '\t' | awk '{{print $1"\tqtl\t"$2"\t"$3}}' > {output.region}
    bash {params.scriptdir}/getuniqscaffold.sh {input.fa} {output.region} {params.length} {params.genome} {threads} {params.datatype} {params.memory} {output.scaffold1} {output.scaffold2} {output.unmap1} {output.unmap2} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """
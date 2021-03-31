rule snpindexRegion:
    input:
        join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates()) + "_snpindex_input.txt")
    params:
        depth = config['snpindex']['min_avgSNPcount_Mb'],
        index = config['snpindex']['min_avgindex'],
        genome = config['ref']['genome'],
        scriptdir = script_dir
    output:
        region = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_snpindex_region.txt"),
        snp = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_snpindex_candicatesnp.txt"),
    message: "\nSNPindex calculation\n"
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates())+ "_snpindexregion.log")
    shell:"""
    set +e
    bash {params.scriptdir}/snpindex2region.sh {input} {output.region} {output.snp} {params.depth} {params.index} {params.genome} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """

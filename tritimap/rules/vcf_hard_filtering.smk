rule vcfHardFiltration:
    input:
        join(dir_path+"/05_vcfout", "_".join(samples.bulk.drop_duplicates()) + "_raw_gatk.vcf")
    params:
        genome = config['ref']['genome'],
        scriptdir = script_dir
    output:
        join(dir_path+"/05_vcfout", "_".join(samples.bulk.drop_duplicates()) + "_hardfiltered_gatk.vcf")
    message:"\nRaw vcf file hard filter\n"
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates()) + "_gatkhardfilter.log")
    shell:"""
    set +e
    bash {params.scriptdir}/gatkhardfilter.sh \
    {input} {output} {params.genome} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """

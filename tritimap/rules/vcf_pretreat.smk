rule vcf2tab:
    input:
        join(dir_path+"/05_vcfout", "_".join(samples.bulk.drop_duplicates()) + "_genofiltered_gatk.vcf")
    params:
        pool1 = pool1Name,
        pool2 = pool2Name,
        genome = config['ref']['genome'],
        snpdepth = config['gatk']['min_SNP_DP'],
        scriptdir = script_dir
    output:
        snpindex = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates()) + "_snpindex_input.txt"),
        qtlseqr = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates()) + "_qtlseqr_input.txt")
    message: "\nTransform genotype filterd VCF file to tab format\n"
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates()) + "_vcf2tab.log")
    shell:"""
    set +e
    bash {params.scriptdir}/vcf2tab.sh {input} {output.snpindex} {output.qtlseqr} {params.pool1} {params.pool2} {params.snpdepth} {params.genome} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """
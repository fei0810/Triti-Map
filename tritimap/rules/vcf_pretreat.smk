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
        indel = temp(join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates()) + "_indel.txt")),
        qtlseqr = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates()) + "_qtlseqr_input.txt")
    message: "\nTransform genotype filterd VCF file to tab format. Input file: {input}\n"
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates()) + "_vcf2tab.log")
    shell:"""
    set +e
    bash {params.scriptdir}/vcf2tab.sh {input} {output.snpindex} {output.qtlseqr} {params.pool1} {params.pool2} {params.snpdepth} {params.genome} {output.indel} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """

rule getRegionIndelandSNPbed:
    input:
        indel = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates()) + "_indel.txt"),
        snp = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_filter_snpinfo.txt"),
        region=join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_filter_region.txt")
    params:
        dirpath = join(dir_path+"/06_regionout")
    output:
        indelinfo = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_filter_indelinfo.txt"),
        indelbed = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_filter_indel.bed"),
        snpbed = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_filter_snp.bed")
    message: "\nGet candidate region indel and snp bed file. Input file: {input.indel}; {input.snp}; {input.region}\n"
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates()) + "_get_region_indel_snp_bed.log")
    shell:"""
    cat {input.region} | grep -v 'start' | cut -f1,3,4 | sort -k1,1 -k2,2n > {params.dirpath}/temp.candidateregion.indel.bed
    bedmap --skip-unmapped --delim '\t' --header --echo {input.indel} {params.dirpath}/temp.candidateregion.indel.bed | cat <(head -n1 {input.indel}) - | cut -f1,3- | sed 's/^#CHROM/CHROM/' > {output.indelinfo} && rm {params.dirpath}/temp.candidateregion.indel.bed
    cat {output.indelinfo} | grep -v '^CHROM'| awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2,$3,$4}}' |sort -k1,1 -k2,2n > {output.indelbed}
    cat {input.snp} | grep -v '^CHROM'| awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2,$3,$4}}' |sort -k1,1 -k2,2n > {output.snpbed}
    """
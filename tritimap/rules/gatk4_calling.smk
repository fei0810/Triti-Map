def get_all_bam(wildcards):
    return expand(dir_path+"/04_GATKout/{bulk}_{bulktype}_final.bam", **wildcards)

contigs = pd.read_csv(config["ref"]["genome"] + ".fai", sep = '\t', header=None, usecols=[0], squeeze=True, dtype=str)

# if len(contigs) <= 23:
#     print("Your genome file have " +  str(len(contigs)) + " contigs.")
# else:
#     print("Your genome file have " + str(len(contigs)) + " contigs.\nYou may need modify your genome file and rebuild index")

rule gatk4HC:
    input:
        expand(dir_path+"/04_GATKout/{bulk}_{bulktype}_final.bam", bulk=bulkname, bulktype=typename)
    params:
        filter = "--dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20.0 --create-output-variant-index false",
        genome = config['ref']['genome']
    output:
        temp(join(dir_path+"/05_vcfout", "_".join(samples.bulk.drop_duplicates()) + "_{contig}.vcf"))
    message: "\nCalling {input} variants with GATK4 HaplotyeCaller by contig\n"
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates()) + "_{contig}.log")
    shell:"""
    gatk HaplotypeCaller \
    -L {wildcards.contig} -R {params.genome} \
    $(for bam in {input}; do echo "-I $bam"; done) {params.filter} \
    -O {output} > {log} 2>&1
    """

rule GATK4MergeVcf:
    input:
        expand(join(dir_path+"/05_vcfout", "_".join(samples.bulk.drop_duplicates()) + "_{contig}.vcf"), contig=contigs)
    params:
        java_p = java_parameter
    output:
        join(dir_path+"/05_vcfout", "_".join(samples.bulk.drop_duplicates()) + "_raw_gatk.vcf")
    message: "\nmerge all separate VCF files\n"
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates()) + "gatkmergevcf.log")
    shell:"""
    gatk {params.java_p} MergeVcfs \
    $(for vcf in {input}; do echo "-I $vcf"; done) -O {output} > {log} 2>&1
    """

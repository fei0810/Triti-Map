def get_merged_fastq(wildcards):
    return expand(dir_path+"/02_mergedata/{bulk}_{bulktype}_merge_fastp_R{group}.fq.gz", bulk=wildcards.bulk, bulktype=wildcards.bulktype,group=[1, 2])

rule dnaBWAmem2Mapping:
    input:
        get_merged_fastq
    params:
        group = r"@RG\tID:{bulk}_{bulktype}\tSM:{bulk}_{bulktype}\tPL:illumina\tLB:lib1\tPU:{bulk}_{bulktype}",
        genome = config['ref']['genome'],
    output:
        raw = dir_path+"/03_mappingout/{bulk}_{bulktype}.sort.bam",
        truemap = dir_path+"/03_mappingout/{bulk}_{bulktype}.truemap.sort.bam",
        errormap = dir_path+"/03_mappingout/{bulk}_{bulktype}.errormap.sort.bam",
        unmap = dir_path+"/03_mappingout/{bulk}_{bulktype}.unmap.bam"
    threads: thread
    message: "mapping {input} with BWA mem2"
    run:
        shell("""
        bwa-mem2 mem -v 1 -t {threads} -M -Y -R '{params.group}' {params.genome} {input} | \
        samtools sort -@ {threads} -o {output.raw} - ;
        samtools view -h -f 3 {output.raw} | egrep -v 'XA:Z|SA:Z' | samtools sort -@ {threads} -o {output.truemap} ;
        samtools view -h {output.raw} | egrep '^@|XA:Z|SA:Z' | samtools sort -@ {threads} -o {output.errormap} ;
        samtools view -f 4 -b {output.raw} > {output.unmap} ;
        samtools index -c {output.raw};
        samtools index -c {output.truemap};
        samtools index -c {output.errormap};
    """)
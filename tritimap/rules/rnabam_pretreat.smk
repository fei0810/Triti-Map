def get_rnabam(wildcards):
    return expand(dir_path+"/03_mappingout/{bulk}_{bulktype}_step2/{bulk}_{bulktype}Aligned.sortedByCoord.out.bam", **wildcards)

rule rnaGATK4ReplaceRG:
    input:
        get_rnabam
    params:
        java_p = java_parameter,
        info = "-LB lib1 -PL illumina -PU {bulk}_{bulktype} -SM {bulk}_{bulktype} -ID {bulk}_{bulktype} -SO coordinate --CREATE_INDEX false"
    output:
        temp(dir_path+"/04_GATKout/{bulk}_{bulktype}_reprg.bam")
    log:
        dir_path+"/logs/step3_{bulk}_{bulktype}_rna_gatkreprg.log"
    message: "Replace ReadsGroups {input} with GATK4"
    shell:("gatk {params.java_p} AddOrReplaceReadGroups -I {input} -O {output} {params.info} > {log} 2>&1")

rule rnaGATK4Markdup:
    input:
        dir_path+"/04_GATKout/{bulk}_{bulktype}_reprg.bam"
    params:
        java_p = java_parameter,
        info = "--REMOVE_DUPLICATES true --USE_JDK_DEFLATER true --USE_JDK_INFLATER true --CREATE_INDEX false --VALIDATION_STRINGENCY SILENT"
    output:
        bam = temp(dir_path+"/04_GATKout/{bulk}_{bulktype}_rmdup.bam"),
        metrics = dir_path+"/04_GATKout/{bulk}_{bulktype}.rmdup.metrics"
    log:
        dir_path+"/logs/{bulk}_{bulktype}_rna_gatkrmdup.log"
    message: "\nRemove Duplicates {input} with GATK4\n"
    shell:("samtools index -c {input} && gatk {params.java_p} MarkDuplicates {params.info} -I {input} -O {output.bam} -M {output.metrics} > {log} 2>&1")

rule rnaFilter2Uniqmap:
    input:
        dir_path+"/04_GATKout/{bulk}_{bulktype}_rmdup.bam"
    output:
        temp(dir_path+"/04_GATKout/{bulk}_{bulktype}_uniqmap.bam")
    log:
        dir_path+"/logs/{bulk}_{bulktype}_rna_uniqmap.log"
    message: "\nFilter {input} to uniqmap and prepore pair with sambamba\n"
    shell:("samtools index -c {input} && samtools view -h {input} |egrep 'NH:i:1[^0-9]|^@' | samtools view -h -f 3 -S -b - | samtools sort -o {output} - ")

rule rnaGATK4FixMateInfo:
    input:
        dir_path+"/04_GATKout/{bulk}_{bulktype}_uniqmap.bam"
    params:
        java_p = java_parameter,
        info = "-SO coordinate --CREATE_INDEX false"
    output:
        temp(dir_path+"/04_GATKout/{bulk}_{bulktype}_fix.bam")
    log:
        dir_path+"/logs/{bulk}_{bulktype}_rna_gatkfix.log"
    message: "\nFix Mate Information {input} with GATK4\n"
    shell:("samtools index -c {input} && gatk {params.java_p} FixMateInformation {params.info} -I {input} -O {output} > {log} 2>&1")

rule rnaGATK4SplitNCigar:
    input:
        dir_path+"/04_GATKout/{bulk}_{bulktype}_fix.bam"
    params:
        java_p = java_parameter,
        genome = config['ref']['genome']
    output:
        dir_path+"/04_GATKout/{bulk}_{bulktype}_final.bam"
    log:
        dir_path+"/logs/{bulk}_{bulktype}_rna_gatksplit.log"
    message: "\nSplit N Cigar Reads {input} with GATK4\n"
    shell:("samtools index -c {input} && gatk {params.java_p} SplitNCigarReads --create-output-bam-index false -R {params.genome} -I {input} -O {output} && samtools index -c {output}")

rule rnaErrorMapReads:
    input:
        dir_path+"/04_GATKout/{bulk}_{bulktype}_final.bam"
    output:
        truemap = dir_path+"/04_GATKout/{bulk}_{bulktype}.truemap.sort.bam",
        errormap = dir_path+"/04_GATKout/{bulk}_{bulktype}.errormap.sort.bam"
    threads: thread
    message: "Get {input} errormap(XA or SA tag) reads and truemap reads"
    run:
        shell("""
        samtools view -h -f 3 {input} | egrep -v 'XA:Z|SA:Z' | samtools sort -@ {threads} -o {output.truemap} ;
        samtools view -h {input} | egrep '^@|XA:Z|SA:Z' | samtools sort -@ {threads} -o {output.errormap} ;
        samtools index -c {output.truemap};
        samtools index -c {output.errormap}
    """)
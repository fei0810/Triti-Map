rule dnaGATK4Markdup:
    input:
        dir_path+"/03_mappingout/{bulk}_{bulktype}.truemap.sort.bam"
    params:
        java_p = java_parameter,
        info = "--REMOVE_DUPLICATES true --USE_JDK_DEFLATER true --USE_JDK_INFLATER true --CREATE_INDEX false --VALIDATION_STRINGENCY SILENT"
    output:
        bam = dir_path+"/04_GATKout/{bulk}_{bulktype}_final.bam",
        metrics = dir_path+"/04_GATKout/{bulk}_{bulktype}.final.metrics"
    log:
        dir_path+"/logs/{bulk}_{bulktype}_dna_gatkrmdup.log"
    message: "\nRemove Duplicates {input} with GATK4\n"
    shell:("samtools index -c {input} && gatk {params.java_p} MarkDuplicates {params.info} -I {input} -O {output.bam} -M {output.metrics} > {log} 2>&1 && samtools index -c {output.bam}")
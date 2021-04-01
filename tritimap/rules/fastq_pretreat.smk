def get_fastq(wildcards):
    fastqs = samples.loc[(wildcards.sample, wildcards.marker, wildcards.bulk, wildcards.bulktype), ["fq1", "fq2"]].dropna()
    return {"fq1": fastqs.fq1, "fq2": fastqs.fq2}

rule Rawfastq2Cleanfastq:
    input:
        unpack(get_fastq)
    params:
        quality = config['trim']['quality'],
        minlen = config['trim']['minlength'],
    output:
        fq1 = dir_path+"/01_cleandata/{sample}_{marker}_{bulk}_{bulktype}_fastp_R1.fq.gz",
        fq2 = dir_path+"/01_cleandata/{sample}_{marker}_{bulk}_{bulktype}_fastp_R2.fq.gz",
        json = dir_path+"/01_cleandata/{sample}_{marker}_{bulk}_{bulktype}.json",
        html = dir_path+"/01_cleandata/{sample}_{marker}_{bulk}_{bulktype}.html"
    log:
        dir_path+"/logs/{sample}_{marker}_{bulk}_{bulktype}_fastp.log"
    threads: thread
    message: "\nProcessing {input.fq1} and {input.fq2} with fastp\n"
    run:
        shell("""
        fastp -w {threads} \
        -j {output.json} \
        -h {output.html} \
        -3 -l {params.minlen} -q {params.quality} \
        -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} > {log} 2>&1 """)

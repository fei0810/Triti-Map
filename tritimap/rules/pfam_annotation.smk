def get_candidatefasta(wildcards):
    return expand(dir_path+"/07_assembleout/{bulk}_candidate_denovo.fasta", **wildcards)
def get_unmapfasta(wildcards):
    return expand(dir_path+"/07_assembleout/{bulk}_unmap_denovo.fasta", **wildcards)

rule PfamAnnoEBI:
    input:
        get_candidatefasta
    params:
        email = config['email'],
        scriptdir = script_dir
    output:
        pfamanno = dir_path+"/07_assembleout/{bulk}_candidate_denovo_pfam_anno.txt",
        pfamseq = dir_path+"/07_assembleout/{bulk}_candidate_denovo_pfam_anno.fasta"
    message: "\nUse EBI interproscan API do annotation: {input}\n"
    log:
        join(dir_path+"/logs", "{bulk}_PfamAnnoEBI.log")
    shell:"""
    set +e
    bash {params.scriptdir}/format_pfamout.sh {input} {params.email} {params.scriptdir}/hmmer3_hmmscan.py {output.pfamanno} {output.pfamseq} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """

rule unmap_PfamAnnoEBI:
    input:
        get_unmapfasta
    params:
        email = config['email'],
        scriptdir = script_dir
    output:
        pfamanno = dir_path+"/07_assembleout/{bulk}_unmap_denovo_pfam_anno.txt",
        pfamseq = dir_path+"/07_assembleout/{bulk}_unmap_denovo_pfam_anno.fasta"
    message: "\nUse EBI interproscan API do annotation: {input}\n"
    log:
        join(dir_path+"/logs", "{bulk}_unmap_PfamAnnoEBI.log")
    shell:"""
    set +e
    bash {params.scriptdir}/format_pfamout.sh {input} {params.email} {params.scriptdir}/hmmer3_hmmscan.py {output.pfamanno} {output.pfamseq} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """
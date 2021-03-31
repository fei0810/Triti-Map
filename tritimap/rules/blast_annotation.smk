def get_pfamcandidatefasta(wildcards):
    return expand(dir_path+"/07_assembleout/{bulk}_candidate_denovo_pfam_anno.fasta", **wildcards)
def get_pfam_unmap_candidatefasta(wildcards):
    return expand(dir_path+"/07_assembleout/{bulk}_unmap_denovo_pfam_anno.fasta", **wildcards)

rule BlastAnnoEBI:
    input:
        get_pfamcandidatefasta
    params:
        email = config['email'],
        database = config['blast_database'],
        scriptdir = script_dir
    output:
        dir_path+"/07_assembleout/{bulk}_candidate_denovo_blast_anno.txt"
    message: "\nUse EBI Blast API do annotation\n"
    log:
        join(dir_path+"/logs", "{bulk}_BlastAnnoEBI.log")
    shell:"""
    set +e
    bash {params.scriptdir}/format_blastout.sh {input} {params.email} scripts/ncbiblast.py scripts/json2tsv.R {params.database} {output} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """

rule unmap_BlastAnnoEBI:
    input:
        get_pfam_unmap_candidatefasta
    params:
        email = config['email'],
        database = config['blast_database'],
        scriptdir = script_dir
    output:
        dir_path+"/07_assembleout/{bulk}_unmap_denovo_blast_anno.txt"
    message: "\nUse EBI Blast API do annotation\n"
    log:
        join(dir_path+"/logs", "{bulk}_unmap_BlastAnnoEBI.log")
    shell:"""
    set +e
    bash {params.scriptdir}/format_blastout.sh {input} {params.email} {params.scriptdir}/ncbiblast.py {params.scriptdir}/json2tsv.R {params.database} {output} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """
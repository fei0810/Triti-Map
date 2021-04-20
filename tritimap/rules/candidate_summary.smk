rule CandidateSummary:
    input:
        blastinfo = dir_path+"/07_assembleout/{bulk}_candidate_denovo_blast_anno.txt",
        pfaminfo = dir_path+"/07_assembleout/{bulk}_candidate_denovo_blast_anno.txt",
        mappinginfo = dir_path+"/07_assembleout/{bulk}_candidate_denovo2ref.info.txt"
    params:
        scriptdir = script_dir
    output:
        dir_path+"/07_assembleout/{bulk}_candidate_denovo_summary_info.txt"
    log:
        dir_path+"/logs/{bulk}_candidate_denovo_summary_info.log"
    threads: thread
    message: "\nProcessing {input.blastinfo} {input.pfaminfo} {input.mappinginfo} get summary information.\n"
    shell:"""
    set +e
    bash {params.scriptdir}/candidate_summary.sh {input.blastinfo} {input.pfaminfo} {input.mappinginfo} {output} > {log} 2>&1
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
        exit 1
    else
        exit 0
    fi
    """
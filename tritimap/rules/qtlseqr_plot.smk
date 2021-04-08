rule QTLseqr_plot:
    input:
        qtlseqr = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates()) + "_qtlseqr_input.txt"),
        snpindex = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates()) + "_snpindex_input.txt")
    params:
        bulk1 = pool1Name,
        bulk2 = pool2Name,
        pop_struc = config['snpindex']['pop_struc'],
        bulksize = config['snpindex']['bulksize'],
        winsize = config['snpindex']['winsize'],
        filterprobs = config['snpindex']['filter_probs'],
        pvalue = config['snpindex']['fisher_p'],
        min_length = config['snpindex']['min_length'],
        dirname = dir_path
    output:
        qtlout = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_output.txt"),
        qtlregion = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_filter_region.txt"),
        qtlrawregion = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_raw_region.txt"),
        qtlsnpinfo = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_filter_snpinfo.txt"),
        qtlcount_p = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_SNPcounts_point.pdf"),
        qtlcount_l = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_SNPcounts_line.pdf"),
        qtlindex_p = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_SNPindex_point.pdf"),
        qtlindex_l = join(dir_path+"/06_regionout", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr_SNPindex_line.pdf")
    message: "\nGet candidate regions and mutations by QTLseqr. Input files: {input.qtlseqr}; {input.snpindex}\n"
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates())+ "_qtlseqr.log")
    script:
         "../scripts/qtlseqr_plot.R"

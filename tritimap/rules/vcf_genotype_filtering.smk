rule vcfGenotypeFiltration:
    input:
        join(dir_path+"/05_vcfout", "_".join(samples.bulk.drop_duplicates()) + "_hardfiltered_gatk.vcf")
    params:
        filtertype = genotypefilterType,
        pool1 = pool1Name,
        pool2 = pool2Name,
        repool = rece_poolName,
        parent1 = parent1Name,
        parent2 = parent2Name,
        reparent = rece_parentName,
        golobe = java_parameter,
        genome = config['ref']['genome'],
    output:
        join(dir_path+"/05_vcfout", "_".join(samples.bulk.drop_duplicates()) + "_genofiltered_gatk.vcf")
    message:"\nRaw vcf file genotype filter. Input file: {input}\n"
    log:
        join(dir_path+"/logs", "_".join(samples.bulk.drop_duplicates()) + "_genofilter.log")
    run:
        if params.filtertype == 'onlypool':
        # homozygous and identical mutations shared by two mixed pools
            shell("""
            gatk {params.golobe} SelectVariants \
            -R {params.genome} \
            -V {input} \
            -O {output} \
            --invertSelect \
            --selectExpressions '(vc.getGenotype("{params.pool1}").isHomVar() && vc.getGenotype("{params.pool2}").isHomVar()) || (vc.getGenotype("{params.pool1}").isHomRef() && vc.getGenotype("{params.pool2}").isHomRef())'
            """)
        #TODO HomRef() or not
        elif params.filtertype == 'recessivepool':
        # homozygous and identical mutations shared by two mixed pools
        # heterozygous mutations shared by recessive mixed pools
            shell("""
            gatk {params.golobe} SelectVariants \
            -R {params.genome} \
            -V {input} \
            -O {output} \
            --invertSelect \
            --selectExpressions '(vc.getGenotype("{params.pool1}").isHomVar() && vc.getGenotype("{params.pool2}").isHomVar()) || vc.getGenotype("{params.repool}").isHet() || (vc.getGenotype("{params.pool1}").isHomRef() && vc.getGenotype("{params.pool2}").isHomRef())'
            """)
        elif params.filtertype == 'poolandparent':
        # homozygous and identical mutations shared by two mixed pools
        # homozygous and identical mutations shared by parents
            shell("""
            gatk {params.golobe} SelectVariants \
            -R {params.genome} \
            -V {input} \
            -O {output} \
            --invertSelect \
            --selectExpressions '(vc.getGenotype("{params.pool1}").isHomVar() && vc.getGenotype("{params.pool2}").isHomVar()) || (vc.getGenotype("{params.parent1}").isHomVar() && vc.getGenotype("{params.parent2}").isHomVar()) || (vc.getGenotype("{params.pool1}").isHomRef() && vc.getGenotype("{params.pool2}").isHomRef()) || (vc.getGenotype("{params.parent1}").isHomRef() && vc.getGenotype("{params.parent2}").isHomRef())'
            """)
        elif params.filtertype == 'recessiveparent':
        # homozygous and identical mutations shared by two mixed pools
        # homozygous and identical mutations shared by parents
        # heterozygous mutations shared by recessive parents or mixed pools
        # mutations from recessive mixed pools that are homozygous and different from the recessive parent
            shell("""
            gatk {params.golobe} SelectVariants \
            -R {params.genome} \
            -V {input} \
            -O {output} \
            --invertSelect \
            --selectExpressions '(vc.getGenotype("{params.pool1}").isHomVar() && vc.getGenotype("{params.pool2}").isHomVar()) || (vc.getGenotype("{params.parent1}").isHomVar() && vc.getGenotype("{params.parent2}").isHomVar()) || vc.getGenotype("{params.reparent}").isHet() || vc.getGenotype("{params.repool}").isHet() || (vc.getGenotype("{params.pool1}").isHomRef() && vc.getGenotype("{params.pool2}").isHomRef()) || (vc.getGenotype("{params.parent1}").isHomRef() && vc.getGenotype("{params.parent2}").isHomRef()) || (vc.getGenotype("{params.repool}").isHomVar() && (vc.getGenotype("{params.reparent}").isHomRef() || vc.getGenotype("{params.reparent}").isHet()))'
            """)
        else:
            exit("Sorry, config file 'poolName's value are not correct")

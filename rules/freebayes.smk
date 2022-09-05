rule freebayes_se:
    input:
        _0 =  "results/snippy_se/{sample}_se/snps.depth.gz",
        i  = "results/snippy_se/{sample}_se/snps.bam",
        _2 = "results/snippy_se/{sample}_se/snps.tab",
        _3 = "results/snippy_se/{sample}_se/snps.consensus.fa",
        ref = REFERENCE
    output:
        o = "results/freebayes_se/{sample}_se/var.vcf"
    conda:
        "../envs/freebayes.yaml"

    params:
        "--min-mapping-quality 20 " 
        "--min-base-quality 20 "  
        "--min-coverage 100 "
        "--min-alternate-count 10 " 
        "--min-alternate-fraction 0.01 "
        "--ploidy 2 "
        "-V " #ver este <_
    shell:
        "freebayes {params} -f {input.ref} {input.i} > {output}"
    
rule freebayes_pe:
    input:
        _0  = "results/snippy_pe/{sample}_pe/snps.depth.gz",
        i   = "results/snippy_pe/{sample}_pe/snps.bam",
        _2  = "results/snippy_pe/{sample}_pe/snps.tab",
        _3  = "results/snippy_pe/{sample}_pe/snps.consensus.fa",
        ref = REFERENCE
    output:
        o = "results/freebayes_pe/{sample}_pe/var.vcf",
        dir = directory("results/freebayes_pe/{sample}_pe/")
    conda:
        "../envs/freebayes.yaml"

    params:
        "--min-mapping-quality 20 " 
        "--min-base-quality 20 "  
        "--min-coverage 100 "
        "--min-alternate-count 10 " 
        "--min-alternate-fraction 0.01 "
    shell:
        "freebayes {params} -f {input.ref} {input.i} > {output.o}"
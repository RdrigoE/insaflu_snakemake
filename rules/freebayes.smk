rule freebayes:
    input:
        _0 = "projects/{project}/samples/snippy/snps.depth.gz",
        i  = "projects/{project}/samples/snippy/snps.bam",
        _2 = "projects/{project}/samples/snippy/snps.tab",
        _3 = "projects/{project}/samples/snippy/snps.consensus.fa",
        ref = REFERENCE
    output:
        o = "projects/{project}/freebayes/{sample}_var.vcf"
    conda:
        "../envs/freebayes.yaml"

    params:
        "--min-mapping-quality 20 " 
        "--min-base-quality 20 "  
        "--min-coverage 100 "
        "--min-alternate-count 10 " 
        "--min-alternate-fraction 0.01 "
        "--ploidy 2 "
        "-V " #ver este <-
    shell:
        "freebayes {params} -f {input.ref} {input.i} > {output}"
    
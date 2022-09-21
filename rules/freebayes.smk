rule freebayes:
    input:
        i  = "align_samples/{sample}/snippy/snps.bam",
        ref = REFERENCE
    output:
        o = "projects/{project}/main_result/freebayes/{sample}_var.vcf",
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
        "freebayes {params} -f {input.ref} {input.i} > {output.o}"
    
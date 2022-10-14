#https://sanjaynagi.github.io/freebayes-parallel/ check this out
#https://github.com/freebayes/freebayes/blob/master/scripts/freebayes-parallel

def get_consensus(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/snps.bam"


rule freebayes:
    input:
        samples = get_consensus,
        ref = REFERENCE
    output:
        o = "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
    # conda:
    #     "../envs/freebayes.yaml"

    params:
        extra = "--min-mapping-quality 20 " 
        "--min-base-quality 20 "  
        "--min-coverage 100 "
        "--min-alternate-count 10 " 
        "--min-alternate-fraction 0.01 "
        "--ploidy 2 "
        "-V " #ver este <-
    # shell:
    #     "freebayes {params} -f {input.ref} {input.i} > {output.o}"
    wrapper:
        "v1.14.0/bio/freebayes"
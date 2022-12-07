# https://sanjaynagi.github.io/freebayes-parallel/ check this out
# https://github.com/freebayes/freebayes/blob/master/scripts/freebayes-parallel


def get_consensus(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/snps.bam"


rule freebayes:
    input:
        samples=get_consensus,
        ref=REFERENCE_FASTA,
    output:
        "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
    threads: 12
    params:
        extra="--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10  --min-alternate-fraction 0.01 --ploidy 2 -V ",  #ver este <-
    retries: 3
    wrapper:
        "v1.14.0/bio/freebayes"


# rule freebayes:
#     input:
#         get_consensus,
#     output:
#         "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
#     conda:
#         "../envs/freebayes.yaml"
#     params:
#         "--min-mapping-quality 20 "
#         "--min-base-quality 20 "
#         "--min-coverage 100 "
#         "--min-alternate-count 10 "
#         "--min-alternate-fraction 0.01 "
#         "--ploidy 2 "
#         "-V ",
#     shell:
#         "freebayes {params} -f {REFERENCE_FASTA} {input} > {output}"

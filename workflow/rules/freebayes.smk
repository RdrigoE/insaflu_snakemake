# rule freebayes:
#     input:
#         samples=get_snps_freebayes,
#         ref=REFERENCE_FASTA,
#     output:
#         "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
#     threads: 12
#     params:
#         extra="--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10  --min-alternate-fraction 0.001 --ploidy 2 -V ",  #ver este <-
#     log:
#         "logs/projects/{project}/{sample}/freebayes_{sample}_var.log",
#     wrapper:
#         "v1.14.0/bio/freebayes"


rule freebayes:
    input:
        samples=get_snps_freebayes,
        ref=REFERENCE_FASTA,
    output:
        "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
    conda:
        "../envs/freebayes.yaml"
    params:
        "--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10  --min-alternate-fraction 0.01 --ploidy 2 -V ",
    log:
        "logs/projects/{project}/{sample}/freebayes_{sample}_var.log",
    benchmark:
        "benchmark/projects/{project}/{sample}/freebayes_{sample}_var.tsv",
    shell:
        "freebayes {params} -f {input.ref} {input.samples} > {output}"

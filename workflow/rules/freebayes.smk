# rule freebayes:
#     input:
#         samples=get_snps_freebayes,
#         ref=REFERENCE_FASTA,
#     output:
#         "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
#     threads: 2
#     params:
#         extra="--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 --min-alternate-fraction 0.01 --ploidy 2 -V ",
#     log:
#         "logs/projects/{project}/{sample}/freebayes_{sample}_var.log",
#     benchmark:
#         "benchmark/projects/{project}/{sample}/freebayes_{sample}_var.tsv"
#     wrapper:
#         "v1.15.0/bio/freebayes"
#
rule freebayes:
    input:
        samples=get_snps_freebayes,
        ref=REFERENCE_FASTA,
    output:
        "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
    conda:
        "../envs/snippy.yaml"
    params:
        "--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10  --min-alternate-fraction 0.01 --ploidy 2 -V ",
    resources:
        mem_mb=memory["freebayes"],
    log:
        "logs/projects/{project}/{sample}/freebayes_{sample}_var.log",
    benchmark:
        "benchmark/projects/{project}/{sample}/freebayes_{sample}_var.tsv"
    shell:
        "freebayes {params} -f {input.ref} {input.samples} > {output}"

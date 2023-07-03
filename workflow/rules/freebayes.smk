rule freebayes:
    input:
        samples=get_snps_freebayes,
        ref=REFERENCE_FASTA,
    output:
        "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
    # conda:
    #     "../envs/freebayes.yaml"
    threads: 8
    params:
        "--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10  --min-alternate-fraction 0.01 --ploidy 2 -V ",
    # resources:
    #     mem_mb=4000,
    log:
        "logs/projects/{project}/{sample}/freebayes_{sample}_var.log",
    benchmark:
        "benchmark/projects/{project}/{sample}/freebayes_{sample}_var.tsv"
    wrapper:
        "v1.14.0/bio/freebayes"

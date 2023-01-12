

rule freebayes:
    input:
        samples=get_consensus_freebayes,
        ref=REFERENCE_FASTA,
    output:
        "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
    threads: 12
    params:
        extra="--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10  --min-alternate-fraction 0.01 --ploidy 2 -V ",  #ver este <-
    # retries: 3
    log:
        "projects/{project}/freebayes_{sample}_var.log",
    wrapper:
        "v1.14.0/bio/freebayes"

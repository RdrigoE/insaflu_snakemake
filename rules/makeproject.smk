rule makeproject:
    input:
        "align_samples/{sample}/snippy/{sample}_consensus.fasta"
    output:
        d3 = "projects/{project}/sample_{sample}/snippy/snps.consensus.fa",
        d4 = "projects/{project}/sample_{sample}/snippy/snps.depth.gz",
        d5 = "projects/{project}/sample_{sample}/snippy/snps.vcf",
        d6 = "projects/{project}/sample_{sample}/snippy/{sample}_consensus.fasta"
    shell:
        "mkdir projects/{wildcards.project}/main_result/depth -p && "
        "mkdir projects/{wildcards.project}/sample_{wildcards.sample}/ -p && "
        " cp -r align_samples/{wildcards.sample}/snippy/ projects/{wildcards.project}/sample_{wildcards.sample}/"
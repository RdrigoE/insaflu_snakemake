rule makeproject:
    input:
        "align_samples/{sample}/snippy/snps.aligned.fa"
    output:
        d2 =directory("projects/{project}/sample_{sample}/snippy/"),
        d3 = "projects/{project}/sample_{sample}/snippy/snps.consensus.fa",
        d4 = "projects/{project}/sample_{sample}/snippy/snps.depth.gz",
        

    shell:
        "mkdir projects/{wildcards.project}/main_result/ -p | "
        "mkdir projects/{wildcards.project}/sample_{wildcards.sample}/snippy/ -p |"
        " cp align_samples/{wildcards.sample}/snippy/snps.* projects/{wildcards.project}/sample_{wildcards.sample}/snippy/"
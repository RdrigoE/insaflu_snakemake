rule makeproject:
    input:
        "align_samples/{sample}/snippy/snps.consensus.fa"
    output:
        #d1 =directory("projects/{project}/main_result/"),
        d2 =directory("projects/{project}/sample_{sample}/snippy/"),
        d3 = "projects/{project}/sample_{sample}/snippy/snps.consensus.fa"

    shell:
        "mkdir projects/{wildcards.project}/main_result/ -p | "
        "mkdir projects/{wildcards.project}/sample_{wildcards.sample}/snippy/ -p |"
        " cp align_samples/{wildcards.sample}/snippy/snps.* projects/{wildcards.project}/sample_{wildcards.sample}/snippy/"
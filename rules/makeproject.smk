def get_consensus(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/{wildcards.sample}_consensus.fasta"


def get_directory(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/*"


rule makeproject:
    input:
        get_consensus,
    output:
        consensus="projects/{project}/sample_{sample}/{sample}_consensus.fasta",
        depth="projects/{project}/sample_{sample}/snps.depth.gz",
        # vcf="projects/{project}/sample_{sample}/snps.vcf",
    params:
        get_directory,
    shell:
        "mkdir projects/{wildcards.project}/main_result/depth -p && "
        "mkdir projects/{wildcards.project}/sample_{wildcards.sample}/ -p && "
        " cp -r {params} projects/{wildcards.project}/sample_{wildcards.sample}/ &&"
        " cp config_user/parameters.yaml projects/{wildcards.project}/"

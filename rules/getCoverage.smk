def get_consensus(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/snps.depth.gz"


rule getCoverage:
    input:
        consensus=get_consensus,
    output:
        coverage=temp("projects/{project}/main_result/{sample}_coverage.csv"),
    conda:
        "../envs/coverage.yaml"
    shell:
        "python software/getCoverage/getCoverage.py -i {input.consensus} -r {REFERENCE} -o {output.coverage}"

def get_consensus(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/snps.depth.gz"


rule getCoverage:
    input:
        consensus=get_consensus,
    output:
        coverage="align_samples/{sample}/{sample}_coverage.csv",
    conda:
        "../envs/coverage.yaml"
    shell:
        "python {coverage_script} -i {input.consensus} -r {REFERENCE_FASTA} -o {output.coverage}"


checkpoint mergeCoverage:
    input:
        expand(
            "align_samples/{sample}/{sample}_coverage.csv",
            sample=config_user["samples"],
        ),
    output:
        coverage_regular=expand(
            "projects/{project}/main_result/coverage.csv",
            project=config_user["project"],
        ),
        coverage_translate=temp(
            expand(
                "projects/{project}/main_result/coverage_translate.csv",
                project=config_user["project"],
            )
        ),
    shell:
        "python {scripts_directory}mergeCoverage.py '{input}' {output.coverage_regular} {output.coverage_translate} {REFERENCE_GB}"

rule getCoverage:
    input:
        depth=get_depth,
    output:
        coverage="align_samples/{sample}/{sample}_coverage.csv",
    conda:
        "../envs/coverage.yaml"
    log:
        "logs/samples/{sample}/coverage.log",
    benchmark:
        "benchmark/samples/{sample}/coverage.tsv"
    shell:
        "python {coverage_script} -i {input.depth} -r {REFERENCE_FASTA} -o {output.coverage}"


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
        coverage_translate=expand(
            "projects/{project}/main_result/coverage_translate.csv",
            project=config_user["project"],
        ),
    conda:
        "../envs/base.yaml"
    log:
        expand(
            "logs/projects/{project}/coverage.log",
            project=config_user["project"],
        ),
    benchmark:
        f"benchmark/projects/{PROJECT_NAME}/coverage.tsv"
    shell:
        "python {scripts_directory}merge_coverage.py '{input}' {output.coverage_regular} {output.coverage_translate} {REFERENCE_GB}"

checkpoint mergeCoverage:
    input:
        expand(
            "projects/{project}/main_result/{sample}_coverage.csv",
            project=config_user["project"],
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
        "python utils/mergeCoverage.py '{input}' {output.coverage_regular} && "
        "python utils/coverage_translate.py '{input}' {output.coverage_translate} {REFERENCE_GB}"

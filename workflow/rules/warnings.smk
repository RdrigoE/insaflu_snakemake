rule warning:
    output:
        "projects/{project}/main_result/warning.txt",
    conda:
        "../envs/base.yaml"
    log:
        "logs/projects/{project}/main_result/warning.log",
    benchmark:
        "benchmark/projects/{project}/main_result/warning.tsv"
    shell:
        "echo 'No sample with all segments above 90% coverage for the value given in the user/parameters.yaml file.' > {output}"

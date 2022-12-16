rule:
    output:
        "projects/{project}/main_result/warning.txt",
    shell:
        "echo 'No sample with all segments above 90% coverage for the value given in the user/parameters.yaml file.' > {output}"

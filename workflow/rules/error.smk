# write rule error that outputs a text file with the message "not enough coverage in any sample to continue the analysis" and exits the pipeline


rule warning_no_coverage:
    conda:
        "../envs/base.yaml"
    output:
        "projects/{project}/main_result/warning.txt",
    log:
        "logs/projects/{project}/main_result/warning.log",
    shell:
        "echo 'Not enough coverage in any sample to continue the analysis' > {output}"

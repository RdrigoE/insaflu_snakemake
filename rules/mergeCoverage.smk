with open('config_user/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

rule mergeCoverage:
    input:
        expand("projects/{project}/main_result/coverage/{sample}_coverage.csv",project=config_user['project'], sample=config_user['samples'])
    output:
        "projects/{project}/main_result/coverage.csv"
    shell:
        "python utils/mergeCoverage.py '{input}' {output}"
with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

locus = get_locus(REFERENCE_GB,config_user['locus'])
if type(get_locus(REFERENCE_GB,config_user['locus'])) == type([1]):
    loop = get_locus(REFERENCE_GB,config_user['locus'])

else:
    loop = 1

checkpoint mergeCoverage:
    input:
        expand("projects/{project}/main_result/coverage/{sample}_coverage.csv",project=config_user['project'], sample=config_user['samples'])
    output:
        coverage_regular = expand("projects/{project}/main_result/coverage.csv",project=config_user['project']),
        coverage_translate = expand("projects/{project}/main_result/coverage_translate.csv",project=config_user['project'])

    shell:
        "python utils/mergeCoverage.py '{input}' {output.coverage_regular} | "
        "python utils/coverage_translate.py '{input}' {output.coverage_translate} {loop}"
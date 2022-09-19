with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

rule snpeff_concat:
    input: 
        expand("projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",sample=config_user['samples'], project=config_user['project']),

    output:
        "projects/{project}/main_result/snpeff.vcf"
    shell:
    	"python utils/snpeff_concat.py '{input}' {output}"
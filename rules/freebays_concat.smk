with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

rule snpeff_concat:
    input: 
        expand("projects/{project}/main_result/freebayes/{sample}_var.vcf",sample=config_user['samples'], project=config_user['project']),

    output:
        "projects/{project}/main_result/validated_minor_iSNVs.csv"
    shell:
    	"python utils/freebays_concat.py '{input}' {output}"
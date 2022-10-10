configfile: "config/parameters.yaml"

with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)
locus = config_user['locus']
samples = config_user['samples']
identification = config_user['identification']
version = config_user['version']
if len(get_locus(REFERENCE_GB)) != 1:
    replace = f"sed -i 's/{locus}/{identification}.{version}/g' "
else:
    replace = 'true '

rule snpeff_sample:
    input:
        i1 = "projects/{project}/sample_{sample}/snippy/snps.vcf",
        i2 = "projects/{project}/main_result/snp_ready.txt"
    output:
        "projects/{project}/sample_{sample}/snippy/{sample}_snpeff.vcf"
    conda:
        "../envs/snpeff.yaml"
    threads: 
        config['snpeff_threads']
    params:
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats -c config/snpeff.config"
    shell:
        "{replace} {input.i1} && "
        "snpEff {params} -v {REFERENCE_NAME} {input.i1}  > {output}" #-t {threads}

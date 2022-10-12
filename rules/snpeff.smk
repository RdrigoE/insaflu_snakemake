configfile: "config/parameters.yaml"

with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)
locus = get_locus(REFERENCE_GB)
samples = config_user['samples']

if len(locus) == 1:
    # print("LOCUS: ", locus[0])
    # print("VERSION: ",get_id_version(REFERENCE_GB))
    replace = f"sed -i 's/{locus[0]}/{get_id_version(REFERENCE_GB)}/g' "
else:
    replace = f'true'


rule prepare_snpeff:
    input:
        ref_gb = REFERENCE_GB,
        ref_fa = REFERENCE,
    output:
        temp("projects/{project}/main_result/snp_ready.txt"),
    conda:
        "../envs/snpeff.yaml"
    shell:
        "python utils/create_snpeff_text.py $CONDA_PREFIX {input.ref_gb} '{locus}' {input.ref_fa} {REFERENCE_NAME} {output} "
        

rule snpeff:
    input:
        i1 = "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
        i2 = "projects/{project}/main_result/snp_ready.txt"
    output:
        o = "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
    conda:
        "../envs/snpeff.yaml"
    threads: 
        config['snpeff_threads']
    params:
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats -c config/snpeff.config"
    shell:
        #not ready for multithreading >< -t {threads}
        "echo '{replace}' && "
        "{replace} {input.i1} &&"
        "snpEff {params} -v {REFERENCE_NAME} {input.i1} > {output.o}"


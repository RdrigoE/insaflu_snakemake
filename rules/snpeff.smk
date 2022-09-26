configfile: "config/parameters.yaml"

with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)
locus = config_user['locus']
samples = config_user['samples']
identification = config_user['identification']
version = config_user['version']
if locus != 'Flu':
    replace = f"sed -i 's/{locus}/{identification}.{version}/g'"
else:
    replace = 'true '


rule prepare_snpeff:
    input:
        ref_gb = REFERENCE_GB,
        ref_fa = REFERENCE,
    output:
        temp("projects/{project}/main_result/snp_ready.txt"),
        # dir("projects/{project}/main_result/")
    conda:
        "../envs/snpeff.yaml"
    shell:
        "python utils/create_snpeff_text.py $CONDA_PREFIX {input.ref_gb} {locus} {REFERENCE_NAME} && "
        "conda activate $CONDA_PREFIX && "
        "mkdir  config/data/{REFERENCE_NAME} -p && "
        "cat {input.ref_gb} > config/data/{REFERENCE_NAME}/genes.gbk && "
        "cat {input.ref_fa} > config/data/{REFERENCE_NAME}/sequences.fa && "
        "snpEff build -genbank {REFERENCE_NAME} -c config/snpeff.config > {output}"

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
        "{replace} {input.i1} &&"
        "snpEff {params} -v {REFERENCE_NAME} {input.i1} > {output.o}"


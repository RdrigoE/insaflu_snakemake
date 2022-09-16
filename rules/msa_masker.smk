with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)
p = config_user["project"]
locus = get_locus(REFERENCE_GB,config_user['locus'])
if type(get_locus(REFERENCE_GB,config_user['locus'])) == type([1]):
    loop = len(get_locus(REFERENCE_GB,config_user['locus']))
else:
    loop = 1


locus = config_user['locus']
rule pre_msa_masker:
    input:
        expand("projects/{p}/main_result/mafft/Alignment_nt_All.fasta",p = config_user["project"]),
    output:
        expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",project = config_user["project"], seg=get_locus(REFERENCE_GB,config_user['locus'])),
        
    shell:
        "python utils/pre_msa_masker.py {input} projects/{p}/main_result {loop} '{locus}'"

rule msa_masker_proteins:
    input:
        all = "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
        i = expand("projects/{project}/main_result/depth/{sample}__{ref}.depth",sample=config_user['samples'], project=config_user['project'], ref=get_locus(run_config["gb_reference"],run_config["locus"])),        

    output:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_masked.fasta"
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c 1"
    shell:
        "python software/msa_masker/msa_masker.py -i {input.all} -df projects/{wildcards.project}/main_result/depth/ -o {output} {params}"

rule msa_masker:
    input:
        all = "projects/{project}/main_result/mafft/Alignment_nt_All_concat.fasta",
        i = expand("projects/{project}/main_result/depth/{sample}.depth",sample=config_user['samples'], project=config_user['project'], ref=get_locus(run_config["gb_reference"],run_config["locus"])),        

    output:
       "projects/{project}/main_result/mafft/Alignment_nt_All_masked.fasta"  
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c 1"
    shell:
        "python software/msa_masker/msa_masker.py -i {input.all} -df projects/{wildcards.project}/main_result/depth/ -o {output} {params}"
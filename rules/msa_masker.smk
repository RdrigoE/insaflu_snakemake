with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)
p = config_user["project"]
locus = get_locus(REFERENCE_GB)
loop = len(locus)


locus = config_user['locus']
rule create_segments:
    input:
        expand("projects/{p}/main_result/Alignment_nt_All_sep.fasta",p = config_user["project"]),
    output:
        expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",project = config_user["project"], seg=get_locus(REFERENCE_GB)),
        
    shell:
        "python utils/pre_msa_masker.py {input} projects/{p}/main_result {loop} '{locus}'"

# rule msa_masker_proteins:
#     input:
#         protein = "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
#         i = expand("projects/{project}/main_result/depth/{sample}__{ref}.depth",sample=config_user['samples'], project=config_user['project'], ref=get_locus(run_config["gb_reference"],run_config["locus"])),        

#     output:
#         temp("projects/{project}/main_result/{seg}/Alignment_nt_{seg}_masked.fasta")
#     conda:
#         "../envs/msa_masker.yaml"
#     params:
#         "--c 1"
#     shell:
#         "python software/msa_masker/msa_masker.py -i {input.protein} -df projects/{wildcards.project}/main_result/depth/ -o {output} {params}"

# rule msa_masker:
#     input:
#         all = "projects/{project}/main_result/Alignment_nt_All_aligned.fasta",
#         i = expand("projects/{project}/main_result/depth/{sample}.depth",sample=config_user['samples'], project=config_user['project'], ref=get_locus(run_config["gb_reference"],run_config["locus"])),        

#     output:
#        "projects/{project}/main_result/Alignment_nt_All.fasta"
#     conda:
#         "../envs/msa_masker.yaml"
#     params:
#         "--c 1"
#     shell:
#         "python software/msa_masker/msa_masker.py -i {input.all} -df projects/{wildcards.project}/main_result/depth/ -o {output} {params}"
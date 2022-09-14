with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

rule pre_msa_masker:
    input:"projects/{project}/main_result/mafft/Alignment_nt_All.fasta"
    output:expand("projects/{project}/main_result/mafft/Alignment_nt_{seg}.fasta",project = config_user["project"], seg=get_locus(REFERENCE_GB,config_user['locus']))
    shell:"python utils.pre_msa_masker.py {input} {seg=get_locus(REFERENCE_GB,config_user['locus'])}"

rule msa_masker:
    input:
        all = "projects/{project}/main_result/mafft/Alignment_nt_All.fasta",
        i = expand("projects/{project}/main_result/depth/{sample}.depth", sample=config_user['samples'], project=config_user['project'])
    output:
        "projects/{project}/main_result/mafft/Alignment_nt_All_masked.fasta"
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c 1"
    shell:
        "python software/msa_masker/msa_masker.py -i {input.all} -df projects/{wildcards.project}/main_result/depth/ -o {output} {params}"
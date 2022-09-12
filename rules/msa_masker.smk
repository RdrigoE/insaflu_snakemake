with open('config_user/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)


rule msa_masker:
    input:
        all = "projects/{project}/main_result/mafft/mafft.fasta",
        i = expand("projects/{project}/main_result/depth/{sample}__SARS_COV_2.depth", sample=config_user['samples'], project=config_user['project'])
    output:
        "projects/{project}/main_result/mafft/mafft_masked.fasta"
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c 1"
    shell:
        "python software/msa_masker/msa_masker.py -i {input.all} -df projects/{wildcards.project}/main_result/depth/ -o {output} {params}"
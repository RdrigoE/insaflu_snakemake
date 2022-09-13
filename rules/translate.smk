with open('config_user/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

p = project=config_user["project"]
r = project=config_user["ref"]

rule translate:
    input:
        ref = REFERENCE_GB,
        i = expand("projects/{project}/main_result/mafft/mafft_masked.fasta",project=config_user["project"])
    output:
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}.fasta",project=config_user["project"],ref=config_user["ref"],protein=config_user["proteins"]),
        dir = directory(f"projects/{p}/main_result/{r}/")
    shell:
        "python utils/translate.py {input.ref} SARS_CoV_2 {input.i} {output.dir}"
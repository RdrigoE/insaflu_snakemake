with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

p = config_user["project"]
locus = get_locus(REFERENCE_GB,config_user["locus"])

rule translate:
    input:
        ref = REFERENCE_GB,
        i = expand("projects/{p}/main_result/mafft/Alignment_nt_All_masked.fasta",p=p)
    output:
        expand("projects/{project}/main_result/{ref}.fasta",ref = locus_protein_alignment, project = config_user["project"]),
        expand("projects/{p}/main_result/{i}",p=p, i = locus),
    shell:
        "python utils/translate.py {input.ref} {wildcards.locus} {input.i} {output.dir}"
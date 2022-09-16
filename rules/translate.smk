with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

rule translate:
    input:
        ref = REFERENCE_GB,
        i = "projects/{project}/main_result/{locus}/Alignment_nt_{locus}_masked.fasta",
        coverage = "projects/{project}/main_result/coverage_translate.csv",
        
    output:
        dir = "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta"
    shell:
        "python utils/translate.py {input.ref} {input.i} {output.dir} '{wildcards.locus}' '{wildcards.gene}' {input.coverage}"
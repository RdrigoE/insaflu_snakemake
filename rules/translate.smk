
rule translate:
    input:
        ref = REFERENCE_GB,
        i = "projects/{project}/main_result/{locus}/Alignment_nt_{locus}.fasta",
        coverage = "projects/{project}/main_result/coverage_translate.csv",
    output:
        temp("projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta")
    shell:
        "python utils/translate.py {input.ref} {input.i} {output} '{wildcards.locus}' '{wildcards.gene}' {input.coverage}"
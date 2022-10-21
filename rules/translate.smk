
rule translate:
    input:
        Alignment_nt = "projects/{project}/main_result/{locus}/Alignment_nt_{locus}.fasta",
        coverage = "projects/{project}/main_result/coverage_translate.csv",
    output:
        temp("projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta")
    shell:
        "python utils/translate.py {REFERENCE_GB} {input.Alignment_nt} {output} '{wildcards.locus}' '{wildcards.gene}' {input.coverage}"
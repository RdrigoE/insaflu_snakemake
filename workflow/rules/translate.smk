rule translate:
    input:
        Alignment_nt="projects/{project}/main_result/{locus}/Alignment_nt_{locus}_mafft.fasta",
        coverage="projects/{project}/main_result/coverage_translate.csv",
    output:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta",
    log:
        "projects/{project}/main_result/{locus}/translate_{locus}_{gene}.log",
    shell:
        "python {scripts_directory}translate.py {REFERENCE_GB} {input.Alignment_nt} {output} '{wildcards.locus}' '{wildcards.gene}' {input.coverage}"

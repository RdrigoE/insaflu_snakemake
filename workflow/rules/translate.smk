rule translate:
    input:
        Alignment_nt="projects/{project}/main_result/{locus}/Alignment_nt_{locus}_mafft.fasta",
        coverage="projects/{project}/main_result/coverage_translate.csv",
    output:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta",
    conda:
        "../envs/base.yaml"
    resources:
        mem_mb=memory["translate"],
    log:
        "logs/projects/{project}/main_result/{locus}/translate_{locus}_{gene}.log",
    benchmark:
        "benchmark/projects/{project}/main_result/{locus}/translate_{locus}_{gene}.tsv"
    shell:
        "python {scripts_directory}translate.py {REFERENCE_GB} {input.Alignment_nt} {output} '{wildcards.locus}' '{wildcards.gene}' {input.coverage}"

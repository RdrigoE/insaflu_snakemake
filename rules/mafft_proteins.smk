rule mafft_proteins:
    input:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta"
    output:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_mafft.fasta"
    conda:
        "../envs/mafft.yaml"
    params:
        "--preservecase --amino"
    shell:
        "mafft {params} {input} > {output}"
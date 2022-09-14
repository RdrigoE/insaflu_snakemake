rule mafft_proteins:
    input:
        "projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}.fasta"
    output:
        "projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_mafft.fasta"
    conda:
        "../envs/mafft.yaml"
    params:
        "--preservecase --amino"
    shell:
        "mafft {params} {input} > {output}"

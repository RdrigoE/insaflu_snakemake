configfile: "config/parameters.yaml"

rule mafft_proteins:
    input:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta"
    output:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_mafft.fasta"
    conda:
        "../envs/mafft.yaml"
    threads:
        config['mafft_threads']
    params:
        "--preservecase --amino"
    shell:
        "mafft --thread {threads} {params} {input} > {output}"
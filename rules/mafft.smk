rule mafft:
    input:
        "projects/{project}/main_result/AllConsensus.fasta"
    output:
        "projects/{project}/main_result/mafft/mafft.fasta" 
    conda:
        "../envs/mafft.yaml"
    params:
        "--preservecase"
    shell:
        "mafft {params} {input} > {output}"


rule mafft:
    input:
        "projects/{project}/main_result/AllConsensus.fasta"
    output:
        "projects/{project}/main_result/mafft/Alignment_nt_All.fasta" 
    conda:
        "../envs/mafft.yaml"
    params:
        "--preservecase"
    shell:
        "mafft {params} {input} > {output}"

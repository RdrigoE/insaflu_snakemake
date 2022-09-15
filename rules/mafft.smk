rule mafft_p_way:
    input:
        "projects/{project}/main_result/AllConsensus_sep.fasta"
    output:
        "projects/{project}/main_result/mafft/Alignment_nt_All.fasta" 
    conda:
        "../envs/mafft.yaml"
    params:
        "--preservecase"
    shell:
        "mafft {params} {input} > {output}"

rule mafft:
    input:
        "projects/{project}/main_result/AllConsensus_concat.fasta"
    output:
        "projects/{project}/main_result/mafft/Alignment_nt_All_concat.fasta" 
    conda:
        "../envs/mafft.yaml"
    params:
        "--preservecase"
    shell:
        "mafft {params} {input} > {output}"

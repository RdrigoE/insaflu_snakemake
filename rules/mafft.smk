configfile: "config/parameters.yaml"


rule mafft_p_way:
    input:
        "projects/{project}/main_result/AllConsensus.fasta"
    output:
        temp("projects/{project}/main_result/Alignment_nt_All_sep.fasta")
    conda:
        "../envs/mafft.yaml"
    threads:
        config['mafft_threads']
    params:
        "--preservecase"
    shell:
        "mafft --thread {threads} {params} {input} > {output}"

rule mafft:
    input:
        "projects/{project}/main_result/All_nt.fasta"
    output:
        temp("projects/{project}/main_result/Alignment_nt_All_aligned.fasta")
    conda:
        "../envs/mafft.yaml"
    threads:
        config['mafft_threads']
    params:
        "--preservecase"
    shell:
        "mafft --thread {threads} {params} {input} > {output}"

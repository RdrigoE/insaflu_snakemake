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
        "projects/{project}/main_result/Alignment_nt_All.fasta"
    conda:
        "../envs/mafft.yaml"
    threads:
        config['mafft_threads']
    params:
        "--preservecase"
    shell:
        "mafft --thread {threads} {params} {input} > {output}"

rule mafft_nt:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta"
    output:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft.fasta"
    conda:
        "../envs/mafft.yaml"
    threads:
        config['mafft_threads']
    params:
        "--preservecase"
    shell:
        "mafft --thread {threads} {params} {input} > {output}"


rule fasttree_nt:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft.fasta"
    output:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_tree.tree"
    conda:
        "../envs/fasttree.yaml"
    params:
        "-gtr -boot 1000"
    shell:
        "fasttree {params} {input} > {output}"
    
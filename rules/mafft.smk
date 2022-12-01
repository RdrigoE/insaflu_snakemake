configfile: "config/parameters.yaml"


rule mafft_p_way:
    input:
        "projects/{project}/main_result/AllConsensus.fasta",
    output:
        temp("projects/{project}/main_result/Alignment_nt_All_sep.fasta"),
    conda:
        "../envs/mafft.yaml"
    threads: config["mafft_threads"]
    params:
        "--preservecase",
    shell:
        "mafft --thread {threads} {params} {input} > {output}"


rule mafft:
    input:
        "projects/{project}/main_result/All_nt.fasta",
    output:
        "projects/{project}/main_result/Alignment_nt_All.fasta",
    conda:
        "../envs/mafft.yaml"
    threads: config["mafft_threads"]
    params:
        "--preservecase",
    shell:
        "mafft --thread {threads} {params} {input} > {output}"


rule mafft_nt:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
    output:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft.fasta",
    conda:
        "../envs/mafft.yaml"
    threads: config["mafft_threads"]
    params:
        "--preservecase",
    shell:
        "mafft --thread {threads} {params} {input} > {output}"


rule fasttree_nt:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft.fasta",
    output:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_tree.tree",
    conda:
        "../envs/fasttree.yaml"
    params:
        "-gtr -boot 1000",
    shell:
        "fasttree {params} {input} > {output}"


rule cp_nt:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
    output:
        "projects/{project}/main_result/Alignment_nt_{seg}.fasta",
    conda:
        "../envs/mafft.yaml"
    shell:
        "cp {input} {output}"


rule cp_nt_tree:
    input:
        tree="projects/{project}/main_result/Tree_ML_All.tree",
        nwk="projects/{project}/main_result/Tree_ML_All.nwk",
    output:
        tree="projects/{project}/main_result/Alignment_nt_{seg}.tree",
        nwk="projects/{project}/main_result/Alignment_nt_{seg}.nwk",
    shell:
        "cp {input.tree} {output.tree} &&"
        "cp {input.nwk} {output.nwk}"

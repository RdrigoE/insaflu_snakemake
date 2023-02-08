configfile: "../config/threads.yaml"


rule align_w_mafft_medaka:
    input:
        align_file="align_samples/{sample}/medaka/medaka_align_{seg}.fasta",
    output:
        aligned_file=temp("align_samples/{sample}/medaka/medaka_aligned_{seg}.fasta"),
    conda:
        "../envs/mafft.yaml"
    threads: config["mafft_threads"]
    params:
        "--preservecase",
    log:
        "logs/align_samples/{sample}/medaka/medaka_align/{seg}.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/medaka_align/{seg}.tsv",
    shell:
        "mafft --thread {threads} {params} {input.align_file} > {output.aligned_file}"


rule align_mafft_snippy:
    input:
        align_file="align_samples/{sample}/snippy/snippy_align_{seg}.fasta",
    output:
        aligned_file=temp("align_samples/{sample}/snippy/snippy_aligned_{seg}.fasta"),
    conda:
        "../envs/mafft.yaml"
    threads: config["mafft_threads"]
    params:
        "--preservecase",
    log:
        "logs/align_samples/{sample}/snippy/snippy_align/{seg}.log",
    benchmark:
        "benchmark/align_samples/{sample}/snippy/snippy_align/{seg}.tsv",
    shell:
        "mafft --thread {threads} {params} {input.align_file} > {output.aligned_file}"


rule mafft_pre_aa:
    input:
        "projects/{project}/main_result/AllConsensus.fasta",
    output:
        temp("projects/{project}/main_result/Alignment_nt_All_sep.fasta"),
    conda:
        "../envs/mafft.yaml"
    threads: config["mafft_threads"]
    params:
        "--preservecase",
    log:
        "logs/projects/{project}/main_result/Alignment_nt_All_sep.log",
    benchmark:
        "benchmark/projects/{project}/main_result/Alignment_nt_All_sep.tsv",
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
    log:
        "logs/projects/{project}/main_result/Alignment_nt_All.log",
    benchmark:
        "benchmark/projects/{project}/main_result/Alignment_nt_All.tsv",
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
    log:
        "logs/projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft.log",
    benchmark:
        "benchmark/projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft.tsv",
    shell:
        "mafft --thread {threads} {params} {input} > {output}"


rule mafft_proteins:
    input:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta",
    output:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_mafft.fasta",
    conda:
        "../envs/mafft.yaml"
    threads: config["mafft_threads"]
    params:
        "--preservecase --amino",
    log:
        "logs/projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_mafft.log",
    benchmark:
        "benchmark/projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_mafft.tsv",
    shell:
        "mafft --thread {threads} {params} {input} > {output}"


rule cp_Alignment_nt:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
    output:
        "projects/{project}/main_result/Alignment_nt_{seg}.fasta",
    conda:
        "../envs/mafft.yaml"
    log:
        "logs/projects/{project}/main_results/copy_Alignment_nt_{seg}.log",
    benchmark:
        "benchmark/projects/{project}/main_results/copy_Alignment_nt_{seg}.tsv",
    shell:
        "cp {input} {output}"

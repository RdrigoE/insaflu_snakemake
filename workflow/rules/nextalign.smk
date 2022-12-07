configfile: "../config/threads.yaml"


rule nextalign_pre_aa:l
    input:
        "projects/{project}/main_result/AllConsensus.fasta",
    output:
        temp("projects/{project}/main_result/Alignment_nt_All_sep.fasta"),
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    params:
        "-j {threads} -n Alignment_nt_All_sep.fasta",
    shell:
        "nextalign run {params} -r {REFERENCE_FASTA} {input}"


# rule nextalign_p_way_1:
#     input:
#         "projects/{project}/main_result/Alignment_nt_All.fasta",
#     output:
#         "projects/{project}/main_result/Alignment_nt_All_sep.fasta",
#     shell:
#         "cp {input} {output}"


rule nextalign:
    input:
        "projects/{project}/main_result/All_nt.fasta",
    output:
        "projects/{project}/main_result/Alignment_nt_All.fasta",
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    params:
        "-j {threads} -n Alignment_nt_All_sep.fasta",
    shell:
        "nextalign run {params} -r {REFERENCE_FASTA} {input}"

rule nextalign_nt:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
    output:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_nextalign.fasta",
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    params:
        "-j {threads} -n Alignment_nt_All_sep.fasta",
    shell:
        "nextalign run {params} -r {REFERENCE_FASTA} {input}"

rule nextalign_proteins:
    input:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta",
    output:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_nextalign.fasta",
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    params:
        "-j {threads} -n Alignment_nt_All_sep.fasta",
    shell:
        "nextalign run {params} -r {REFERENCE_FASTA} {input}"

rule cp_Alignment_nt:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
    output:
        "projects/{project}/main_result/Alignment_nt_{seg}.fasta",
    conda:
        "../envs/nextalign.yaml"
    shell:
        "cp {input} {output}"

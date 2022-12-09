configfile: "../config/threads.yaml"


# rule align_nextAlign_medaka:
#     input:
#         align_file="align_samples/{sample}/medaka/medaka_align_{seg}.fasta",
#     output:
#         aligned_file=temp("align_samples/{sample}/medaka/medaka_aligned_{seg}.fasta"),
#     conda:
#         "../envs/mafft.yaml"
#     threads: 8
#     params:
#         "--preservecase",
#     shell:
#         "nextalign run -r {REFERENCE_FASTA} {params} --output-fasta {output} {input}"


# rule align_nextAlign_snippy:
#     input:
#         align_file="align_samples/{sample}/snippy/snippy_align_{seg}.fasta",
#     output:
#         aligned_file=temp("align_samples/{sample}/snippy/snippy_aligned_{seg}.fasta"),
#     conda:
#         "../envs/mafft.yaml"
#     threads: config["mafft_threads"]
#     params:
#         "--preservecase",
#     shell:
#         "nextalign run -r {REFERENCE_FASTA} {params} --output-fasta {output} {input}"


rule nextalign_pre_aa:
    input:
        "projects/{project}/main_result/AllConsensus.fasta",
    output:
        temp("projects/{project}/main_result/Alignment_nt_All_sep.fasta"),
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    params:
        "-j {threads}",
    shell:
        "nextalign run -r {REFERENCE_FASTA} {params} --output-fasta {output} {input}"


# rule nextalign_p_way_1:
#     input:
#         "projects/{project}/main_result/Alignment_nt_All.fasta",
#     output:
#         "projects/{project}/main_result/Alignment_nt_All_sep.fasta",
#     shell:
#         "cp {input} {output}"


# rule nextalign:
#     input:
#         "projects/{project}/main_result/All_nt.fasta",
#     output:
#         "projects/{project}/main_result/Alignment_nt_All.fasta",
#     conda:
#         "../envs/nextalign.yaml"
#     threads: config["nextalign_threads"]
#     params:
#         "-j {threads}",
#     shell:
#         "nextalign run -r {REFERENCE_FASTA} {params} --output-fasta {output} {input}"


rule nextAlign_pre_aa:
    input:
        "projects/{project}/main_result/AllConsensus.fasta",
    output:
        temp("projects/{project}/main_result/Alignment_nt_All_sep.fasta"),
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    params:
        "-j {threads} ",
    shell:
        "nextalign run -r {REFERENCE_FASTA} {params} --output-fasta {output} {input}"


rule nextalign_nt:
    input:
        "projects/{project}/main_result/AllConsensus.fasta",
    output:
        temp("projects/{project}/main_result/Alignment_nt_All_sep.fasta"),
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    params:
        "-j {threads} ",
    shell:
        "nextalign run -r {REFERENCE_FASTA} {params} --output-fasta {output} {input}"


rule nextalign_proteins:
    input:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta",
    output:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_nextalign.fasta",
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    params:
        "-j {threads}",
    shell:
        "nextalign run -r {REFERENCE_FASTA} {params} --output-fasta {output} {input}"


# rule cp_Alignment_nt:
#     input:
#         "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
#     output:
#         "projects/{project}/main_result/Alignment_nt_{seg}.fasta",
#     conda:
#         "../envs/nextalign.yaml"
#     shell:
#         "cp {input} {output}"

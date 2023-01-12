rule seqret_all_nt:
    input:
        "projects/{project}/main_result/All_nt.fasta",
    output:
        "projects/{project}/main_result/All_nt.nex",
    conda:
        "../envs/seqret.yaml"
    params:
        "-sformat fasta -osformat2 nexusnon",
    log:
        "projects/{project}/main_result/All_nt_to_nex.log",
    shell:
        "seqret {params} -sequence {input} -outseq {output}"


rule seqret_all_nt_only_90:
    input:
        "projects/{project}/main_result/All_nt_only_90plus.fasta",
    output:
        "projects/{project}/main_result/All_nt_only_90plus.nex",
    conda:
        "../envs/seqret.yaml"
    params:
        "-sformat fasta -osformat2 nexusnon",
    log:
        "projects/{project}/main_result/All_nt_only_90plus_to_nex.log",
    shell:
        "seqret {params} -sequence {input} -outseq {output}"


rule seqret_all_consensus:
    input:
        "projects/{project}/main_result/AllConsensus.fasta",
    output:
        "projects/{project}/main_result/AllConsensus.nex",
    conda:
        "../envs/seqret.yaml"
    params:
        "-sformat fasta -osformat2 nexusnon",
    log:
        "projects/{project}/main_result/AllConsensus_to_nex.log",
    shell:
        "seqret {params} -sequence {input} -outseq {output}"


rule seqret_alignment_aa_gene:
    input:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_mafft.fasta",
    output:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_mafft.nex",
    conda:
        "../envs/seqret.yaml"
    params:
        "-sformat fasta -osformat2 nexusnon",
    log:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_mafft_to_nex.log",
    shell:
        "seqret {params} -sequence {input} -outseq {output}"


rule alignment_nt_seg_mafft_seqret:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft.fasta",
    output:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft.nex",
    conda:
        "../envs/seqret.yaml"
    params:
        "-sformat fasta -osformat2 nexusnon",
    log:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft_to_nex.log",
    shell:
        "seqret {params} -sequence {input} -outseq {output}"


rule alignment_nt_seg_seqret:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
    output:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.nex",
    conda:
        "../envs/seqret.yaml"
    params:
        "-sformat fasta -osformat2 nexusnon",
    log:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_to_nex.log",
    shell:
        "seqret {params} -sequence {input} -outseq {output}"


rule alignment_nt_All_seqret:
    input:
        "projects/{project}/main_result/Alignment_nt_All.fasta",
    output:
        "projects/{project}/main_result/Alignment_nt_All.nex",
    conda:
        "../envs/seqret.yaml"
    params:
        "-sformat fasta -osformat2 nexusnon",
    log:
        "projects/{project}/main_result/Alignment_nt_All_to_nex.log",
    shell:
        "seqret {params} -sequence {input} -outseq {output}"


rule alignment_nt_specific_seqret:
    input:
        "projects/{project}/main_result/Alignment_nt_{seg}_mafft.fasta",
    output:
        "projects/{project}/main_result/Alignment_nt_{seg}_mafft.nex",
    conda:
        "../envs/seqret.yaml"
    params:
        "-sformat fasta -osformat2 nexusnon",
    log:
        "projects/{project}/main_result/Alignment_nt_{seg}_mafft_to_nex.log",
    shell:
        "seqret {params} -sequence {input} -outseq {output}"

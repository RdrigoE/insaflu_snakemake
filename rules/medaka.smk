rule medaka_consensus:
    input:
        # i="samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
        i="user_data/{sample}.fastq.gz",
        ref=REFERENCE,
    output:
        out="align_samples/{sample}/medaka/consensus.fasta",
        i="align_samples/{sample}/medaka/calls_to_draft.bam",
        i2="align_samples/{sample}/medaka/snps.bam",
        hdf="align_samples/{sample}/medaka/consensus_probs.hdf",
        ref_out="align_samples/{sample}/medaka/ref.fasta",
    conda:
        "../envs/medaka_1_4_4.yaml"
    shell:
        "cp {input.ref} {output.ref_out} && medaka_consensus -i {input.i} -d {output.ref_out} -o align_samples/{wildcards.sample}/medaka -t 4 -m r941_min_high_g360"
        " && cp {output.i} {output.i2}"


rule medaka_depth:
    input:
        i="align_samples/{sample}/medaka/calls_to_draft.bam",
    output:
        only_depth="align_samples/{sample}/medaka/snps.depth.gz",
        depth="align_samples/{sample}/medaka/snps.depth",
        ind="align_samples/{sample}/medaka/snps.depth.gz.tbi",
    conda:
        "../envs/medaka_1_2_1.yaml"
    params:
        "-aa ",
    shell:
        "samtools depth {params} {input.i} | bgzip -c > {output.only_depth} "
        "&& tabix -p vcf {output.only_depth} && gunzip -c {output.only_depth}  > {output.depth} "


rule medaka_depth_follow:
    input:
        "align_samples/{sample}/medaka/snps.depth",
    output:
        "align_samples/{sample}/medaka/{seg}.depth",
    shell:
        "python3 utils/split_depth_file.py {input} {REFERENCE_GB}"


rule medaka_vfc:
    input:
        hdf="align_samples/{sample}/medaka/consensus_probs.hdf",
        ref=REFERENCE,
        consensus="align_samples/{sample}/medaka/consensus.fasta",
    output:
        vcf="align_samples/{sample}/medaka/round_1.vcf",
    conda:
        "../envs/medaka_1_4_4.yaml"
    shell:
        "medaka variant --verbose {input.ref} {input.hdf} {output.vcf}"


rule medaka_annotate:
    input:
        ref=REFERENCE,
        vcf="align_samples/{sample}/medaka/round_1.vcf",
        bam="align_samples/{sample}/medaka/calls_to_draft.bam",
    output:
        snps="align_samples/{sample}/medaka/snps_ann.vcf",
        vcf_zipped="align_samples/{sample}/medaka/snps_ann.vcf.gz",
    conda:
        "../envs/medaka_1_2_1.yaml"
    shell:
        "medaka tools annotate {input.vcf} {input.ref} {input.bam} {output.snps} && bgzip {output.snps} -c > {output.vcf_zipped} && tabix {output.vcf_zipped}"


def get_add_freq_medaka(software_parameters):
    return f'{software_parameters["mincov_medaka"]} {software_parameters["min_prop_for_variant_evidence"]}'


rule mask_between_top_and_50:
    input:
        vcf_file="align_samples/{sample}/medaka/snps_ann.vcf",
        file_coverage="align_samples/{sample}/medaka/snps.depth.gz",
        normal_reference_fasta="align_samples/{sample}/medaka/ref.fasta",
    output:
        vcf_file_out="align_samples/{sample}/medaka/snps.vcf",
        vcf_file_out_removed_by_filter="align_samples/{sample}/medaka/snps_removed_by_filter.vcf",
        temp_file=temp("align_samples/{sample}/medaka/temp_file.txt"),
        normal_reference_fasta_removed="align_samples/{sample}/medaka/ref_masked.fasta",
        vcf_file_out_compr="align_samples/{sample}/medaka/snps.vcf.gz",
    conda:
        "../envs/filter_medaka.yaml"
    params:
        get_add_freq_medaka,
    shell:
        "touch {output.temp_file} && "
        "python utils/add_freq_medaka_consensus.py {input.normal_reference_fasta} {input.vcf_file} {input.file_coverage} "
        "{output.vcf_file_out} {output.vcf_file_out_removed_by_filter} {output.temp_file} {output.normal_reference_fasta_removed} {output.temp_file} {params} &&"
        "bgzip -c -f {output.vcf_file_out} > {output.vcf_file_out_compr} && tabix {output.vcf_file_out_compr}"


rule bcf_consensus:
    input:
        vcf_ziped="align_samples/{sample}/medaka/snps.vcf.gz",
        ref="align_samples/{sample}/medaka/ref_masked.fasta",
    output:
        temp("align_samples/{sample}/medaka/first_consensus.fasta"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools consensus -s SAMPLE -f {REFERENCE} {input.vcf_ziped} -o {output}"


rule create_align_file:
    input:
        first_consensus="align_samples/{sample}/medaka/first_consensus.fasta",
    output:
        align_file=temp("align_samples/{sample}/medaka/medaka_align_{seg}.fasta"),
    shell:
        "python utils/mask_consensus_by_deep.py {REFERENCE} {input.first_consensus} {output.align_file} {wildcards.seg}"


rule align_mafft:
    input:
        align_file="align_samples/{sample}/medaka/medaka_align_{seg}.fasta",
    output:
        aligned_file=temp("align_samples/{sample}/medaka/medaka_aligned_{seg}.fasta"),
    conda:
        "../envs/mafft.yaml"
    threads: 8
    params:
        "--preservecase",
    shell:
        "mafft --thread {threads} {params} {input.align_file} > {output.aligned_file}"


rule msa_masker_medaka:
    input:
        align_file="align_samples/{sample}/medaka/medaka_aligned_{seg}.fasta",
        depth="align_samples/{sample}/medaka/{seg}.depth",
    output:
        temp("align_samples/{sample}/medaka/consensus_aligned_{seg}.fasta"),
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c " + str(software_parameters["mincov_medaka"]),
    shell:
        "python software/msa_masker/msa_masker.py -i {input.align_file} -df {input.depth} -o {output} {params}"


rule get_masked_consensus_medaka:
    input:
        lambda wildcards: expand(
            "align_samples/{sample}/medaka/consensus_aligned_{seg}.fasta",
            sample=wildcards.sample,
            seg=get_locus(REFERENCE_GB),
        ),
    output:
        final_consensus="align_samples/{sample}/medaka/pre_{sample}_consensus.fasta",
    shell:
        "python utils/get_consensus_medaka.py '{input}' {output}"


rule mask_regions_consensus_medaka:
    input:
        consensus="align_samples/{sample}/medaka/pre_{sample}_consensus.fasta",
    output:
        final_consensus="align_samples/{sample}/medaka/{sample}_consensus.fasta",
    params:
        mask_regions_parameters(software_parameters),
    shell:
        "python utils/mask_regions.py {input} {output} {params}"


# rule copy_to_compare:
#     input:
#         "align_samples/Demo_Sample_103/medaka/Demo_Sample_103_consensus.fasta",
#     output:
#         "test_results/AllConsensus_no_ref.fasta",
#     conda:
#         "test_insaflu"
#     shell:
#         "cp {input} {output} && "
#         "sed -i 's/Demo_Sample_/Sample/g' {output} && "
#         "python test_results/analyse_results.py test_results/AllConsensus.fasta {output} test_results/SARS_CoV_2_COVID_19_Wuhan_Hu_1_MN908947.gbk "
#         " test_results/show_me.csv && cat test_results/show_me.csv"

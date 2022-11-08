rule align_pe:
    input:
        reads_1="samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        reads_2="samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",
    output:
        bam="align_samples/{sample}/iVar/snps.bam",
    conda:
        "../envs/ivar.yaml"
    shell:
        "mkdir align_samples/{wildcards.sample}/reference -p && "
        "cp {REFERENCE} align_samples/{wildcards.sample}/reference/ && "
        "bwa index align_samples/{wildcards.sample}/{REFERENCE} && "
        "bwa mem align_samples/{wildcards.sample}/{REFERENCE} {input.reads_1} {input.reads_2} | "
        "samtools view -u -T align_samples/{wildcards.sample}/{REFERENCE} -q 20 | "
        "samtools sort --reference align_samples/{wildcards.sample}/{REFERENCE} > {output}"


rule align_se:
    input:
        reads_1="samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz",
    output:
        "align_samples/{sample}/iVar/snps.bam",
    conda:
        "../envs/ivar.yaml"
    shell:
        "bwa index align_samples/{wildcards.sample}/{REFERENCE} && "

        "bwa mem align_samples/{wildcards.sample}/{REFERENCE} {input.reads_1} | "
        "samtools view -u -T align_samples/{wildcards.sample}/{REFERENCE} -q 20 | "
        "samtools sort --reference align_samples/{wildcards.sample}/{REFERENCE} > {output} && touch snps.vcf"


ruleorder: align_pe > align_se

# rule conver_gbk_to_gff:
#     input:"reference/{REFERENCE_NAME}.gbk"
#     output:"reference/{REFERENCE_NAME}.gff"
#     conda:"../envs/perl.yaml"
#     shell:"perl software/bp_genbank2gff3.pl {input}"

rule call_variant:
    input:
        bam="align_samples/{sample}/iVar/snps.bam",
        gff = "reference/SARS_CoV_2_COVID_19_Wuhan_Hu_1_MN908947.gff3"
    output:
        vcf_file="align_samples/{sample}/iVar/snps.tsv",
    conda:
        "../envs/ivar.yaml"
    params:
        "snps",
    shell:
        "samtools mpileup -A -d 600000 -B -Q 0 {input.bam} |"
        "ivar variants -p {params} -q 20 -t 0.51 -m 10 -r align_samples/{wildcards.sample}/{REFERENCE} -g {input.gff} "


rule filter_variants:
    input:
        bam="align_samples/{sample}/iVar/snps.bam",
        vcf_file="align_samples/{sample}/iVar/snps.tsv",
    output:
        vcf_file="align_samples/{sample}/iVar/snps_filtered.tsv",
    conda:
        "../envs/ivar.yaml"
    params:
        "snps_filtered",
    shell:
        "ivar filtervariants -t 0.51 -p {params} {output.vcf_file}"


rule generate_consensus:
    input:
        bam="align_samples/{sample}/iVar/snps.bam",
    output:
        consensus="align_samples/{sample}/iVar/snps.consensus.fa",
    conda:
        "../envs/ivar.yaml"
    params:
        "snps.consensus",
    shell:
        "samtools mpileup -aa -A -Q 0 {input.bam} | ivar consensus -p align_samples/{wildcards.sample}/iVar/{params} -q 10 -t 0.51"
        # "&& python3 utils/change_id_name.py {output.consensus} SARS_CoV_2"


rule get_depth:
    input:
        i="align_samples/{sample}/iVar/snps.bam",
        consensus="align_samples/{sample}/iVar/snps.consensus.fa",
    output:
        consensus="align_samples/{sample}/iVar/new_snps.consensus.fa",
        only_depth="align_samples/{sample}/iVar/snps.depth.gz",
        depth="align_samples/{sample}/iVar/snps.depth",
        ind="align_samples/{sample}/iVar/snps.depth.gz.tbi",
    conda:
        "../envs/medaka_1_4_4.yaml"
    params:
        "-aa -q 20",
    shell:
        "samtools depth {params} {input.i} | bgzip -c > {output.only_depth} "
        "&& tabix -p vcf {output.only_depth} && gunzip -c {output.only_depth}  > {output.depth} &&"
        "python utils/split_consensus.py {input.consensus} {output.depth} {REFERENCE_GB} {output.consensus}"


rule iVar_depth_step_2:
    input:
        zipped="align_samples/{sample}/iVar/snps.depth",
    output:
        unzipped="align_samples/{sample}/iVar/{seg}.depth",
    shell:
        "python utils/split_depth_file.py {input.zipped} {REFERENCE_GB}"


rule create_align_file_iVar:
    input:
        first_consensus="align_samples/{sample}/iVar/new_snps.consensus.fa",
    output:
        align_file=temp("align_samples/{sample}/iVar/iVar_align_{seg}.fasta"),
    shell:
        "python utils/mask_consensus_by_deep.py align_samples/{wildcards.sample}/{REFERENCE} {input.first_consensus} {output.align_file} {wildcards.seg}"


rule align_mafft_iVar:
    input:
        align_file="align_samples/{sample}/iVar/iVar_align_{seg}.fasta",
    output:
        aligned_file=temp("align_samples/{sample}/iVar/iVar_aligned_{seg}.fasta"),
    conda:
        "../envs/mafft.yaml"
    threads: 12
    params:
        "--preservecase",
    shell:
        "mafft --thread {threads} {params} {input.align_file} > {output.aligned_file}"


rule msa_masker_iVar:
    input:
        align_file="align_samples/{sample}/iVar/iVar_aligned_{seg}.fasta",
        depth="align_samples/{sample}/iVar/{seg}.depth",
    output:
        temp("align_samples/{sample}/iVar/consensus_aligned_{seg}.fasta"),
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c " + str(software_parameters["msa_masker"]),
    shell:
        "python software/msa_masker/msa_masker.py -i {input.align_file} -df {input.depth} -o {output} {params}"


rule get_masked_consensus_iVar:
    input:
        lambda wildcards: expand(
            "align_samples/{sample}/iVar/consensus_aligned_{seg}.fasta",
            sample=wildcards.sample,
            seg=get_locus(REFERENCE_GB),
        ),
    output:
        final_consensus="align_samples/{sample}/iVar/pre_{sample}_consensus.fasta",
    shell:
        "python utils/get_consensus_medaka.py '{input}' {output}"


rule mask_regions_consensus_iVar:
    input:
        consensus="align_samples/{sample}/iVar/pre_{sample}_consensus.fasta",
    output:
        final_consensus="align_samples/{sample}/iVar/{sample}_consensus.fasta",
    params:
        mask_regions_parameters(software_parameters),
    shell:
        "python utils/mask_regions.py {input} {output} {params}"


# rule check_equals:
#     input:
#         "consensus.fasta",
#     output:
#         "zzz.csv",
#     conda:
#         "test_insaflu"
#     shell:
#         "python3 utils/change_id_name.py {input} SARSCOV2_1 && "
#         "python3 analyse_results.py {input} AllConsensus.fasta reference.gbk {output}"

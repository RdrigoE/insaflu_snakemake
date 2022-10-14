from utils.get_software_parameters import get_snippy_parameters,mask_regions_parameters
configfile: "config/parameters.yaml"

rule snippy_pe:
    input:
        r1="samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        r2="samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",
        ref=REFERENCE
    output:
        _0  = "align_samples/{sample}/snippy/snps.depth.gz",
        _1 = "align_samples/{sample}/snippy/snps.bam",
        _2 = "align_samples/{sample}/snippy/snps.tab",
        _3 = "align_samples/{sample}/snippy/snps.consensus.fa",
    conda:
        "../envs/snippy.yaml"
    threads: 
        config['snippy_threads']
    params:
        get_snippy_parameters(software_parameters)
    shell:
        "rm -r align_samples/{wildcards.sample}/snippy/ && "
        "snippy --cpus {threads} --pe1 {input.r1} --pe2 {input.r2} --ref {input.ref} --outdir align_samples/{wildcards.sample}/snippy/ {params} "


rule snippy_se:
    input:
        r1="samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz",
        ref=REFERENCE
    output:
        _0  = "align_samples/{sample}/snippy/snps.depth.gz",
        _1 = "align_samples/{sample}/snippy/snps.bam",
        _2 = "align_samples/{sample}/snippy/snps.tab",
        _3 = "align_samples/{sample}/snippy/snps.consensus.fa",

    conda:
        "../envs/snippy.yaml"
    params:
        get_snippy_parameters(software_parameters)
    shell:
        "rm -r align_samples/{wildcards.sample}/snippy/ && " 
        "snippy --se {input.r1} --ref {input.ref} --outdir align_samples/{wildcards.sample}/snippy/ {params}"


ruleorder: snippy_pe > snippy_se


rule snippy_depth_step_1:
    input:
        zipped = "align_samples/{sample}/snippy/snps.depth.gz",
    output:
        out = "align_samples/{sample}/snippy/depth/snps.depth",
    shell:
        "gunzip -k -c {input.zipped} > {output}"


rule snippy_depth_step_2:
    input:
        zipped = "align_samples/{sample}/snippy/depth/snps.depth",
    output:
        out = "align_samples/{sample}/snippy/depth/{seg}.depth",
    shell:
        "python utils/split_depth_file.py align_samples/{wildcards.sample}/snippy/depth/snps.depth {REFERENCE_GB}"


rule create_align_file_snippy:
    input:
        first_consensus = "align_samples/{sample}/snippy/snps.consensus.fa",
        ref = REFERENCE
    output:
        align_file = temp("align_samples/{sample}/snippy/snippy_align_{seg}.fasta"),
    shell:
        "python utils/mask_consensus_by_deep.py {input.ref} {input.first_consensus} {output.align_file} {wildcards.seg}"


rule align_mafft_snippy:
    input:
        align_file = "align_samples/{sample}/snippy/snippy_align_{seg}.fasta"
    output:
        aligned_file = temp("align_samples/{sample}/snippy/snippy_aligned_{seg}.fasta")
    conda:
        "../envs/mafft.yaml"
    threads:
        config['mafft_threads']
    params:
        "--preservecase"
    shell:
        "mafft --thread {threads} {params} {input.align_file} > {output.aligned_file}"


rule msa_masker_snippy:
    input:
        align_file = "align_samples/{sample}/snippy/snippy_aligned_{seg}.fasta",
        depth = "align_samples/{sample}/snippy/depth/{seg}.depth"
    output:
       temp("align_samples/{sample}/snippy/consensus_aligned_{seg}.fasta")
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c 1"
    shell:
        "python software/msa_masker/msa_masker.py -i {input.align_file} -df {input.depth} -o {output} {params}"


rule get_masked_consensus_snippy:
    input:
        lambda wildcards:
            expand("align_samples/{sample}/snippy/consensus_aligned_{seg}.fasta", sample = wildcards.sample, seg = get_locus(REFERENCE_GB))
    output:
       final_consensus = "align_samples/{sample}/snippy/pre_{sample}_consensus.fasta"
    shell:
        "python utils/get_consensus_medaka.py '{input}' {output}"


rule mask_regions_consensus_snippy:
    input:
        consensus = "align_samples/{sample}/snippy/pre_{sample}_consensus.fasta"
    output:
        final_consensus = "align_samples/{sample}/snippy/{sample}_consensus.fasta"
    params:
        mask_regions_parameters(software_parameters)
    shell: 
        "python utils/mask_regions.py {input} {output} {params}"
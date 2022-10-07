rule medaka_consensus:
    input:
        i = "samples/{sample}/nano_trimmed_reads/{sample}.trimmed.fastq.gz",
        ref = REFERENCE
    output:
        out = "align_samples/{sample}/medaka/consensus.fasta",
        i = "align_samples/{sample}/medaka/calls_to_draft.bam",
        hdf = "align_samples/{sample}/medaka/consensus_probs.hdf"
    conda:
        "../envs/medaka_1_4_4.yaml"
    shell:
        "medaka_consensus -i {input.i} -d {input.ref} -o align_samples/{wildcards.sample}/medaka -t 4 -m r941_min_high_g360"


rule medaka_depth:
    input:
        i = "align_samples/{sample}/medaka/calls_to_draft.bam",
    output:
        depth = "align_samples/{sample}/medaka/depth/{seg}.depth.gz",
        depth_unzipped = "align_samples/{sample}/medaka/depth/{seg}.depth",
        depth_indexed = "align_samples/{sample}/medaka/depth/{seg}.depth.gz.tbi"
    conda:
        "../envs/medaka_1_4_4.yaml"
    params:
        "-aa -q 10"
    shell:
        "samtools depth {params} {input.i} | bgzip -c > {output.depth} "
        "&& tabix -p vcf {output.depth} && gunzip -k {output.depth} && python utils/split_depth_file.py {output.depth_unzipped} {REFERENCE_GB}"

rule medaka_vfc:
    input:
        hdf = "align_samples/{sample}/medaka/consensus_probs.hdf",
        ref = REFERENCE,
        consensus = "align_samples/{sample}/medaka/consensus.fasta",
        depth = "align_samples/{sample}/medaka/depth/snps.depth.gz.tbi",
    output:
        vcf = "align_samples/{sample}/medaka/round_1.vcf",
    conda:
        "../envs/medaka_1_4_4.yaml"
    shell:
        "medaka variant --verbose {input.ref} {input.hdf} {output.vcf}"

rule medaka_annotate:
    input:
        ref = REFERENCE,
        vcf = "align_samples/{sample}/medaka/round_1.vcf",
        bam = "align_samples/{sample}/medaka/calls_to_draft.bam"
    output:
        snps = "align_samples/{sample}/medaka/snps.vcf",
        vcf_zipped = "align_samples/{sample}/medaka/snps.vcf.gz"
    conda:
        "../envs/medaka_1_4_4.yaml"
    shell:
        "medaka tools annotate {input.vcf} {input.ref} {input.bam} {output.snps} && bgzip {output.snps} -c > {output.vcf_zipped} && tabix {output.vcf_zipped}"


rule bcf_consensus:
    input:
        vcf_ziped = "align_samples/{sample}/medaka/snps.vcf.gz",
        ref = REFERENCE
    output:
        out = temp("align_samples/{sample}/medaka/first_consensus.fasta"),
    conda:
        "../envs/medaka_1_4_4.yaml"
    shell:
        "bcftools consensus -s SAMPLE -f {input.ref} {input.vcf_ziped} -o {output.out}"

# rule mask first consensus
rule create_align_file:
    input:
        first_consensus = "align_samples/{sample}/medaka/first_consensus.fasta",
        ref = REFERENCE
    output:
        align_file = temp("align_samples/{sample}/medaka/medaka_align_{seg}.fasta"),
    shell:
        "python utils/mask_consensus_by_deep.py {input.ref} {input.first_consensus} {output.align_file} {wildcards.seg}"

rule align_mafft:
    input:
        align_file = "align_samples/{sample}/medaka/medaka_align_{seg}.fasta"
    output:
        aligned_file = temp("align_samples/{sample}/medaka/medaka_aligned_{seg}.fasta")
    conda:
        "../envs/mafft.yaml"
    threads:
        config['mafft_threads']
    params:
        "--preservecase"
    shell:
        "mafft --thread {threads} {params} {input.align_file} > {output.aligned_file}"

rule msa_masker_medaka:
    input:
        align_file = "align_samples/{sample}/medaka/medaka_aligned_{seg}.fasta",
        depth = "align_samples/{sample}/medaka/depth/{seg}.depth"
    output:
       temp("align_samples/{sample}/medaka/consensus_aligned_{seg}.fasta")
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c 1"
    shell:
        "python software/msa_masker/msa_masker.py -i {input.align_file} -df {input.depth} -o {output} {params}"

rule get_masked_consensus:
    input:
        lambda wildcards:
            expand("align_samples/{sample}/medaka/consensus_aligned_{seg}.fasta", sample = wildcards.sample, seg = get_locus(REFERENCE_GB))
    output:
       final_consensus = "align_samples/{sample}/medaka/{sample}_consensus.fasta"
    shell:
        "python utils/get_consensus.py '{input}' {output}"


# rule add freq to vcf 

# rule add snpeff
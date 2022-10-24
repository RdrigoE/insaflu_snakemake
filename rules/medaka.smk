rule medaka_consensus:
    input:
        i = "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
        ref = REFERENCE
    output:
        out = "align_samples/{sample}/medaka/consensus.fasta",
        i = "align_samples/{sample}/medaka/calls_to_draft.bam",
        i2 = "align_samples/{sample}/medaka/snps.bam",
        hdf = "align_samples/{sample}/medaka/consensus_probs.hdf"
    conda:
        "../envs/medaka_1_4_4.yaml"
    shell:
        "medaka_consensus -i {input.i} -d {input.ref} -o align_samples/{wildcards.sample}/medaka -t 4 -m r941_min_high_g360"
        " && cp {output.i} {output.i2}"


rule medaka_depth:
    input:
        i = "align_samples/{sample}/medaka/calls_to_draft.bam",
    output:
        only_depth =  "align_samples/{sample}/medaka/snps.depth.gz",
        depth = "align_samples/{sample}/medaka/snps.depth",
        ind = "align_samples/{sample}/medaka/snps.depth.gz.tbi",
    conda:
        "../envs/medaka_1_4_4.yaml"
    params:
        "-aa -q 20"
    shell:
        "samtools depth {params} {input.i} | bgzip -c > {output.only_depth} "
        # "samtools depth {input.i} | bgzip -c > {output.only_depth} "

        "&& tabix -p vcf {output.only_depth} && gunzip {output.only_depth} > {output.depth} "

rule medaka_depth_follow:
    input: "align_samples/{sample}/medaka/snps.depth"
    output: "align_samples/{sample}/medaka/{seg}.depth"
    shell: "python3 utils/split_depth_file.py {input} {REFERENCE_GB}" 

rule medaka_vfc:
    input:
        hdf = "align_samples/{sample}/medaka/consensus_probs.hdf",
        ref = REFERENCE,
        consensus = "align_samples/{sample}/medaka/consensus.fasta",
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
        vcf_ziped = "align_samples/{sample}/medaka/snps.vcf.gz"
    output:
        temp("align_samples/{sample}/medaka/first_consensus.fasta"),
    conda:
        "../envs/medaka_1_4_4.yaml"
    shell:
        "bcftools consensus -s SAMPLE -f {REFERENCE} {input.vcf_ziped} -o {output}"


rule create_align_file:
    input:
        first_consensus = "align_samples/{sample}/medaka/first_consensus.fasta"
    output:
        align_file = temp("align_samples/{sample}/medaka/medaka_align_{seg}.fasta"),
    shell:
        "python utils/mask_consensus_by_deep.py {REFERENCE} {input.first_consensus} {output.align_file} {wildcards.seg}"

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
        depth = "align_samples/{sample}/medaka/{seg}.depth",
    output:
       temp("align_samples/{sample}/medaka/consensus_aligned_{seg}.fasta")
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c " + str(software_parameters['mincov_medaka'])
    shell:
        "python software/msa_masker/msa_masker.py -i {input.align_file} -df {input.depth} -o {output} {params}"



rule get_masked_consensus_medaka:
    input:
        lambda wildcards:
            expand("align_samples/{sample}/medaka/consensus_aligned_{seg}.fasta", sample = wildcards.sample, seg = get_locus(REFERENCE_GB))
    output:
       final_consensus = "align_samples/{sample}/medaka/pre_{sample}_consensus.fasta"
    shell:
        "python utils/get_consensus_medaka.py '{input}' {output}"


    
rule mask_regions_consensus_medaka:
    input:
        consensus = "align_samples/{sample}/medaka/pre_{sample}_consensus.fasta"
    output:
        final_consensus = "align_samples/{sample}/medaka/{sample}_consensus.fasta"
    params:
        mask_regions_parameters(software_parameters)
    shell:
        "python utils/mask_regions.py {input} {output} {params}"

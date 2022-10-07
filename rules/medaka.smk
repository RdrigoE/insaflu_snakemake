rule medaka_consensus:
    input:
        i = "samples/{sample}/nano_trimmed_reads/{sample}.trimmed.fastq.gz",
        ref = REFERENCE
    output:
        dir = directory("align_samples/{sample}/medaka"),
        out = "align_samples/{sample}/medaka/consensus.fasta",
        i = "align_samples/{sample}/medaka/calls_to_draft.bam",
        hdf = "align_samples/{sample}/medaka/consensus_probs.hdf"
    conda:
        "../envs/medaka_1_4_4.yaml"
    shell:
        "medaka_consensus -i {input.i} -d {input.ref} -o {output.dir} -t 4 -m r941_min_high_g360"


rule medaka_depth:
    input:
        i = "align_samples/{sample}/medaka/calls_to_draft.bam",
    output:
        depth = "align_samples/{sample}/medaka/snps.depth.gz",
        depth_indexed = "align_samples/{sample}/medaka/snps.depth.gz.tbi"
    conda:
        "../envs/medaka_1_4_4.yaml"
    params:
        "-aa -q 10"
    shell:
        "samtools depth {params} {input.i} | bgzip -c > {output.depth} "
        "&& tabix -p vcf {output.depth}"

rule medaka_vfc:
    input:
        hdf = "align_samples/{sample}/medaka/consensus_probs.hdf",
        ref = REFERENCE,
        consensus = "align_samples/{sample}/medaka/consensus.fasta",
        depth = "align_samples/{sample}/medaka/snps.depth.gz.tbi",
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
        out = "align_samples/{sample}/medaka/final_consensus.fasta",
    conda:
        "../envs/medaka_1_4_4.yaml"
    shell:
        "bcftools consensus -s SAMPLE -f {input.ref} {input.vcf_ziped} -o {output.out}"

# rule mask final consensus
# rule bcf_consensus:
#     input:
#         vcf_ziped = "align_samples/{sample}/medaka/snps.vcf.gz",
#         ref = REFERENCE
#     output:
#         out = "align_samples/{sample}/medaka/final_consensus.fasta",
#     conda:
#         "../envs/medaka_1_4_4.yaml"
#     shell:
#         "bcftools consensus -s SAMPLE -f {input.ref} {input.vcf_ziped} -o {output.out}"

# rule mask_consensus:
#     input:
#         consensus = "align_samples/{sample}/medaka/final_consensus.fasta",
#         depth = "align_samples/{sample}/medaka/snps.depth"
#     output:
#         consensus = 
#     conda:

#     shell:

# rule add freq to vcf 

# rule add snpeff
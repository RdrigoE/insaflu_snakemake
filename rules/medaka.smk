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
        "../envs/medaka_1_7.yaml"
    shell:
        "medaka_consensus -i {input.i} -d {input.ref} -o {output.dir} -t 8 -m r941_min_high_g360"
        #email => Medaka_consensus -i filtered_fastq -d reference_fasta -o output_folder -t threads –m r941_min_high_g360
rule medaka_depth:
    input:
        i = "align_samples/{sample}/medaka/calls_to_draft.bam",
    output:
        depth = "align_samples/{sample}/medaka/snps.depth.gz"
    conda:
        "../envs/medaka_1_2.yaml"
    params:
        "-aa -q 10"
    shell:
        "samtools depth {params} {input.i} | bgzip -c > {output.depth}; tabix {output.depth}"

rule tabix:
    input:
        depth = "align_samples/{sample}/medaka/snps.depth.gz",
    output:
        depth = "align_samples/{sample}/medaka/snps.depth.gz.tbi"
    conda:
        "../envs/tabix.yaml"
    shell:
        "tabix {input.depth}"

rule medaka_vfc:
    input:
        hdf = "align_samples/{sample}/medaka/consensus_probs.hdf",
        ref = REFERENCE,
        bam = "align_samples/{sample}/medaka/calls_to_draft.bam",
        depth = "align_samples/{sample}/medaka/snps.depth.gz.tbi"
    output:
        vcf = "align_samples/{sample}/medaka/snps.vfc",
        dir = 'align_samples/{sample}/medaka/medaka_variant'
    conda:
        "../envs/medaka_1_2.yaml"

    shell:
        "medaka_variant -i {input.bam} -f {input.ref} -m {input.hdf} -o {output.dir} -p {output.vcf}"



# rule medaka_annotate:
#     input:
#         ref = REFERENCE,
#         vcf = "align_samples/{sample}/medaka/snps.vfc",
#         bam = "align_samples/{sample}/medaka/calls_to_draft.bam"
#     output:
#         vfc = "align_samples/{sample}/medaka/snps.vfc"
#     conda:
#         "../envs/medaka_1_2.yaml"
#     shell:
#         "medaka_tools annotate {input.vcf} -f {input.ref} {input.bam} {output}"

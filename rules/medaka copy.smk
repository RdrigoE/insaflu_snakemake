rule medaka_consensus:
    input:
        i = "samples/{sample}/nano_trimmed_reads/nano_{sample}.trimmed.fastq.gz",
        ref = REFERENCE
    output:
        dir = directory("align_samples/{sample}/medaka"),
        out = "align_samples/{sample}/medaka/consensus.fasta",
        i = "align_samples/{sample}/medaka/calls_to_draft.bam",
        hdf = "align_samples/{sample}/medaka/consensus_probs.hdf"
    conda:
        "../envs/medaka_1_7.yaml"
    shell:
        "source ./testONT/software/medaka/bin/activate && "
        "./testONT/software/medaka/bin/medaka_consensus -i {input.i} -d {input.ref} -o {output.dir} -t 4 -m r941_min_high_g360"
        #email => Medaka_consensus -i filtered_fastq -d reference_fasta -o output_folder -t threads –m r941_min_high_g360
rule medaka_depth:
    input:
        i = "align_samples/{sample}/medaka/calls_to_draft.bam",
    output:
        depth = "align_samples/{sample}/medaka/snps.depth.gz",
        depth_indexed = "align_samples/{sample}/medaka/snps.depth.gz.tbi"
    conda:
        "../envs/medaka_1_2.yaml"
    params:
        "-aa -q 10"
    shell:
        "samtools depth {params} {input.i} | bgzip -c > {output.depth} "
        "&& tabix -p vcf {output.depth}"

# rule tabix:
#     input:
#         depth = "align_samples/{sample}/medaka/snps.depth.gz",
#     output:
#         depth = "align_samples/{sample}/medaka/snps.depth.gz.tbi"
#     conda:
#         "../envs/tabix.yaml"
#     shell:
#         "tabix {input.depth}"

rule medaka_vfc:
    input:
        hdf = "align_samples/{sample}/medaka/consensus_probs.hdf",
        ref = REFERENCE,
        bam = "align_samples/{sample}/medaka/calls_to_draft.bam",
        depth = "align_samples/{sample}/medaka/snps.depth.gz.tbi"
    output:
        vcf = "align_samples/{sample}/medaka/round_1_phased.vcf",
        dir = "align_samples/{sample}/medaka/variant"
    conda:
        "../envs/medaka_1_x.yaml"

    shell:
        "source ./testONT/software/medaka/bin/activate && "
        "./testONT/software/medaka/bin/medaka_variant -i {input.bam} -f {input.ref} -s {input.hdf} -p {output.dir}"


rule clair3_vfc:
    input:
        hdf = "align_samples/{sample}/medaka/consensus_probs.hdf",
        ref = REFERENCE,
        bam = "align_samples/{sample}/medaka/calls_to_draft.bam",
        depth = "align_samples/{sample}/medaka/snps.depth.gz.tbi"
    output:
        vcf = "align_samples/{sample}/medaka/clair3/merge_output.vcf.gz",
        dir = directory("align_samples/{sample}/medaka/clair3/")
    conda:
        "../envs/clair3.yaml"
    shell:
        "run_clair3.sh --bam_fn={input.bam} --ref_fn={input.ref} --threads=8 --platform='ont' --model_path=$CONDA_PREFIX/bin/models/ont --output={output.dir}"


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
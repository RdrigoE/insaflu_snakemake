rule medaka:
    input:
        i = "samples/{sample}/nano_trimmed_reads/{sample}.trimmed.fastq.gz",
        ref = REFERENCE
    output:
        dir = directory("align_samples/{sample}/medaka"),
        out = "align_samples/{sample}/medaka/consensus.fa"
    conda:
        "../envs/medaka.yaml"
    shell:
        "medaka_consensus -i {input.i} -d {input.ref} -m r941_min_high_g360 -o {output.dir} -t 8"

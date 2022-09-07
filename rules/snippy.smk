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
        dir = directory("align_samples/{sample}/snippy/")
    conda:
        "../envs/snippy.yaml"
    params:
        "--mapqual 20 --mincov 10 --minfrac 0.51"
    shell:
        "rm -r {output.dir}|"
        "snippy --R1 {input.r1} --R2 {input.r2} --ref {input.ref} --outdir {output.dir} {params}"

rule snippy_se:
    input:
        r1="samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz",
        ref=REFERENCE
    output:
        _0  = "align_samples/{sample}/snippy/snps.depth.gz",
        _1 = "align_samples/{sample}/snippy/snps.bam",
        _2 = "align_samples/{sample}/snippy/snps.tab",
        _3 = "align_samples/{sample}/snippy/snps.consensus.fa",
        dir = directory("align_samples/{sample}/snippy/")
    conda:
        "../envs/snippy.yaml"
    params:
        "--mapqual 20 --mincov 10 -minfrac 0.51"
    shell:
        "rm -r {output.dir}|"
        "snippy --se {input.r1} --ref {input.ref} --outdir {output.dir} {params}"

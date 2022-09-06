rule snippy_pe:
    input:
        r1="results/trimmed_reads_pe/{sample}_1.trimmed.fastq.gz",
        r2="results/trimmed_reads_pe/{sample}_2.trimmed.fastq.gz",
        ref=REFERENCE
    output:
        _0  =  "results/snippy_pe/{sample}_pe/snps.depth.gz",
        _1 = "results/snippy_pe/{sample}_pe/snps.bam",
        _2 = "results/snippy_pe/{sample}_pe/snps.tab",
        _3 = "results/snippy_pe/{sample}_pe/snps.consensus.fa",
        dir = directory("results/snippy_pe/{sample}_pe/")
    conda:
        "../envs/snippy.yaml"
    params:
        "--mapqual 20 --mincov 10 --minfrac 0.51"
    shell:
        "rm -r {output.dir}|"
        "snippy --R1 {input.r1} --R2 {input.r2} --ref {input.ref} --outdir {output.dir} {params}"

rule snippy_se:
    input:
        r1="results/trimmed_reads_se/{sample}.trimmed.fastq.gz",
        ref=REFERENCE
    output:
        _0  =  "results/snippy_se/{sample}_se/snps.depth.gz",
        _1 = "results/snippy_se/{sample}_se/snps.bam",
        _2 = "results/snippy_se/{sample}_se/snps.tab",
        _3 = "results/snippy_se/{sample}_se/snps.consensus.fa",
        dir = directory("results/snippy_se/{sample}_se")
    conda:
        "../envs/snippy.yaml"
    params:
        "--mapqual 20 --mincov 10 -minfrac 0.51"
    shell:
        "rm -r {output.dir}|"
        "snippy --se {input.r1} --ref {input.ref} --outdir {output.dir} {params}"
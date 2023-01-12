
rule raw_nanostat:
    input:
        get_raw_input_ont,
    output:
        dir=directory("samples/{sample}/raw_nanostat/"),
        stats="samples/{sample}/raw_nanostat/{sample}_stats.txt",
    conda:
        "../envs/nanostat.yaml"
    log:
        "logs/raw_nanostat/{sample}.log",
    shell:
        "mkdir {output.dir} -p && NanoStat --fastq {input} --outdir {output.dir} -n {wildcards.sample}_stats.txt"


rule trimmed_nanostat:
    input:
        "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
    output:
        dir=directory("samples/{sample}/nano_trimmed_fastqc/"),
        stats="samples/{sample}/nano_trimmed_fastqc/{sample}_stats.txt",
    conda:
        "../envs/nanostat.yaml"
    log:
        "logs/nano_trimmed_fastqc/{sample}.log",
    shell:
        "NanoStat --fastq {input} --outdir {output.dir} -n {wildcards.sample}_stats.txt"

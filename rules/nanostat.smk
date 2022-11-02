def get_raw_input_ont(wildcards):
    return f"user_data/{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz"


rule raw_nanostat:
    input:
        get_raw_input_ont,
    output:
        dir=directory("samples/{sample}/raw_nanostat/"),
        stats="samples/{sample}/raw_nanostat/{sample}_stats.txt",
    conda:
        "../envs/nanostat.yaml"
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
    shell:
        "NanoStat --fastq {input} --outdir {output.dir} -n {wildcards.sample}_stats.txt"

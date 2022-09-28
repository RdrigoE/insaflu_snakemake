rule raw_nanostat_se:
    input:
        i = "user_data/{sample}.fastq.gz"
    output:
        dir=directory("samples/{sample}/raw_nanostat/"),
        name = "samples/{sample}/raw_nanostat/{sample}_stats.txt"
    conda:
        "../envs/nanostat.yaml"
    shell:
        "mkdir {output.dir} -p && NanoStat --fastq {input.i} --outdir {output.dir} -n {wildcards.sample}_stats.txt"


rule trimmed_nanostat_se:
    input:
        "samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz"
    output:
        o="samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html",
        dir=directory("samples/{sample}/trimmed_fastqc/")
    conda:
        "../envs/nanostat.yaml"
    shell:
        "NanoStat --fastq {input.i} --outdir -n {output.dir}"

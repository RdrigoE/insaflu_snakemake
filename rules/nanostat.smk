rule raw_nanostat_se:
    input:
        i = "user_data/{sample}.fastq.gz"
    output:
        dir=directory("samples/{sample}/raw_nanofilt")
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc --fastq {input.i} --outdir {output.dir}"


rule trimmed_nanostat_se:
    input:
        "samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz"
    output:
        o="samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html",
        dir=directory("samples/{sample}/trimmed_fastqc/")
    conda:
        "../envs/nanostat.yaml"
    shell:
        "fastqc --fastq {input.i} --outdir {output.dir}"

rule raw_fastqc_se:
    input:
        i = "samples/{sample}." + ext
    output:
        o = "results/raw_fastqc_se/{sample}_fastqc.html",
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input.i} -o results/raw_fastqc_se/ --nogroup"

rule raw_fastqc_pe:
    input:
        i1 = "samples/{sample}_1."+ext,
        i2 = "samples/{sample}_2."+ext
    output:
        o1 = "results/raw_fastqc_pe/{sample}_1_fastqc.html",
        o2 = "results/raw_fastqc_pe/{sample}_2_fastqc.html"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input.i1} {input.i2} -o results/raw_fastqc_pe/ --nogroup"

rule trimmed_fastqc_pe:
    input:
        i1 = "results/trimmed_reads_pe/{sample}_1.trimmed.fastq.gz",
        i2 = "results/trimmed_reads_pe/{sample}_2.trimmed.fastq.gz"
    output:
        o1 = "results/trimmed_fastqc_pe/{sample}_1.trimmed_fastqc.html",
        o2 = "results/trimmed_fastqc_pe/{sample}_2.trimmed_fastqc.html",
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input.i1} {input.i2} -o results/trimmed_fastqc_pe/ --nogroup"       

rule trimmed_fastqc_se:
    input:
        "results/trimmed_reads_se/{sample}.trimmed.fastq.gz"
    output:
        "results/trimmed_fastqc_se/{sample}.trimmed_fastqc.html"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input} -o results/trimmed_fastqc_se/ --nogroup"


rule raw_fastqc_se:
    input:
        i = "user_data/{sample}." + ext
    output:
        o = "samples/{sample}/raw_fastqc/{sample}_fastqc.html",
        dir=directory("samples/{sample}/raw_fastqc")
    conda:
        "../envs/fastqc.yaml"
    params:
        "-o {output.dir} --nogroup"
    shell:
        "fastqc {input.i} {params}"

rule raw_fastqc_pe:
    input:
        i1 = "user_data/{sample}_1."+ext,
        i2 = "user_data/{sample}_2."+ext
    output:
        o1 = "samples/{sample}/{sample}_1_fastqc.html",
        o2 = "samples/{sample}/{sample}_2_fastqc.html",
        dir=directory("samples/{sample}/raw_fastqc")
    conda:
        "../envs/fastqc.yaml"
    params:
        "-o {output.dir} --nogroup"
    shell:
        "fastqc {input.i1} {input.i2}"

rule trimmed_fastqc_pe:
    input:
        i1 = "samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        i2 = "samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz"
    output:
        o1 = "samples/{sample}/trimmed_fastqc/{sample}_1.trimmed_fastqc.html",
        o2 = "samples/{sample}/trimmed_fastqc/{sample}_2.trimmed_fastqc.html",
        dir=directory("samples/{sample}/trimmed_fastqc")
    conda:
        "../envs/fastqc.yaml"
    params:
        "-o {output.dir} --nogroup"
    shell:
        "fastqc {input.i1} {input.i2} {params}"       

rule trimmed_fastqc_se:
    input:
        "samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz"
    output:
        o="samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html",
        dir="samples/{sample}/trimmed_fastqc/"
    conda:
        "../envs/fastqc.yaml"
    params:
        "-o {output.dir} --nogroup"
    shell:
        "fastqc {input} {params}"


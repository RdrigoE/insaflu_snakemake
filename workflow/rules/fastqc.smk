
rule raw_fastqc_se:
    input:
        get_raw_input_fastq_se,
    output:
        html="samples/{sample}/raw_fastqc/{sample}_fastqc.html",
        dir=directory("samples/{sample}/raw_fastqc"),
    conda:
        "../envs/fastqc.yaml"
    params:
        "--nogroup",
    threads: 1
    resources:
        mem_mb=memory["raw_fastqc_se"],
    log:
        "logs/samples/{sample}/raw_fastqc_se.log",
    benchmark:
        "benchmark/samples/{sample}/raw_fastqc_se.tsv"
    shell:
        "fastqc {input} -o {output.dir} {params} && python3 {scripts_directory}move_fastqc_output.py {wildcards.sample}"


rule raw_fastqc_pe:
    input:
        get_raw_input_fastq_pe,
    output:
        html_1="samples/{sample}/raw_fastqc/{sample}_1_fastqc.html",
        html_2="samples/{sample}/raw_fastqc/{sample}_2_fastqc.html",
        dir=directory("samples/{sample}/raw_fastqc"),
    conda:
        "../envs/fastqc.yaml"
    params:
        "--nogroup",
    threads: 1
    resources:
        mem_mb=memory["raw_fastqc_pe"],
    log:
        "logs/samples/{sample}/raw_fastqc_pe.log",
    benchmark:
        "benchmark/samples/{sample}/raw_fastqc_pe.tsv"
    shell:
        "fastqc {input} -o {output.dir} {params} && python3 {scripts_directory}move_fastqc_output.py {wildcards.sample}"


rule trimmed_fastqc_pe:
    input:
        read_1="samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        read_2="samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",
    output:
        html_1="samples/{sample}/trimmed_fastqc/{sample}_1.trimmed_fastqc.html",
        html_2="samples/{sample}/trimmed_fastqc/{sample}_2.trimmed_fastqc.html",
        dir=directory("samples/{sample}/trimmed_fastqc"),
    conda:
        "../envs/fastqc.yaml"
    params:
        "--nogroup",
    resources:
        mem_mb=memory["trimmed_fastqc_pe"],
    log:
        "logs/samples/{sample}/trimmed_fastqc_pe.log",
    benchmark:
        "benchmark/samples/{sample}/trimmed_fastqc_pe.tsv"
    threads: 2
    shell:
        "fastqc {input.read_1} {input.read_2} -o {output.dir} {params} -t {threads}"


rule trimmed_fastqc_se:
    input:
        read="samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz",
    output:
        html="samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html",
        dir=directory("samples/{sample}/trimmed_fastqc/"),
    conda:
        "../envs/fastqc.yaml"
    params:
        "--nogroup",
    threads: 2
    resources:
        mem_mb=memory["trimmed_fastqc_se"],
    log:
        "logs/samples/{sample}/trimmed_fastqc_se.log",
    benchmark:
        "benchmark/samples/{sample}/trimmed_fastqc_se.tsv"
    shell:
        "fastqc {input.read} -o {output.dir} {params} -t {threads}"

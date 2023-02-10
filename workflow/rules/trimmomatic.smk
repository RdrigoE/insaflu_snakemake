configfile: "../config/threads.yaml"


rule trimme_reads_PE:
    input:
        get_raw_input_trimmomatic,
    output:
        read_1="samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        read_2="samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",
        read_un1="samples/{sample}/trimmed_reads/{sample}_1.untrimmed.fastq.gz",
        read_un2="samples/{sample}/trimmed_reads/{sample}_2.untrimmed.fastq.gz",
    conda:
        "../envs/trimmomatic.yaml"
    params:
        get_trimmomatic_parameters(software_parameters),
    threads: config["trimmomatic_threads"]
    log:
        "logs/samples/{sample}/trimmomatic.log",
    benchmark:
        "benchmark/samples/{sample}/trimmomatic.tsv"
    shell:
        "trimmomatic PE "
        "{input} "
        "{output.read_1} "
        "{output.read_un1} "
        "{output.read_2} "
        "{output.read_un2} "
        "-threads {threads} "
        "{params}"


rule trimme_reads_SE:
    input:
        get_raw_input_trimmomatic_se,
    output:
        "samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz",
    conda:
        "../envs/trimmomatic.yaml"
    threads: config["trimmomatic_threads"]
    params:
        get_trimmomatic_parameters(software_parameters),
    log:
        "logs/samples/{sample}/trimmomatic.log",
    benchmark:
        "benchmark/samples/{sample}/trimmomatic.tsv"
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "{input} "
        "{output} "
        "{params}"

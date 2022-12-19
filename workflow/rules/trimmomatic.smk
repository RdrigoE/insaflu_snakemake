configfile: "../config/threads.yaml"


def get_raw_input_trimmomatic_se(wildcards):
    return f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz"


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
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "{input} "
        "{output} "
        "{params}"


def get_raw_input_trimmomatic(wildcards):
    return [
        f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz",
        f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq2']}.fastq.gz",
    ]


rule trimme_reads_PE:
    input:
        get_raw_input_trimmomatic,
    threads: config["trimmomatic_threads"]
    output:
        read_1="samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        read_2="samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",
        read_un1="samples/{sample}/trimmed_reads/{sample}_1.untrimmed.fastq.gz",
        read_un2="samples/{sample}/trimmed_reads/{sample}_2.untrimmed.fastq.gz",
    conda:
        "../envs/trimmomatic.yaml"
    params:
        get_trimmomatic_parameters(software_parameters),
    shell:
        "trimmomatic PE "
        "{input} "
        "{output.read_1} "
        "{output.read_un1} "
        "{output.read_2} "
        "{output.read_un2} "
        "-threads {threads} "
        "{params}"

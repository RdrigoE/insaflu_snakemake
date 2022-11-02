configfile: "config/parameters.yaml"


def get_raw_input_ont(wildcards):
    return f"user_data/{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz"


rule nanofilt_SE:
    input:
        get_raw_input_ont,
    output:
        "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
    conda:
        "../envs/nanofilt.yaml"
    params:
        get_nanofilt_parameters(software_parameters),
    shell:
        "gunzip -c {input} | NanoFilt {params} | gzip > {output}"

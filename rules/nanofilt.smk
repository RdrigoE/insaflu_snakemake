configfile: "config/parameters.yaml"

def get_nanofilt_parameters(software_parameters):
    nanofilt_parameters = f' -q {software_parameters["QUALITY"]}' if software_parameters["QUALITY"] else ''
    nanofilt_parameters += f' -l {software_parameters["LENGTH"]}' if software_parameters["LENGTH"] else ''
    nanofilt_parameters += f' --headcrop {software_parameters["HEADCROP"]}' if software_parameters["HEADCROP"] else ''
    nanofilt_parameters += f' --tailcrop {software_parameters["TAILCROP"]}' if software_parameters["TAILCROP"] else ''
    nanofilt_parameters += f' --maxlength {software_parameters["MAXLENGTH"]}' if software_parameters["MAXLENGTH"] else ''
    return nanofilt_parameters

def get_raw_input_ont(wildcards):
    return f"user_data/{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz"

rule nanofilt_SE:
    input:
        get_raw_input_ont
    output:
        "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz"
    conda:
        "../envs/nanofilt.yaml"
    params:
        get_nanofilt_parameters(software_parameters)
    shell:
        "gunzip -c {input} | NanoFilt {params} | gzip > {output}"
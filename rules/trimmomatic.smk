configfile: "config/parameters.yaml"

def get_trimmomatic_parameters(software_parameters):
    trimmomatic_params = f'ILLUMINACLIP {software_parameters["ILLUMINACLIP"]}' if software_parameters['ILLUMINACLIP'] != None else ''
    trimmomatic_params += f' HEADCROP:{software_parameters["HEADCROP"]}' if software_parameters["HEADCROP"] != None else '' 
    trimmomatic_params += f' CROP:{software_parameters["CROP"]}' if software_parameters["CROP"] else ''
    trimmomatic_params += f' SLIDINGWINDOW:{software_parameters["SLIDINGWINDOW1"]}:{software_parameters["SLIDINGWINDOW2"]}' if software_parameters["SLIDINGWINDOW2"] else ''
    trimmomatic_params += f' LEADING:{software_parameters["LEADING"]}' if software_parameters["LEADING"] else ''
    trimmomatic_params += f' TRAILING:{software_parameters["TRAILING"]}' if software_parameters["TRAILING"] else ''
    trimmomatic_params += f' MINLEN:{software_parameters["MINLEN"]}' if software_parameters["MINLEN"] else ''
    trimmomatic_params += ' TOPHRED33'
    return trimmomatic_params
def get_raw_input_trimmomatic_se(wildcards):
    return f"user_data/{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz"


rule trimme_reads_SE:
    input:
        get_raw_input_trimmomatic_se
    output:
        "samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz"
    conda:
        "../envs/trimmomatic.yaml"
    threads: 
        config['trimmomatic_threads']
    params:
        get_trimmomatic_parameters(software_parameters)
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "{input} "
        "{output} " 
        "{params}"

def get_raw_input_trimmomatic(wildcards):
    return [f"user_data/{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz",
            f"user_data/{config_user['samples'][wildcards.sample]['fastq2']}.fastq.gz"]

rule trimme_reads_PE:
    input:
        get_raw_input_trimmomatic
    threads:
        config['trimmomatic_threads']
    output:
        o1="samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        o2="samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",
        
        o_un1="samples/{sample}/trimmed_reads/{sample}_1.untrimmed.fastq.gz",
        o_un2="samples/{sample}/trimmed_reads/{sample}_2.untrimmed.fastq.gz"
    conda:
        "../envs/trimmomatic.yaml"
    params:
        get_trimmomatic_parameters(software_parameters)
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "{input} "
        "{output.o1} "
        "{output.o_un1} "
        "{output.o2} "
        "{output.o_un2} "
        "{params}"
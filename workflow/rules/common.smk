# FASTQC
def get_raw_input_fastq_se(wildcards):
    return f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz"


def get_raw_input_fastq_pe(wildcards):
    return [
        f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz",
        f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq2']}.fastq.gz",
    ]


# TRIMMOMATIC
def get_raw_input_trimmomatic_se(wildcards):
    return f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz"


def get_raw_input_trimmomatic(wildcards):
    return [
        f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz",
        f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq2']}.fastq.gz",
    ]


# GET COVERAGE
def get_depth(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/snps.depth.gz"


# NANOSTAT
def get_raw_input_ont(wildcards):
    return f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz"


# NANOFILT
# def get_raw_input_ont(wildcards):
#     return f"{dic_directory['samples']}{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz"


# MEDAKA
def get_add_freq_medaka(software_parameters):
    params = f'{software_parameters["mincov_medaka"]} {software_parameters["min_prop_for_variant_evidence"]}'
    return params


# MAKE PROJECT
def get_consensus_project(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/{wildcards.sample}_consensus.fasta"


def get_vcf_project(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/snps.vcf"


def get_directory(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/*"


# FREEBAYES
def get_snps_freebayes(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/snps.bam"

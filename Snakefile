import pandas

class Data:
    def __init__(self,file):
        self.df = pandas.read_csv(file,na_values=None)
    
    def get_name(self):
        return list(self.df["sample_name"])
    
    def get_sample_1(self):
        return list(self.df["sample_1"])
    
    def get_sample_2(self):
        if self.df["sample_2"].isnull().values.any():
            return []
        return list(self.df["sample_2"])
    
    def is_single_end(self):
        if self.get_sample_2() == []:
            return True
        return False
    
    def is_unique_file(self):
        if len(self.get_sample_1()) == 1 and self.get_sample_2() == []: 
            return True
        return False

test = Data("./config_user/sample_info_1.csv")

import yaml

def get_output_files_se(SAMPLES, PROJECT):
    SAMPLES = SAMPLES
    PROJECT = PROJECT
    with open('config_user/config_run.yaml', 'w') as file:
        documents = yaml.dump({'samples':SAMPLES, 'project':PROJECT}, file)
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_fastqc.html", sample=SAMPLES),
        expand("samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html", sample=SAMPLES),
        expand("samples/{sample}/spades/contigs.fasta", sample=SAMPLES),
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=SAMPLES),
        expand("align_samples/{sample}/snippy/snps.consensus.fa",project=PROJECT, sample=SAMPLES),
        expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/AllConsensus.fasta", project=PROJECT),
        expand("projects/{project}/main_result/coverage/{sample}_coverage.tab", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/coverage.csv",project=PROJECT),
        expand("projects/{project}/main_result/freebayes/{sample}_var.vcf", sample=SAMPLES, project=PROJECT),
        #expand("projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/mafft/mafft.fasta", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/fasttre/tree", sample=SAMPLES, project=PROJECT), 
    )

def get_output_files_pe(SAMPLES, PROJECT):
    SAMPLES = SAMPLES
    PROJECT = PROJECT
    with open('config_user/config_run.yaml', 'w') as file:
        documents = yaml.dump({'samples':SAMPLES, 'project':PROJECT}, file)
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html", sample=SAMPLES,direction=["1","2"]), #generalziar
        expand("samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html", sample=SAMPLES,direction=["1","2"]),
        expand("samples/{sample}/spades/contigs.fasta", sample=SAMPLES),
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=SAMPLES),
        expand("align_samples/{sample}/snippy/snps.consensus.fa",project=PROJECT, sample=SAMPLES),
        expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/AllConsensus.fasta", project=PROJECT),
        expand("projects/{project}/main_result/coverage/{sample}_coverage.csv", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/coverage.csv",project=PROJECT),
        expand("projects/{project}/main_result/freebayes/{sample}_var.vcf", sample=SAMPLES, project=PROJECT),
        #expand("projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/mafft/mafft.fasta", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/fasttre/tree", sample=SAMPLES, project=PROJECT),  
    )    


import os

ext = "fastq.gz"
REFERENCE = "reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.fasta"
REFERENCE_GB = "reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.gb"
REFERENCE_GFF3 = "reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.gff3"


include: "rules/fastqc.smk"
include: "rules/trimmomatic.smk"
include: "rules/spades.smk"
include: "rules/abricate.smk"
include: "rules/snippy.smk"
include: "rules/makeproject.smk"
include: "rules/getCoverage.smk"
include: "rules/mergeCoverage.smk"
include: "rules/freebayes.smk"
include: "rules/snpeff.smk"
include: "rules/concat.smk"
include: "rules/mafft.smk"
include: "rules/move_depth.smk"
include: "rules/msa_masker.smk"
include: "rules/fasttree.smk"
include: "rules/seqret.smk"


if test.get_sample_2() == []:
    get_output = get_output_files_se
else:
    get_output = get_output_files_pe


rule all:
    input:
        get_output(test.get_name(),"insaflu_comp_1")


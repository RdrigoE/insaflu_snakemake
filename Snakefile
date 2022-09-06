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
test = Data("test.csv")


def get_output_files_se(SAMPLES):
    return(
        expand("results/raw_fastqc_se/{sample}_fastqc.html", sample=SAMPLES),
        expand("results/trimmed_fastqc_se/{sample}.trimmed_fastqc.html", sample=SAMPLES),
        expand("results/abricate/{sample}/abricate_{sample}.csv", sample=SAMPLES),
        expand("results/snippy_se/{sample}_se/snps.tab", sample=SAMPLES),
        expand("results/coverage_se/{sample}_coverage.tab", sample=SAMPLES),
        expand("results/freebayes_se/{sample}_se/var.vcf", sample=SAMPLES)

    )


def get_output_files_pe(SAMPLES):
    SAMPLES = SAMPLES
    return(
        expand("results/raw_fastqc_pe/{sample}_{direction}_fastqc.html", sample=SAMPLES,direction=["1","2"]),
        expand("results/trimmed_reads_pe/{sample}_{direction}.trimmed.fastq.gz",sample=SAMPLES,direction=["1","2"]),
        expand("results/trimmed_fastqc_pe/{sample}_{direction}.trimmed_fastqc.html", sample=SAMPLES,direction=["1","2"]),
        expand("results/snippy_pe/{sample}_pe/snps.tab", sample=SAMPLES),
        expand("results/coverage_pe/{sample}_coverage.tab", sample=SAMPLES),
        expand("results/freebayes_pe/{sample}_pe/var.vcf", sample=SAMPLES),
        "results/concat/multifile.fasta",
        "results/mafft_pe/mafft.fasta",
        "results/fasttre/tree"
    )    


import os

ext = "fastq.gz"
REFERENCE = "config/SARS_CoV_2_Wuhan_Hu_1_MN908947.fasta"
REFERENCE_GB = "config/SARS_CoV_2_Wuhan_Hu_1_MN908947.gb"



include: "rules/fastqc.smk"
include: "rules/trimmomatic.smk"
include: "rules/spades.smk"
include: "rules/abricate.smk"
include: "rules/snippy.smk"
include: "rules/getCoverage.smk"
include: "rules/freebayes.smk"
include: "rules/concat.smk"
include: "rules/mafft.smk"
include: "rules/fasttree.smk"
include: "rules/seqret.smk"


if test.get_sample_2() == []:
    get_output = get_output_files_se
else:
    get_output = get_output_files_pe


rule all:
    input:
        get_output(test.get_name())


from utils.get_gene_bank import get_genes
from utils.import_user_data import Data
from utils.import_user_data import read_yaml
from utils.get_locus import get_locus_protein, get_locus
import yaml

sample_data = Data("./config_user/sample_info_1.csv")
run_config = read_yaml('./config_user/config_user1.yaml')
locus_protein_alignment = get_locus_protein(run_config["gb_reference"],run_config["locus"])

def get_output_sample_se():
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_fastqc.html", sample=config_user['samples']),
        expand("samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html", sample=config_user['samples']),
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=config_user['samples']),
    )

def get_output_sample_pe():
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html", sample=config_user['samples'],direction=["1","2"]), #generalizar
        expand("samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html", sample=config_user['samples'],direction=["1","2"]),
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=config_user['samples']),
        )    

def get_output_files_se():
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_fastqc.html", sample=config_user['samples']),
        expand("samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html", sample=config_user['samples']),
        #expand("samples/{sample}/spades/contigs.fasta", sample=config_user['samples']),
        #expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=config_user['samples']),
        #expand("align_samples/{sample}/snippy/snps.consensus.fa",project=config_user['project'], sample=config_user['samples']),
        #expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta", sample=config_user['samples'], project=config_user['project']),
        #expand("projects/{project}/main_result/AllConsensus.fasta", project=config_user['project']),
        #expand("projects/{project}/main_result/coverage/{sample}_coverage.csv", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/coverage.csv",project=config_user['project']),
        #expand("projects/{project}/main_result/freebayes/{sample}_var.vcf", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",project=config_user['project'], seg=get_locus(run_config["gb_reference"],run_config["locus"])),
        #expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}_masked.fasta",project=config_user['project'], seg=get_locus(run_config["gb_reference"],run_config["locus"])),
        # expand("projects/{project}/main_result/depth/{sample}__{ref}.depth",sample=config_user['samples'], project=config_user['project'], ref=get_locus(run_config["gb_reference"],run_config["locus"])),        
        # expand("projects/{project}/main_result/depth/{sample}.depth",sample=config_user['samples'], project=config_user['project']),        
        expand("projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/mafft/Alignment_nt_All_concat.fasta", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/{locus_protein_alignment_file}_trans.fasta", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
        expand("projects/{project}/main_result/{locus_protein_alignment_file}_mafft.fasta", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
        expand("projects/{project}/main_result/fasttre/tree", sample=config_user['samples'], project=config_user['project']), 
        expand("projects/{project}/main_result/{locus_protein_alignment_file}_tree.tree", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
        
    )

def get_output_files_pe():
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html", sample=config_user['samples'],direction=["1","2"]), #generalizar
        expand("samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html", sample=config_user['samples'],direction=["1","2"]),
        #expand("samples/{sample}/spades/contigs.fasta", sample=config_user['samples']),
        #expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=config_user['samples']),
        #expand("align_samples/{sample}/snippy/snps.consensus.fa",project=config_user['project'], sample=config_user['samples']),
        #expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta", sample=config_user['samples'], project=config_user['project']),
        #expand("projects/{project}/main_result/AllConsensus.fasta", project=config_user['project']),
        #expand("projects/{project}/main_result/coverage/{sample}_coverage.csv", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/coverage.csv",project=config_user['project']),
        #expand("projects/{project}/main_result/freebayes/{sample}_var.vcf", sample=config_user['samples'], project=config_user['project']),
        #expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",project=config_user['project'], seg=get_locus(run_config["gb_reference"],run_config["locus"])),
        #expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}_masked.fasta",project=config_user['project'], seg=get_locus(run_config["gb_reference"],run_config["locus"])),
        # expand("projects/{project}/main_result/depth/{sample}__{ref}.depth",sample=config_user['samples'], project=config_user['project'], ref=get_locus(run_config["gb_reference"],run_config["locus"])),        
        # expand("projects/{project}/main_result/depth/{sample}.depth",sample=config_user['samples'], project=config_user['project']),        
        expand("projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/mafft/Alignment_nt_All_concat.fasta", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/{locus_protein_alignment_file}_trans.fasta", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
        expand("projects/{project}/main_result/{locus_protein_alignment_file}_mafft.fasta", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
        expand("projects/{project}/main_result/fasttre/tree", sample=config_user['samples'], project=config_user['project']), 
        expand("projects/{project}/main_result/{locus_protein_alignment_file}_tree.tree", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
        
    ) 


def prepare_run(settings):
    if settings['only_samples'] == True:
        if sample_data.get_sample_2() == []:
            return get_output_sample_se
        else:
            return get_output_sample_pe
    else:
        if sample_data.get_sample_2() == []:
            return get_output_files_se
        else:
            return get_output_files_pe


REFERENCE_GB =run_config['gb_reference'] 
REFERENCE  =run_config['fasta_reference']



get_output = prepare_run(run_config)
config_user = {'samples':sample_data.get_sample_names(), 
               'project':run_config['project_name'], 
               'locus': run_config['locus'], 
               'proteins':get_genes(run_config['gb_reference'])}


with open('config/config_run.yaml', 'w') as file:
    documents = yaml.dump(config_user, file)
    file.close()


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
include: "rules/translate.smk"
include: "rules/move_depth.smk"
include: "rules/msa_masker.smk"
include: "rules/fasttree.smk"
include: "rules/seqret.smk"
include: "rules/mafft_proteins.smk"
include: "rules/fastree_proteins.smk"

rule all:
    input:
        get_output()
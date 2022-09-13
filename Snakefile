from utils.get_gene_bank import get_genes
from utils.import_user_data import Data
from utils.import_user_data import read_yaml
import yaml

sample_data = Data("./config_user/sample_info_1.csv")
run_config = read_yaml('./config_user/config_user.yaml')


REFERENCE = "reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.fasta"
REFERENCE_GB = "reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.gb"

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
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=config_user['samples']),
        #expand("align_samples/{sample}/snippy/snps.consensus.fa",project=config_user['project'], sample=config_user['samples']),
        #expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta", sample=config_user['samples'], project=config_user['project']),
        #expand("projects/{project}/main_result/AllConsensus.fasta", project=config_user['project']),
        #expand("projects/{project}/main_result/coverage/{sample}_coverage.tab", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/coverage.csv",project=config_user['project']),
        #expand("projects/{project}/main_result/freebayes/{sample}_var.vcf", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/mafft/mafft.fasta", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}.fasta",project=config_user['project'],ref=config_user["ref"],protein=config_user["proteins"]),
        expand("projects/{project}/main_result/fasttre/tree", sample=config_user['samples'], project=config_user['project']), 
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_tree.tree", project=config_user['project'],ref=config_user["ref"],protein=config_user["proteins"]),

    )

def get_output_files_pe():
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html", sample=config_user['samples'],direction=["1","2"]), #generalizar
        expand("samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html", sample=config_user['samples'],direction=["1","2"]),
        #expand("samples/{sample}/spades/contigs.fasta", sample=config_user['samples']),
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=config_user['samples']),
        #expand("align_samples/{sample}/snippy/snps.consensus.fa",project=config_user['project'], sample=config_user['samples']),
        #expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta", sample=config_user['samples'], project=config_user['project']),
        #expand("projects/{project}/main_result/AllConsensus.fasta", project=config_user['project']),
        #expand("projects/{project}/main_result/coverage/{sample}_coverage.csv", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/coverage.csv",project=config_user['project']),
        #expand("projects/{project}/main_result/freebayes/{sample}_var.vcf", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/mafft/mafft.fasta", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}.fasta",project=config_user['project'],ref=config_user["ref"],protein=config_user["proteins"]),
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_mafft.fasta",project=config_user['project'],ref=config_user["ref"],protein=config_user["proteins"]),
        expand("projects/{project}/main_result/fasttre/tree", sample=config_user['samples'], project=config_user['project']), 
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_tree.tree", project=config_user['project'],ref=config_user["ref"],protein=config_user["proteins"]),
        
    ) 


def prepare_run(settings):
    if settings['only_samples'] == False:
        sample = False
        project_name = settings['project_name']
        if os.path.isfile(settings['fasta_reference']):
            REFERENCE = settings['fasta_reference']
        if os.path.isfile(settings['gb_reference']):
            REFERENCE_GB = settings['gb_reference']
        if sample_data.get_sample_2() == []:
            return get_output_sample_se
        else:
            return get_output_sample_pe
    else:
        if sample_data.get_sample_2() == []:
            return get_output_files_se
        else:
            return get_output_files_pe


get_output = prepare_run(run_config)
config_user = {'samples':sample_data.get_sample_names(), 'project':run_config['project_name'], 'ref':run_config['species'], 'proteins':get_genes(REFERENCE_GB)}
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
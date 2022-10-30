from utils.get_gene_bank import get_genes
from utils.import_user_data import Data
from utils.import_user_data import read_yaml
from utils.get_locus import get_locus, get_id_version
from utils.get_software_parameters import get_nanofilt_parameters, mask_regions_parameters, get_trimmomatic_parameters

import yaml
import csv
from Bio import SeqIO

class Checkpoint_Alignment_aa:
    def __init__(self,prefix,sufix,genbank_file,possible_name, coverage_file):
        self.prefix = prefix
        self.sufix = sufix
        self.genbank_file = genbank_file
        self.possible_name = possible_name
        self.coverage_file = coverage_file

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.mergeCoverage.get(**w)

        # the magic, such as it is, happens here: we create the
        # information used to expand the pattern, using arbitrary
        # Python code.

        pattern = self.get_output(self.genbank_file,self.possible_name,self.coverage_file)
        return pattern
    
    def get_output(self,genbank_file,possible_name, coverage_file):
        locus_protein = []
        valide_locus = self.get_locus_w_coverage(coverage_file)
        handle_gb = open(genbank_file)
        for record in SeqIO.parse(handle_gb, "genbank"):
            for features in record.features:
                if (features.type == 'CDS'):
                    try:
                        a = features.qualifiers["locus_tag"]
                        if int(features.qualifiers['locus_tag'][0][-5:]) in valide_locus:
                            locus_protein.append(f"{self.prefix}{int(features.qualifiers['locus_tag'][0][-5:])}/Alignment_aa_{int(features.qualifiers['locus_tag'][0][-5:])}_{features.qualifiers['gene'][0]}{self.sufix}")
                    except:
                        locus_protein.append(f"{self.prefix}{possible_name}/Alignment_aa_{possible_name}_{features.qualifiers['gene'][0]}{self.sufix}")
        handle_gb.close()
        return locus_protein


    def get_locus_w_coverage(self,coverage_file):
        with open(coverage_file, newline='') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            coverage_list = list(csv_reader)
            size = len(coverage_list[0])
            final_output = []
            for times in range(1,size):
                for sample in coverage_list:
                    print(float(sample[times]))
                    if float(sample[times]) >= 90.0:
                        final_output.append(times)
                        break
        return final_output

class Checkpoint_Seg:
    def __init__(self,prefix,sufix,genbank_file,possible_name, coverage_file):
        self.prefix = prefix
        self.sufix = sufix
        self.genbank_file = genbank_file
        self.possible_name = possible_name
        self.coverage_file = coverage_file

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.mergeCoverage.get(**w)

        # the magic, such as it is, happens here: we create the
        # information used to expand the pattern, using arbitrary
        # Python code.

        pattern = self.get_output(self.genbank_file,self.possible_name,self.coverage_file)
        return pattern
    
    def get_output(self,genbank_file,possible_name, coverage_file):
        locus_protein = []
        valide_locus = self.get_locus_w_coverage(coverage_file)
        handle_gb = open(genbank_file)
        for record in SeqIO.parse(handle_gb, "genbank"):
            for features in record.features:
                if (features.type == 'CDS'):
                    try:
                        a = features.qualifiers["locus_tag"]
                        if int(features.qualifiers['locus_tag'][0][-5:]) in valide_locus:
                            locus_protein.append(f"{self.prefix}{int(features.qualifiers['locus_tag'][0][-5:])}/Alignment_nt_{int(features.qualifiers['locus_tag'][0][-5:])}{self.sufix}")
                    except:
                        locus_protein.append(f"{self.prefix}{possible_name}/Alignment_nt_{possible_name}{self.sufix}")
        handle_gb.close()
        return locus_protein


    def get_locus_w_coverage(self,coverage_file):
        with open(coverage_file, newline='') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            coverage_list = list(csv_reader)
            size = len(coverage_list[0])
            final_output = []
            for times in range(1,size):
                for sample in coverage_list:
                    print(float(sample[times]))
                    if float(sample[times]) >= 90.0:
                        final_output.append(times)
                        break
        return final_output

sample_data = Data("./config_user/sample_info.csv")

paired_illumina,single_illumina,ont_samples, sample_info_dic=sample_data.get_options()

paired_illumina_keys = paired_illumina.keys() 
single_illumina_keys = single_illumina.keys()
ont_samples_keys = ont_samples.keys()

run_config = read_yaml('./config_user/config_user.yaml')

def get_output_sample():
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html",sample = paired_illumina.keys(),direction=["1","2"]), #generalizar
        expand("samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html",sample = paired_illumina.keys(),direction=["1","2"]),
        expand("samples/{sample}/spades_pe/contigs.fasta", sample = paired_illumina.keys()),
        expand("samples/{sample}/abricate_pe/abricate_{sample}.csv", sample = paired_illumina.keys()),
   
        expand("samples/{sample}/raw_fastqc/{sample}_fastqc.html", sample=single_illumina.keys()),
        expand("samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html", sample=single_illumina.keys()),
        expand("samples/{sample}/spades_se/contigs.fasta", sample=single_illumina.keys()),
        expand("samples/{sample}/abricate_se/abricate_{sample}.csv", sample=single_illumina.keys()),

        expand("samples/{sample}/raw_nanostat/{sample}_stats.txt", sample = ont_samples.keys()), #generalizar
        expand("samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz", sample = ont_samples.keys()),
        expand("samples/{sample}/nano_trimmed_fastqc/{sample}_stats.txt",sample = ont_samples.keys()),
        # expand("samples/{sample}/rabbitqc/rabbit.html", sample =  ont_samples.keys()),

    )


def get_output_project():
    return(
        # Analyse Illumina Sample Paired-End
        expand("samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html",sample = paired_illumina.keys(),direction=["1","2"]), #generalizar
        expand("samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html",sample = paired_illumina.keys(),direction=["1","2"]),
        expand("samples/{sample}/spades_pe/contigs.fasta", sample = paired_illumina.keys()),
        expand("samples/{sample}/abricate_pe/abricate_{sample}.csv", sample = paired_illumina.keys()),
   
        # Analyse Illumina Sample Single-End
        expand("samples/{sample}/raw_fastqc/{sample}_fastqc.html", sample=single_illumina.keys()),
        expand("samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html", sample=single_illumina.keys()),
        expand("samples/{sample}/spades_se/contigs.fasta", sample=single_illumina.keys()),
        expand("samples/{sample}/abricate_se/abricate_{sample}.csv", sample=single_illumina.keys()),
        
    
        
        # Snippy for Single and Paired End Sample
        # expand("align_samples/{sample}/snippy/depth/{seg}.depth",sample = single_illumina.keys(), seg = SEGMENTS),
        # expand("align_samples/{sample}/snippy/snippy_align_{seg}.fasta",sample = single_illumina.keys(), seg = SEGMENTS),
        # expand("align_samples/{sample}/snippy/snippy_aligned_{seg}.fasta",sample = single_illumina.keys(), seg = SEGMENTS),
        # expand("align_samples/{sample}/snippy/consensus_aligned_{seg}.fasta",sample = single_illumina.keys(), seg = SEGMENTS),
        expand("align_samples/{sample}/snippy/{sample}_consensus.fasta",sample = single_illumina.keys(), seg = SEGMENTS),

        # expand("align_samples/{sample}/snippy/depth/{seg}.depth",sample = paired_illumina.keys(), seg = SEGMENTS),
        # expand("align_samples/{sample}/snippy/snippy_align_{seg}.fasta",sample = paired_illumina.keys(), seg = SEGMENTS),
        # expand("align_samples/{sample}/snippy/snippy_aligned_{seg}.fasta",sample = paired_illumina.keys(), seg = SEGMENTS),
        # expand("align_samples/{sample}/snippy/consensus_aligned_{seg}.fasta",sample = paired_illumina.keys(), seg = SEGMENTS),
        expand("align_samples/{sample}/snippy/{sample}_consensus.fasta",sample = paired_illumina.keys(), seg = SEGMENTS),
        # Analyse ONT Sample 
        expand("samples/{sample}/raw_nanostat/{sample}_stats.txt", sample = ont_samples.keys()), #generalizar
        expand("samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz", sample = ont_samples.keys()),
        expand("samples/{sample}/nano_trimmed_fastqc/{sample}_stats.txt",sample = ont_samples.keys()),
        expand("samples/{sample}/abricate_ont/abricate_{sample}.csv",sample = ont_samples.keys()),
        #expand("samples/{sample}/rabbitqc/rabbit.html", sample=config_user['samples']),
        expand("align_samples/{sample}/medaka/consensus.fasta", sample = ont_samples.keys()),
        expand("align_samples/{sample}/medaka/snps.depth.gz", sample = ont_samples.keys()),
        expand("align_samples/{sample}/medaka/{seg}.depth", sample = ont_samples.keys(), seg = SEGMENTS),
        expand("align_samples/{sample}/medaka/snps.depth.gz.tbi", sample = ont_samples.keys()),
        expand("align_samples/{sample}/medaka/round_1.vcf", sample = ont_samples.keys()),
        expand("align_samples/{sample}/medaka/snps.vcf", sample = ont_samples.keys()),
        expand("align_samples/{sample}/medaka/snps.vcf.gz", sample = ont_samples.keys()),
        expand("align_samples/{sample}/medaka/medaka_align_{seg}.fasta",sample =  ont_samples.keys(), seg=SEGMENTS),
        expand("align_samples/{sample}/medaka/medaka_aligned_{seg}.fasta", sample = ont_samples.keys(), seg=SEGMENTS),
        expand("align_samples/{sample}/medaka/consensus_aligned_{seg}.fasta",sample =  ont_samples.keys(), seg=SEGMENTS),
        expand("align_samples/{sample}/medaka/{sample}_consensus.fasta",sample =  ont_samples.keys()),
                
        # Run project

        expand("projects/{project}/main_result/coverage.csv",project=config_user['project']),
        expand("projects/{project}/main_result/coverage_translate.csv",project=config_user['project']),
        expand("projects/{project}/main_result/depth/{sample}__{ref}.depth",sample=config_user['samples'], project=config_user['project'], ref=get_locus(run_config["gb_reference"])),        
        expand("projects/{project}/main_result/validated_minor_iSNVs.csv",project=config_user['project']),
        expand("projects/{project}/main_result/validated_variants.csv",project=config_user['project']),
        expand("projects/{project}/main_result/validated_minor_iSNVs_inc_indels.csv",project=config_user['project']),
        expand("projects/{project}/main_result/proportions_iSNVs_graph.csv",project=config_user['project']),
        # expand("projects/{project}/main_result/proportions_iSNVs_graph.png",project=config_user['project']),
        expand("projects/{project}/main_result/Alignment_nt_All.fasta", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/All_nt_only_90plus.fasta", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/AllConsensus.fasta", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/All_nt.fasta", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/All_nt.nex", sample=config_user['samples'], project=config_user['project']),

        expand("projects/{project}/main_result/AllConsensus.nex", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/Alignment_nt_All.nex", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/All_nt_only_90plus.nex", sample=config_user['samples'], project=config_user['project']),


        Checkpoint_Alignment_aa(f'projects/{run_config["project_name"]}/main_result/',"_trans.fasta",run_config['gb_reference'],run_config["locus"],f"projects/{config_user['project']}/main_result/coverage_translate.csv"),
        Checkpoint_Alignment_aa(f'projects/{run_config["project_name"]}/main_result/',"_mafft.fasta",run_config['gb_reference'],run_config["locus"],f"projects/{config_user['project']}/main_result/coverage_translate.csv"),
        Checkpoint_Alignment_aa(f'projects/{run_config["project_name"]}/main_result/',"_mafft.nex",run_config['gb_reference'],run_config["locus"],f"projects/{config_user['project']}/main_result/coverage_translate.csv"),
        Checkpoint_Alignment_aa(f'projects/{run_config["project_name"]}/main_result/',"_tree.tree",run_config['gb_reference'],run_config["locus"], f"projects/{config_user['project']}/main_result/coverage_translate.csv"),
        Checkpoint_Seg(f'projects/{run_config["project_name"]}/main_result/', "_tree.tree" ,run_config['gb_reference'],run_config["locus"], f"projects/{config_user['project']}/main_result/coverage_translate.csv"),
        # expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta", project=config_user['project'], seg = SEGMENTS),
        # expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}.nex", project=config_user['project'], seg = SEGMENTS),
        expand("projects/{project}/main_result/snp_ready.txt",project=config_user['project']),
        expand("projects/{project}/main_result/Tree_ML_All.tree", sample=config_user['samples'], project=config_user['project']), 
        expand("projects/{project}/main_result/lineage_report.csv", project=config_user['project'])
)


def prepare_run(settings):
    if settings['only_samples'] == True:
        return get_output_sample
    else:
        return get_output_project


get_output = prepare_run(run_config)

REFERENCE_GB =run_config['gb_reference'] 
REFERENCE  =run_config['fasta_reference']
x = re.findall("(?<=/)(.*?)(?=.fasta)",REFERENCE)
REFERENCE_NAME = x[0]
SEGMENTS = get_locus(REFERENCE_GB)

if len(get_locus(REFERENCE_GB))  == 1:
    version_id = get_id_version(REFERENCE_GB).split('.')
    identification = version_id[0]
    version = version_id[1]
else:
    identification = ''
    version = ''

config_user = {'samples':sample_info_dic, 
               'project':run_config['project_name'], 
               'locus': run_config['locus'], 
               'proteins':get_genes(run_config['gb_reference']),
               'identification': identification,
               'version': version,
               'sample_type': sample_data.get_sample_type()}


PROJECT_NAME = config_user["project"]

with open('config/config_run.yaml', 'w') as file:
    documents = yaml.dump(config_user, file)
    file.close()

with open('config_user/parameters.yaml') as file:
    software_parameters = yaml.load(file, Loader=yaml.FullLoader)



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
include: "rules/fasttree_proteins.smk"
include: "rules/minor_iSNVs.smk"
include: "rules/snpeff_sample.smk"
include: "rules/snp_variant_validated.smk"
include: "rules/nanostat.smk"
include: "rules/nanofilt.smk"
include: "rules/rabbitqc.smk"
include: "rules/medaka.smk"
include: "rules/pangolin.smk"


rule all:
    input:
        get_output()
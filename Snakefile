from utils.get_gene_bank import get_genes
from utils.import_user_data import Data
from utils.import_user_data import read_yaml
from utils.get_locus import get_locus, get_id_version, get_locus_and_genes
from utils.get_software_parameters import (
    get_nanofilt_parameters,
    mask_regions_parameters,
    get_trimmomatic_parameters,
)

import yaml
import csv
from Bio import SeqIO


class Checkpoint_Alignment_aa:
    def __init__(self, prefix, sufix, genbank_file, coverage_file, coverage_limit):
        self.prefix = prefix
        self.sufix = sufix
        self.genbank_file = genbank_file
        self.coverage_file = coverage_file
        self.coverage_limit = coverage_limit

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.mergeCoverage.get(**w)

        # the magic, such as it is, happens here: we create the
        # information used to expand the pattern, using arbitrary
        # Python code.

        pattern = self.get_output(
            self.genbank_file, self.coverage_file, self.coverage_limit
        )
        print("Checkpoint_Alignment_aa")
        return pattern

    def get_locus_w_coverage(self, coverage_file, genbank_file, coverage_limit):
        with open(coverage_file, newline="") as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=",")
            coverage_list = list(csv_reader)
            chrom = get_locus(genbank_file)
            final_output = []
            for idx, _ in enumerate(chrom):
                chrom_pass = False
                for sample in coverage_list:
                    if not chrom_pass:
                        if float(sample[idx + 1]) >= coverage_limit:
                            chrom_pass = True
                            final_output.append(chrom_pass)
                if not chrom_pass:
                    final_output.append(chrom_pass)
        return final_output

    def get_output(self, genbank_file, coverage_file, coverage_limit):
        locus_protein = []
        valide_locus = self.get_locus_w_coverage(
            coverage_file, genbank_file, coverage_limit
        )
        segments = get_locus_and_genes(genbank_file)
        print(segments)
        for idx, seg in enumerate(segments, start=0):
            for gene in segments[seg]:
                if valide_locus[idx]:
                    locus_protein.append(
                        f"{self.prefix}{seg}/Alignment_aa_{seg}_{gene}{self.sufix}"
                    )
        with open("out.txt", "w") as f:
            f.writelines(locus_protein)
        return locus_protein


class Checkpoint_Seg:
    def __init__(self, prefix, sufix, genbank_file, coverage_file, coverage_limit):
        self.prefix = prefix
        self.sufix = sufix
        self.genbank_file = genbank_file
        self.coverage_file = coverage_file
        self.coverage_limit = coverage_limit

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.mergeCoverage.get(**w)

        # the magic, such as it is, happens here: we create the
        # information used to expand the pattern, using arbitrary
        # Python code.

        pattern = self.get_output(
            self.genbank_file, self.coverage_file, self.coverage_limit
        )
        print("Checkpoint_Seg")
        # print(pattern)
        return pattern

    def get_locus_w_coverage(self, coverage_file, genbank_file, coverage_limit):
        with open(coverage_file, newline="") as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=",")
            coverage_list = list(csv_reader)
            chrom = get_locus(genbank_file)
            final_output = []
            for idx, _ in enumerate(chrom):
                chrom_pass = False
                for sample in coverage_list:
                    if not chrom_pass:
                        if float(sample[idx + 1]) >= coverage_limit:
                            chrom_pass = True
                            final_output.append(chrom_pass)
                if not chrom_pass:
                    final_output.append(chrom_pass)
        # print(final_output)
        return final_output

    def get_output(self, genbank_file, coverage_file, coverage_limit):
        locus_protein = []
        valide_locus = self.get_locus_w_coverage(
            coverage_file, genbank_file, coverage_limit
        )
        segments = get_locus_and_genes(genbank_file)
        print(segments)

        for idx, seg in enumerate(segments, start=0):
            for _ in segments[seg]:
                if valide_locus[idx]:
                    locus_protein.append(
                        f"{self.prefix}{seg}/Alignment_nt_{seg}{self.sufix}"
                    )
        return locus_protein


class Checkpoint_Main:
    def __init__(self, genbank_file, coverage_file, coverage_limit):
        self.genbank_file = genbank_file
        self.coverage_file = coverage_file
        self.coverage_limit = coverage_limit

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.mergeCoverage.get(**w)

        # the magic, such as it is, happens here: we create the
        # information used to expand the pattern, using arbitrary
        # Python code.

        pattern = self.get_output(
            self.genbank_file, self.coverage_file, self.coverage_limit
        )
        print("Checkpoint_Main")
        return pattern

    def get_locus_w_coverage(self, coverage_file, genbank_file, coverage_limit):
        with open(coverage_file, newline="") as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=",")
            coverage_list = list(csv_reader)
            chrom = get_locus(genbank_file)
            final_output = []
            for idx, _ in enumerate(chrom):
                chrom_pass = False
                for sample in coverage_list:
                    if not chrom_pass:
                        if float(sample[idx + 1]) >= coverage_limit:
                            chrom_pass = True
                            final_output.append(chrom_pass)
                if not chrom_pass:
                    final_output.append(chrom_pass)
        return final_output

    def get_output(self, genbank_file, coverage_file, coverage_limit):
        return_list = []
        valide_locus = self.get_locus_w_coverage(
            coverage_file, genbank_file, coverage_limit
        )
        segments = get_locus_and_genes(genbank_file)
        leave = False
        print(segments)
        for idx, seg in enumerate(segments, start=0):

            for _ in segments[seg]:
                if valide_locus[idx]:
                    output_files = [
                        expand(
                            "projects/{project}/main_result/validated_minor_iSNVs.csv",
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/validated_variants.csv",
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/validated_minor_iSNVs_inc_indels.csv",
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/proportions_iSNVs_graph.csv",
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/Alignment_nt_All.fasta",
                            sample=config_user["samples"],
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/All_nt_only_90plus.fasta",
                            sample=config_user["samples"],
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/AllConsensus.fasta",
                            sample=config_user["samples"],
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/All_nt.fasta",
                            sample=config_user["samples"],
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/All_nt.nex",
                            sample=config_user["samples"],
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/AllConsensus.nex",
                            sample=config_user["samples"],
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/Alignment_nt_All.nex",
                            sample=config_user["samples"],
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/All_nt_only_90plus.nex",
                            sample=config_user["samples"],
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
                            project=config_user["project"],
                            seg=SEGMENTS,
                        ),
                        expand(
                            "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.nex",
                            project=config_user["project"],
                            seg=SEGMENTS,
                        ),
                        expand(
                            "projects/{project}/main_result/snp_ready.txt",
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/Tree_ML_All.tree",
                            sample=config_user["samples"],
                            project=config_user["project"],
                        ),
                        expand(
                            "projects/{project}/main_result/lineage_report.csv",
                            project=config_user["project"],
                        ),
                    ]
                    for i in output_files:
                        return_list.append(i[0])
                    leave = True
                if leave:
                    break
            if leave:
                break

        return return_list


sample_data = Data("./config_user/2_test_i_se_mp.csv")

(
    paired_illumina,
    single_illumina,
    ont_samples,
    sample_info_dic,
    ALIGNER,
) = sample_data.get_options()

paired_illumina_keys = paired_illumina.keys()
single_illumina_keys = single_illumina.keys()

if ALIGNER == "snippy":
    paired_illumina_snippy = paired_illumina
    single_illumina_snippy = single_illumina
    paired_illumina_iVar = {}
    single_illumina_iVar = {}
elif ALIGNER == "iVar":
    paired_illumina_iVar = paired_illumina
    single_illumina_iVar = single_illumina
    paired_illumina_snippy = {}
    single_illumina_snippy = {}

ont_samples_keys = ont_samples.keys()

run_config = read_yaml("./config_user/2_test_i_se_mp.yaml")

REFERENCE_GB = run_config["gb_reference"]
REFERENCE = run_config["fasta_reference"]
x = re.findall("(?<=/)(.*?)(?=.fasta)", REFERENCE)
REFERENCE_NAME = x[0]
SEGMENTS = get_locus(REFERENCE_GB)


def get_output_sample():
    return (
        expand(
            "samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html",
            sample=paired_illumina.keys(),
            direction=["1", "2"],
        ),  # generalizar
        expand(
            "samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html",
            sample=paired_illumina.keys(),
            direction=["1", "2"],
        ),
        expand(
            "samples/{sample}/spades_pe/contigs.fasta", sample=paired_illumina.keys()
        ),
        expand(
            "samples/{sample}/abricate_pe/abricate_{sample}.csv",
            sample=paired_illumina.keys(),
        ),
        expand(
            "samples/{sample}/abricate_pe/abricate_{sample}.yaml",
            sample=paired_illumina.keys(),
        ),
        expand(
            "samples/{sample}/raw_fastqc/{sample}_fastqc.html",
            sample=single_illumina.keys(),
        ),
        expand(
            "samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html",
            sample=single_illumina.keys(),
        ),
        expand(
            "samples/{sample}/spades_se/contigs.fasta", sample=single_illumina.keys()
        ),
        expand(
            "samples/{sample}/abricate_se/abricate_{sample}.csv",
            sample=single_illumina.keys(),
        ),
        expand(
            "samples/{sample}/abricate_se/abricate_{sample}.yaml",
            sample=single_illumina.keys(),
        ),
        expand(
            "samples/{sample}/raw_nanostat/{sample}_stats.txt",
            sample=ont_samples.keys(),
        ),  # generalizar
        expand(
            "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
            sample=ont_samples.keys(),
        ),
        expand(
            "samples/{sample}/nano_trimmed_fastqc/{sample}_stats.txt",
            sample=ont_samples.keys(),
        ),
        expand(
            "samples/{sample}/abricate_ont/abricate_{sample}.csv",
            sample=ont_samples.keys(),
        ),
        expand(
            "samples/{sample}/abricate_ont/abricate_{sample}.yaml",
            sample=ont_samples.keys(),
        ),
        # expand("samples/{sample}/rabbitqc/rabbit.html", sample =  ont_samples.keys()),
    )


def get_output_project():
    return (
        # Analyse Illumina Sample Paired-End
        expand(
            "samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html",
            sample=paired_illumina.keys(),
            direction=["1", "2"],
        ),  # generalizar
        expand(
            "samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html",
            sample=paired_illumina.keys(),
            direction=["1", "2"],
        ),
        # Trying to get the same consensus
        expand(
            "samples/{sample}/spades_pe/contigs.fasta", sample=paired_illumina.keys()
        ),
        expand(
            "samples/{sample}/abricate_pe/abricate_{sample}.csv",
            sample=paired_illumina.keys(),
        ),
        # Analyse Illumina Sample Single-End
        expand(
            "samples/{sample}/raw_fastqc/{sample}_fastqc.html",
            sample=single_illumina.keys(),
        ),
        expand(
            "samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html",
            sample=single_illumina_iVar.keys(),
        ),
        # Trying to get the same consensus
        expand(
            "samples/{sample}/spades_se/contigs.fasta", sample=single_illumina.keys()
        ),
        expand(
            "samples/{sample}/abricate_se/abricate_{sample}.csv",
            sample=single_illumina.keys(),
        ),
        expand(
            "align_samples/{sample}/iVar/{sample}_consensus.fasta",
            sample=single_illumina_iVar.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/iVar/snps_filtered.tsv",
            sample=single_illumina_iVar.keys(),
            seg=SEGMENTS,
        ),
        # Snippy for Single and Paired End Sample
        expand(
            "align_samples/{sample}/snippy/depth/{seg}.depth",
            sample=single_illumina.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/snippy/snippy_align_{seg}.fasta",
            sample=single_illumina.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/snippy/snippy_aligned_{seg}.fasta",
            sample=single_illumina.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/snippy/consensus_aligned_{seg}.fasta",
            sample=single_illumina.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/snippy/{sample}_consensus.fasta",
            sample=single_illumina_snippy.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/snippy/depth/{seg}.depth",
            sample=paired_illumina.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/snippy/snippy_align_{seg}.fasta",
            sample=paired_illumina.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/snippy/snippy_aligned_{seg}.fasta",
            sample=paired_illumina.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/snippy/consensus_aligned_{seg}.fasta",
            sample=paired_illumina.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/snippy/{sample}_consensus.fasta",
            sample=paired_illumina.keys(),
            seg=SEGMENTS,
        ),
        # Analyse ONT Sample
        expand(
            "samples/{sample}/raw_nanostat/{sample}_stats.txt",
            sample=ont_samples.keys(),
        ),  # generalizar
        expand(
            "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
            sample=ont_samples.keys(),
        ),
        expand(
            "samples/{sample}/nano_trimmed_fastqc/{sample}_stats.txt",
            sample=ont_samples.keys(),
        ),
        # Trying to get the same consensus
        expand(
            "samples/{sample}/abricate_ont/abricate_{sample}.csv",
            sample=ont_samples.keys(),
        ),
        # expand("samples/{sample}/rabbitqc/rabbit.html", sample=config_user["samples"]),
        expand(
            "align_samples/{sample}/medaka/consensus.fasta", sample=ont_samples.keys()
        ),
        expand(
            "align_samples/{sample}/medaka/snps.depth.gz", sample=ont_samples.keys()
        ),
        expand(
            "align_samples/{sample}/medaka/{seg}.depth",
            sample=ont_samples.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/medaka/snps.depth.gz.tbi", sample=ont_samples.keys()
        ),
        expand("align_samples/{sample}/medaka/round_1.vcf", sample=ont_samples.keys()),
        expand("align_samples/{sample}/medaka/snps_ann.vcf", sample=ont_samples.keys()),
        expand(
            "align_samples/{sample}/medaka/snps_ann.vcf.gz", sample=ont_samples.keys()
        ),
        expand(
            "align_samples/{sample}/medaka/medaka_align_{seg}.fasta",
            sample=ont_samples.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/medaka/medaka_aligned_{seg}.fasta",
            sample=ont_samples.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/medaka/consensus_aligned_{seg}.fasta",
            sample=ont_samples.keys(),
            seg=SEGMENTS,
        ),
        expand(
            "align_samples/{sample}/medaka/{sample}_consensus.fasta",
            sample=ont_samples.keys(),
        ),
        # Run project
        expand(
            "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
            sample=config_user["samples"],
            project=config_user["project"],
        ),
        expand(
            "projects/{project}/main_result/coverage.csv",
            project=config_user["project"],
        ),
        expand(
            "projects/{project}/main_result/coverage_translate.csv",
            project=config_user["project"],
        ),
        expand(
            "projects/{project}/main_result/depth/{sample}__{ref}.depth",
            sample=config_user["samples"],
            project=config_user["project"],
            ref=get_locus(REFERENCE_GB),
        ),
        # Trying to get the same consensus
        # Trying to get the same consensus
        Checkpoint_Alignment_aa(
            f'projects/{config_user["project"]}/main_result/',
            "_trans.fasta",
            run_config["gb_reference"],
            f"projects/{config_user['project']}/main_result/coverage_translate.csv",
            software_parameters["min_coverage_consensus"],
        ),
        Checkpoint_Alignment_aa(
            f'projects/{config_user["project"]}/main_result/',
            "_mafft.fasta",
            run_config["gb_reference"],
            f"projects/{config_user['project']}/main_result/coverage_translate.csv",
            software_parameters["min_coverage_consensus"],
        ),
        Checkpoint_Alignment_aa(
            f'projects/{config_user["project"]}/main_result/',
            "_mafft.nex",
            run_config["gb_reference"],
            f"projects/{config_user['project']}/main_result/coverage_translate.csv",
            software_parameters["min_coverage_consensus"],
        ),
        Checkpoint_Alignment_aa(
            f'projects/{config_user["project"]}/main_result/',
            "_tree.tree",
            run_config["gb_reference"],
            f"projects/{config_user['project']}/main_result/coverage_translate.csv",
            software_parameters["min_coverage_consensus"],
        ),
        Checkpoint_Seg(
            f'projects/{config_user["project"]}/main_result/',
            "_tree.tree",
            run_config["gb_reference"],
            f"projects/{config_user['project']}/main_result/coverage_translate.csv",
            software_parameters["min_coverage_consensus"],
        ),
        Checkpoint_Main(
            run_config["gb_reference"],
            f"projects/{config_user['project']}/main_result/coverage_translate.csv",
            software_parameters["min_coverage_consensus"],
        ),
    )


def prepare_run(settings):
    if settings["only_samples"] == True:
        return get_output_sample
    else:
        return get_output_project


get_output = prepare_run(run_config)


if len(get_locus(REFERENCE_GB)) == 1:
    version_id = get_id_version(REFERENCE_GB).split(".")
    identification = version_id[0]
    version = version_id[1]
else:
    identification = ""
    version = ""

config_user = {
    "samples": sample_info_dic,
    "project": run_config["project_name"],
    "locus": get_locus(run_config["gb_reference"]),
    "proteins": get_genes(run_config["gb_reference"]),
    "identification": identification,
    "version": version,
    "sample_type": sample_data.get_sample_type(),
}


PROJECT_NAME = config_user["project"]

with open("config/config_run.yaml", "w") as file:
    documents = yaml.dump(config_user, file)
    file.close()

with open("config_user/parameters.yaml") as file:
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
include: "rules/iVar.smk"


rule all:
    input:
        get_output(),

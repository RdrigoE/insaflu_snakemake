import csv
from Bio import SeqIO


def get_locus(genbank_file):
    """
    The get_locus function takes a genbank file and returns the locus number of that record.
    If there is only one record in the genbank file, it will return the possible_name. If there are multiple records,
    it will return a list of all locus numbers.

    :param genbank_file: Open the genbank file, and then parse it using seqio
    :param possible_name: Determine if the file is a single genbank file or not
    :return: A list of locus numbers
    :doc-author: Trelent
    """
    locus = []
    handle_gb = open(genbank_file)
    for record in SeqIO.parse(handle_gb, "genbank"):
        locus.append(record.name)
    return locus


def get_locus_and_genes(genbank_file: str) -> dict[str, list[str]]:
    """
    The get_locus_and_genes function takes a genbank file as input and returns a dictionary of locus names (keys)
    and gene names (values). The function is used to create the locus_gene dictionary that will be used in the
    get_locus_and_genes function. This function is called by get_locus2gene.

    :param genbank_file:str: Specify the file that contains the information about the genes in a genome
    :return: A dictionary with locus names as keys and a list of gene names as values
    """
    locus_gene: dict[str, list[str]] = {}
    handle_gb = open(genbank_file, encoding="utf-8")
    for record in SeqIO.parse(handle_gb, "genbank"):
        locus_gene[record.name] = []
        for feat in record.features:
            if feat.type == "CDS":
                locus_gene[record.name].append(feat.qualifiers["gene"][0])
    return locus_gene


def get_locus_w_coverage(coverage_file, genbank_file, coverage_limit):
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
    print(final_output)
    return final_output


def get_output(genbank_file, coverage_file, coverage_limit):
    locus_protein = []
    valide_locus = get_locus_w_coverage(coverage_file, genbank_file, coverage_limit)
    segments = get_locus_and_genes(genbank_file)
    for idx, seg in enumerate(segments, start=0):
        for _ in segments[seg]:
            if valide_locus[idx]:
                locus_protein.append(
                    f"{self.prefix}{seg}/Alignment_nt_{seg}{self.sufix}"
                )
    return locus_protein


print(
    get_output(
        genbank_file="/home/reusebio/tese/insaflu_snakemake/reference/SARS_CoV_2_COVID_19_Wuhan_Hu_1_MN908947.gbk",
        coverage_file="/home/reusebio/tese/insaflu_snakemake/projects/test_ont_covid/main_result/coverage_translate.csv",
        coverage_limit=90,
    )
)

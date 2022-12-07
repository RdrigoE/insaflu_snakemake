import sys
import csv
from Bio import SeqIO
import extract_gb_info as ggb
from yaml_io import read_yaml


def write_fasta(dictionary, filename):
    """
    The write_fasta function takes a dictionary and a filename as input.
    It then writes the sequences in the dictionary to that file in FASTA format.

    :param dictionary: Store the key value pairs of the sequence names and sequences
    :param filename: Specify the name of the file that will be created
    :return: The name of the file created
    :doc-author: Trelent
    """
    import textwrap

    with open(filename, "w") as fasta:
        for key, value in dictionary.items():
            fasta.write(f">{key}\n")
            fasta.write("\n".join(textwrap.wrap(str(value), 70)))
            fasta.write("\n")


def get_ref_adjusted_positions(alignment, positions, locus, gene):
    """
    The get_ref_adjusted_positions function takes a fasta file of aligned sequences and the positions
    of interest for each gene, and returns the adjusted positions in reference to the first sequence.


    :param alignment: Get the reference sequence
    :param positions: Get the positions of each gene group in the alignment
    :param locus: Identify the gene that is being used to get the positions
    :param gene: Determine which gene is being analyzed
    :return: The positions of the genes in the reference sequence
    :doc-author: Trelent
    """
    references = []
    new_positions = []
    ref = list(SeqIO.parse(alignment, "fasta"))[0]

    index_list = []
    for idx in range(len(ref)):
        if ref[idx] != "-":
            index_list.append(idx)
    for gene_group in positions:
        if gene_group[0] == gene:
            for group in gene_group[1]:
                group[0] = index_list[group[0]]
                if group[1] >= len(
                    index_list
                ):  # if this is == it will throw an error cuz index list é mais pequena que o numero no group[1 ]
                    group[1] = index_list[-1]
                else:
                    # print(len(index_list),group[1])
                    group[1] = index_list[group[1]]
            new_positions.append(gene_group)
    return new_positions


def get_coverage_to_translate_matrix(filename):
    """
    The get_coverage_to_translate_matrix function takes a csv file as an argument and returns a dictionary.
    The csv file is expected to have two columns, the first being the name of the gene and each subsequent column
    being one sample's coverage for that gene. The function then creates a dictionary where each key is one of
    the samples' names and its value is another dictionary with keys being genes (from the first column) and values
    being that sample's coverage for that gene.

    :param filename: Specify the name of the file to be read
    :return: A dictionary with the coverage as keys and a list of the number of times each amino acid was translated for that coverage
    :doc-author: Trelent
    """
    with open(filename, newline="") as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=",")
        coverage_dic = {}
        for row in csv_reader:
            coverage_dic[row[0]] = row[1:]
    return coverage_dic


def get_position_in_list(reference, locus):
    list_of_locus = ggb.get_locus(reference)
    return list_of_locus.index(locus)


def write_fast_aa(reference, alignment, output, locus, gene, coverage):
    """
    The write_fast_aa function takes a reference sequence, an alignment file, and an output file as arguments.
    It then parses the alignment to find all of the positions that are not gaps in the reference sequence.
    For each position it finds that is not a gap it translates those bases into amino acids and writes them to
    the output fasta file.

    :param reference: Get the reference sequence
    :param alignment: Specify the alignment file
    :param output: Specify the output file name
    :param locus: Specify the reference sequence number
    :param gene: Specify which gene to write the consensus sequence for
    :param coverage: Filter out sequences that have a coverage below 90% for the locus of interest
    :return: A dictionary with the amino acid sequences of each gene in a reference genome
    :doc-author: Trelent
    """
    position = get_position_in_list(reference, locus)
    coverage_dic = get_coverage_to_translate_matrix(coverage)
    reference_id = str(locus)

    dic_directory = read_yaml("../config/constants.yaml")

    software_parameters = read_yaml(dic_directory["software_parameters"])

    coverage_value = software_parameters["min_coverage_consensus"]
    positions = ggb.get_positions_gb(reference)
    positions = get_ref_adjusted_positions(alignment, positions, locus, gene)
    new_consensus = {}
    for gene in positions:
        new_consensus[gene[0]] = {}
        for pos in gene[1]:
            # print(f"This is {gene[0]} with the pos {pos[0], pos[1]}") #This is orf1ab with the pos (265, 13468)
            for record in SeqIO.parse(alignment, "fasta"):
                identifier = record.id
                if (
                    identifier != reference_id
                    and float(
                        coverage_dic[
                            identifier[: identifier.index(f"__{reference_id}")]
                        ][position]
                    )
                    >= coverage_value
                ):
                    try:
                        new_consensus[gene[0]][record.id] += (
                            record.seq[pos[0] : pos[1]]
                            .replace("-", "")
                            .translate(table=11, to_stop=False)
                        )
                    except:
                        new_consensus[gene[0]][record.id] = (
                            record.seq[pos[0] : pos[1]]
                            .replace("-", "")
                            .translate(table=11, to_stop=False)
                        )
                elif record.id == reference_id:
                    try:
                        new_consensus[gene[0]][record.id] += (
                            record.seq[pos[0] : pos[1]]
                            .replace("-", "")
                            .translate(table=11, to_stop=False)
                        )
                    except:
                        new_consensus[gene[0]][record.id] = (
                            record.seq[pos[0] : pos[1]]
                            .replace("-", "")
                            .translate(table=11, to_stop=False)
                        )
    for gene in new_consensus:
        write_fasta(new_consensus[gene], output)


if __name__ == "__main__":
    reference = sys.argv[1]
    alignment = sys.argv[2]
    output = sys.argv[3]
    locus = sys.argv[4]
    gene = sys.argv[5]
    coverage = sys.argv[6]

    write_fast_aa(reference, alignment, output, locus, gene, coverage)

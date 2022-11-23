from Bio import SeqIO


def get_positions_gb(genbank_file):
    """
    The get_positions_gb function takes a genbank file as input and returns a list of lists.
    Each sublist contains the name of the gene, and its start and end positions in that particular
    sequence.
    The get_positions_gb function is called by other functions to retrieve this information.

    :param genbank_file: Open the genbank file and parse it
    :return: A list of lists
    :doc-author: Trelent
    """
    positions = []

    handle_gb = open(genbank_file)
    for record in SeqIO.parse(handle_gb, "genbank"):
        for features in record.features:
            if features.type == "CDS":
                if features.qualifiers.get("gene", None) == None:
                    positions.append(
                        [features.qualifiers.get("locus_tag")[0], features.location]
                    )
                else:
                    positions.append(
                        [features.qualifiers.get("gene")[0], features.location]
                    )
    positions_clean = []

    for idx, gene in enumerate(positions):
        positions_clean.append([gene[0], []])
        for part in gene[1].parts:
            positions_clean[idx][1].append([int(part.start), int(part.end)])
    handle_gb.close()
    return positions_clean


def get_genes(genbank_file):
    """
    The get_genes function takes a genbank file as input and returns a list of all the genes in that file.
    The function is used to create a list of all the genes in each genome.

    :param genbank_file: Specify the name of the file that contains all of the genes in a genome
    :return: A list of all the genes in the genbank file
    :doc-author: Trelent
    """
    genes = []
    handle_gb = open(genbank_file)
    for record in SeqIO.parse(handle_gb, "genbank"):
        for features in record.features:
            if features.type == "CDS":
                try:
                    genes.append(features.qualifiers["gene"][0])
                except:
                    genes.append(features.qualifiers["locus_tag"][0])
    handle_gb.close()
    return genes

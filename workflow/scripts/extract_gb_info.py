"""Extract Genbank info"""
from Bio import SeqIO


def get_locus(genbank_file):
    """
    The get_locus function takes a genbank file and returns the locus number of that record.
    If there is only one record in the genbank file, it will return the possible_name.
    If there are multiple records, it will return a list of all locus numbers.

    :param genbank_file: Open the genbank file, and then parse it using seqio
    :param possible_name: Determine if the file is a single genbank file or not
    :return: A list of locus numbers
    :doc-author: Trelent
    """
    locus = []
    with open(genbank_file, encoding="utf-8") as handle_gb:
        for record in SeqIO.parse(handle_gb, "genbank"):
            locus.append(record.name)
    return locus


def get_locus_len(genbank_file):
    """
    The get_locus function takes a genbank file and returns the locus number of that record.
    If there is only one record in the genbank file, it will return the possible_name.
    If there are multiple records, it will return a list of all locus numbers.

    :param genbank_file: Open the genbank file, and then parse it using seqio
    :param possible_name: Determine if the file is a single genbank file or not
    :return: A list of locus numbers
    :doc-author: Trelents
    """
    locus = []
    with open(genbank_file, encoding="utf-8") as handle_gb:
        for record in SeqIO.parse(handle_gb, "genbank"):
            locus.append(len(record.seq))
    return locus


def get_locus_and_genes(genbank_file):
    """
    The get_locus_and_genes function takes a genbank file as input and returns a dictionary of
    locus tags and genes.
    The function iterates through the features in each record, checking if they are CDSs. If so,
    it checks if the gene qualifier is present or not. If it is present, then that value is used
    for the gene name; otherwise, the locus_tag
    is used instead.

    :param genbank_file:str: Specify the file that contains the genbank information
    :return: A dictionary with the locus tag as key and a list of gene names as value
    :doc-author: Trelent
    """

    locus_gene = {}
    with open(genbank_file, encoding="utf-8") as handle_gb:
        for record in SeqIO.parse(handle_gb, "genbank"):
            locus_gene[record.name] = []
            for feat in record.features:
                if feat.type == "CDS":
                    if feat.qualifiers.get("locus_tag", False):
                        locus_gene[record.name].append(
                            feat.qualifiers["locus_tag"][0].replace(" ", "_")
                        )
                    elif feat.qualifiers.get("gene", False):
                        locus_gene[record.name].append(
                            feat.qualifiers["gene"][0].replace(" ", "_")
                        )
                    elif feat.qualifiers.get("note", False):
                        locus_gene[record.name].append(
                            feat.qualifiers["note"][0].replace(" ", "_")
                        )
                    elif feat.qualifiers.get("product", False):
                        locus_gene[record.name].append(
                            feat.qualifiers["product"][0].replace(" ", "_")
                        )

    return locus_gene


def get_gb_name(genbank_file):
    with open(genbank_file, encoding="utf-8") as handle_gb:
        gbk = list(SeqIO.parse(handle_gb, "genbank"))[0]
        return gbk.name


def get_id_version(genbank_file):
    """
    The get_locus function takes a genbank file and returns the locus number of that record.
    If there is only one record in the genbank file, it will return the possible_name.
    If there are multiple records, it will return a list of all locus numbers.

    :param genbank_file: Open the genbank file, and then parse it using seqio
    :param possible_name: Determine if the file is a single genbank file or not
    :return: A list of locus numbers
    :doc-author: Trelent
    """
    with open(genbank_file, encoding="utf-8") as handle_gb:
        for line in handle_gb:
            if "VERSION" in line:
                return line.split(" ")[-1].strip()
    return ""


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
    with open(genbank_file, encoding="utf-8") as handle_gb:
        for record in SeqIO.parse(handle_gb, "genbank"):
            for feat in record.features:
                if feat.type == "CDS":
                    if feat.qualifiers.get("locus_tag", False):
                        positions.append(
                            [
                                feat.qualifiers.get("locus_tag")[0].replace(
                                    " ", "_"
                                ),
                                feat.location,
                            ]
                        )
                    elif feat.qualifiers.get("gene", False):
                        positions.append(
                            [
                                feat.qualifiers["gene"][0].replace(" ", "_"),
                                feat.location,
                            ]
                        )
                    elif feat.qualifiers.get("note", False):
                        positions.append(
                            [
                                feat.qualifiers["note"][0].replace(" ", "_"),
                                feat.location,
                            ]
                        )
                    elif feat.qualifiers.get("product", False):
                        positions.append(
                            [
                                feat.qualifiers["product"][0].replace(" ", "_"),
                                feat.location,
                            ]
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
    with open(genbank_file, encoding="utf-8") as handle_gb:
        for record in SeqIO.parse(handle_gb, "genbank"):
            for features in record.features:
                if features.type == "CDS":
                    if features.qualifiers.get("gene"):
                        genes.append(
                            features.qualifiers["gene"][0].replace(" ", "_")
                        )
                    elif features.qualifiers.get("locus_tag"):
                        genes.append(
                            features.qualifiers["locus_tag"][0].replace(
                                " ", "_"
                            )
                        )
                    elif features.qualifiers.get("note"):
                        genes.append(
                            features.qualifiers["note"][0].replace(" ", "_")
                        )
                    elif features.qualifiers.get("product"):
                        genes.append(
                            features.qualifiers["product"][0].replace(" ", "_")
                        )
        handle_gb.close()
    return genes


def get_identification_version(segments, reference_gb):
    """
    The get_identification_version function takes a list of SeqRecord objects and the path to a reference GenBank file.
    It returns the identification and version numbers for that GenBank file.

    :param segments:list: Determine if the reference is a single segment or not
    :param reference_gb:str: Specify the reference genome in genbank format
    :return: A tuple with the identification and version of the reference genome
    :doc-author: Trelent
    """
    if len(segments) == 1:
        version_id = get_id_version(reference_gb).split(".")
        identification = version_id[0] if len(version_id) > 0 else ""
        version = version_id[1] if len(version_id) > 1 else ""
        if len(version_id) == 1:
            new_id = get_gb_name(reference_gb)
            return new_id, ""
    else:
        identification = ""
        version = ""
    return identification, version


def get_identification_version_string(segments, reference_gb):
    """
    The get_identification_version function takes a list of SeqRecord objects and the path to a reference GenBank file.
    It returns the identification and version numbers for that GenBank file.

    :param segments:list: Determine if the reference is a single segment or not
    :param reference_gb:str: Specify the reference genome in genbank format
    :return: A tuple with the identification and version of the reference genome
    :doc-author: Trelent
    """
    if len(segments) == 1:
        version_id = get_id_version(reference_gb).split(".")
        identification = version_id[0] if len(version_id) > 0 else ""
        version = version_id[1] if len(version_id) > 1 else ""
        if len(version_id) == 1:
            new_id = get_gb_name(reference_gb)
            return new_id
    else:
        identification = ""
        version = ""
    return f"{identification}.{version}"

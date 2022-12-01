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
                if feat.qualifiers.get("gene", None) == None:
                    locus_gene[record.name].append(feat.qualifiers.get("locus_tag")[0])
                else:
                    locus_gene[record.name].append(feat.qualifiers["gene"][0])

    return locus_gene


def get_id_version(genbank_file):
    """
    The get_locus function takes a genbank file and returns the locus number of that record.
    If there is only one record in the genbank file, it will return the possible_name. If there are multiple records,
    it will return a list of all locus numbers.

    :param genbank_file: Open the genbank file, and then parse it using seqio
    :param possible_name: Determine if the file is a single genbank file or not
    :return: A list of locus numbers
    :doc-author: Trelent
    """
    handle_gb = open(genbank_file)
    for line in handle_gb:
        if "VERSION" in line:
            return line.split(" ")[-1].strip()

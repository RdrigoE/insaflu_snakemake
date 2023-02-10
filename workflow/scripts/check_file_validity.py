"""Checks the valdity of FASTA and GBK files"""
from Bio import SeqIO


def is_fasta(filename: str) -> bool:
    """
    The is_fasta function checks if a file is in FASTA format.

    :param filename: str: Specify the file to be checked
    :return: True if the file is a fasta file, and false otherwise
    :doc-author: Trelent
    """
    with open(filename, "r", encoding="utf-8") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)


def is_genbank(filename: str) -> bool:
    """
    The is_gbk function checks if a file is in GenBank format.

    :param filename: str: Specify the file to be checked
    :return: True if the file is a genbank file, and false otherwise
    :doc-author: Trelent
    """
    with open(filename, "r", encoding="utf-8") as handle:
        gbk = SeqIO.parse(handle, "genbank")
        return any(gbk)


def same_identifiers(fasta_file: str, gbk_file: str) -> bool:
    """
    The same_identifiers function checks whether the IDs in a FASTA file match those
    in a GenBank file. It does this by opening both files and parsing them with
    the SeqIO module from Biopython. The function returns True if the lists of IDs
    are identical, and False otherwise.

    :param fasta_file: str: Specify the path to the fasta file
    :param gbk_file: str: Specify the path to the genbank file
    :return: True if the order of the ids in both files are equal
    :doc-author: Trelent
    """
    with open(fasta_file, "r", encoding="utf-8") as handle:
        fasta = SeqIO.parse(handle, "fasta")  # type: ignore
        fasta_ids: list[str] = [record.id for record in fasta]
    with open(gbk_file, "r", encoding="utf-8") as handle:
        gbk = SeqIO.parse(handle, "genbank")
        gbk_locus: list[str] = [record.name for record in gbk]
    return sorted(fasta_ids) == sorted(gbk_locus)

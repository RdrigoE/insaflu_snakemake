"""Manipulate fasta files by adding and removing reference sequence"""
from Bio import SeqIO


def detach_reference(input_file, ref):
    """
    The detach_reference function takes a fasta file and removes the first record from it,
    and saves that record to a new file. The function returns the filename of the original fasta
    file without its first sequence.

    :param input_file: Specify the file to be parsed
    :param ref: Specify the name of the file that will be created to store the reference sequence
    :return: The reference sequence
    :doc-author: Trelent
    """
    with open(input_file, "r", encoding="utf8") as handle:
        record_list = list(SeqIO.parse(handle, "fasta"))
        reference = record_list.pop(0)
    with open(input_file, "w", encoding="utf8") as handle:
        SeqIO.write(record_list, handle, "fasta")

    with open(ref, "w", encoding="utf8") as handle:
        SeqIO.write(reference, handle, "fasta")


def attach_reference(input_file, ref):
    """
    The attach_reference function takes a fasta file and attaches the reference sequence to it.
    It is used in the main function to attach the reference sequence to each of the input files.

    :param input_file: Specify the name of the file that contains the sequence records
    :param ref: Specify the reference sequence
    :return: The input file with the reference sequence at the beginning of the list
    :doc-author: Trelent
    """
    with open(ref, "r", encoding="utf8") as handle:
        reference = list(SeqIO.parse(handle, "fasta"))[0]
    with open(input_file, "r", encoding="utf8") as handle:
        records = list(SeqIO.parse(handle, "fasta"))[0]
        records.insert(0, reference)
    with open(input_file, "w", encoding="utf8") as handle:
        SeqIO.write(records, handle, "fasta")

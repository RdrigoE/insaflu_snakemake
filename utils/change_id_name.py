"""Change the id in fasta file"""
import sys
from Bio import SeqIO


def change_name(file, name):
    """
    The change_name function changes the name of a fasta file.
    :param file: Specify the file that is to be read
    :param name: Change the name of the sequence
    :return: The name of the sequence
    :doc-author: Trelent
    """
    # with open(file, "r", encoding='UTF8') as handle_fasta:
    consensus_file = SeqIO.parse(file, "fasta")
    output_list = []
    for record in consensus_file:
        record.id = name
        record.description = ""
        print(record)
        output_list.append(record)
    with open(file, "w", encoding="UTF8") as handle_fasta_out_align:
        SeqIO.write(output_list, handle_fasta_out_align, "fasta")


if __name__ == "__main__":
    FILE = sys.argv[1]
    NEW_ID = sys.argv[2]
    print(FILE, NEW_ID)

    change_name(file=FILE, name=NEW_ID)

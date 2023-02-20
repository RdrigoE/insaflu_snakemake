"""Convert primers"""
import sys
import os
import csv
from time import sleep
from Bio import SeqIO


def create_alignment(reference_fasta: str, target_fasta: str) -> str:
    """
    The create_alignment function creates an alignment file from two input fasta files.
    The function takes as input a reference fasta file and a target fasta file, concatenates the two,
    and then runs MAFFT on the resulting concatenated FASTA to create an alignment. The output of this function is
    the name of the resulting alignment FASTA.

    :param reference_fasta: str: Specify the path to the reference fasta file
    :param target_fasta: str: Specify the target fasta file that we want to align
    :return: The name of the alignment file
    :doc-author: Trelent
    """
    concat_file = "concat.fasta"

    alignment_file = "alignment.fasta"
    os.system(f"cat {reference_fasta} {target_fasta} > {concat_file}")
    os.system(
        f"mafft --thread 12 --preservecase {concat_file} > {alignment_file}"
    )
    return alignment_file


def add_one_to_str_float(value: str):
    """
    The add_one_to_str_float function takes a string representation of a float and adds one to the right most decimal place.
    It then returns the new stringified float.

    :param value: Store the value that will be passed to the function
    :return: The float number with one added to the decimal part
    :doc-author: Trelent
    """
    left, rigth = value.split(".")
    result = int(rigth) + 1 // 100 + 1
    stringify = f"{left}.{result}"
    return stringify


def generate_new_coordinates(alignment_file: str) -> list[str]:
    """
    The generate_new_coordinates function takes a fasta alignment file as input and outputs a tsv file with the new coordinates for each sequence in the alignment. The function also returns the name of the second sequence in order to be used later on.

    :param alignment_file: str: Specify the path to the alignment file
    :return: A list of two elements
    :doc-author: Trelent
    """
    identification = []

    description = []
    seqs = []
    with open(alignment_file, "r", encoding="UTF-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            identification.append(record.id)
            description.append(record.description)
            seqs.append(str(record.seq))
    pos_dict = {}
    for idx, sequence in enumerate(seqs):
        pos_dict[description[idx]] = []
        counter = 0
        for pos, nuc in enumerate(sequence):
            if nuc != "-":
                counter += 1
                pos_dict[description[idx]].append(counter)
            else:
                if pos != 0 and isinstance(
                    pos_dict[description[idx]][pos - 1], str
                ):
                    pos_dict[description[idx]].append(
                        add_one_to_str_float(
                            pos_dict[description[idx]][pos - 1]
                        )
                    )
                elif pos == 0:
                    pos_dict[description[idx]].append(str(0.1))
                else:
                    pos_dict[description[idx]].append(str(counter + 0.1))
    coordinates_file = "coordinates_file.tsv"
    with open(coordinates_file, "w", encoding="UTF-8") as handle:
        write = csv.writer(handle, delimiter="\t")
        write.writerow(["alignment_position", description[0], description[1]])
        for idx in range(len(seqs[0])):
            write.writerow(
                [
                    idx + 1,
                    pos_dict[description[0]][idx],
                    pos_dict[description[1]][idx],
                ]
            )
    return [coordinates_file, identification[1]]


def generate_new_primers(
    coordinates_file: str, reference_primers: str, target_primers: str, id: str
) -> str:
    """
    The generate_new_primers function takes a file with the coordinates of the primers,
    a reference file with all the primers and their positions in a reference genome,
    and an id. It returns a new primer file that has been modified to have all positions
    in relation to one specific genome.

    :param coordinates_file: str: Specify the file containing the coordinates of each primer
    :param reference_primers: str: Specify the path to the file containing
    :param target_primers: str: Specify the output file
    :param id: str: Identify the primer pairs in the output file
    :return: A string
    :doc-author: Trelent
    """
    dict_pos = {}
    with open(coordinates_file, "r", encoding="UTF-8") as handler:
        lines = handler.readlines()
        clean_lines = []
        for line in lines[1::]:
            clean_lines.append(line.strip().split("\t"))
            dict_pos[clean_lines[-1][1]] = clean_lines[-1][2]
    with open(reference_primers, "r", encoding="UTF-8") as handler:
        primers = handler.readlines()
        clean_primers = []
        for line in primers:
            clean_primers.append(line.strip().split("\t"))
            clean_primers[-1][0] = id
            clean_primers[-1][1] = int(str(dict_pos[clean_primers[-1][1]]))
            clean_primers[-1][2] = int(str(dict_pos[clean_primers[-1][2]]))
    with open(target_primers, "w", encoding="UTF-8") as handler:
        write = csv.writer(handler, delimiter="\t")
        write.writerows(clean_primers)
    return "new_primers"


def clean_up():
    """
    The clean_up function removes the files that are created during the
    execution of this script. This function is called at the end of main().

    :return: Nothing
    :doc-author: Trelent
    """
    os.remove("concat.fasta")
    os.remove("alignment.fasta")
    os.remove("coordinates_file.tsv")


def get_new_primers(
    reference_fasta: str,
    target_fasta: str,
    reference_primers: str,
    target_primers: str,
):
    """
    The get_new_primers function takes in a reference fasta file, a target fasta file,
    and two primer files. It then creates an alignment between the two sequences using
    the mafft program and generates new coordinates for the primers based on this alignment.
    It then uses these new coordinates to generate new primers that are specific to the target sequence.

    :param reference_fasta: str: Specify the path to the reference genome fasta file
    :param target_fasta: str: Specify the path to the target
    :param reference_primers: str: Specify the path to a file containing the reference primers
    :param target_primers: str: Specify the path to the target primers file
    :param : Specify the name of the alignment file that is generated by primer3
    :return: The name of the file containing the new primers
    :doc-author: Trelent
    """
    alignment_file: str = create_alignment(reference_fasta, target_fasta)
    sleep(2)
    coordinates_file, identification = generate_new_coordinates(alignment_file)
    generate_new_primers(
        coordinates_file, reference_primers, target_primers, identification
    )
    clean_up()


if __name__ == "__main__":
    REFERENCE_FASTA = sys.argv[1]
    TARGET_FASTA = sys.argv[2]
    REFERENCE_PRIMERS = sys.argv[3]
    TARGET_PRIMERS = sys.argv[4]
    get_new_primers(
        REFERENCE_FASTA, TARGET_FASTA, REFERENCE_PRIMERS, TARGET_PRIMERS
    )

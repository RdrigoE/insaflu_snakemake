import sys
import os
from Bio import SeqIO
import csv
from time import sleep


def create_alignment(reference_fasta: str, target_fasta: str) -> str:
    concat_file = "concat.fasta"
    alignment_file = "alignment.fasta"
    os.system(f"cat {reference_fasta} {target_fasta} > {concat_file}")
    os.system(
        f"mafft --thread 12 --preservecase {concat_file} > {alignment_file}"
    )
    return alignment_file


def add_one_to_str_float(x):
    left, rigth = str(x).split(".")
    result = int(rigth) + 1 // 100 + 1
    stringify = f"{left}.{result}"
    return stringify


def generate_new_coordinates(alignment_file: str) -> list[str]:
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
    os.remove("concat.fasta")
    os.remove("alignment.fasta")
    os.remove("coordinates_file.tsv")


def get_new_primers(
    reference_fasta: str,
    target_fasta: str,
    reference_primers: str,
    target_primers: str,
):
    alignment_file: str = create_alignment(reference_fasta, target_fasta)
    sleep(2)
    coordinates_file, id = generate_new_coordinates(alignment_file)
    generate_new_primers(
        coordinates_file, reference_primers, target_primers, id
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

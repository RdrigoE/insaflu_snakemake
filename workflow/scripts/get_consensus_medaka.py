"""tbd"""
import sys
import re
from Bio import SeqIO

if __name__ == "__main__":
    consensus_files = sys.argv[1].split(" ")
    output_file = sys.argv[2]

    final_consensus_file = []
    for consensus_file in consensus_files:
        aligned_sequence = list(SeqIO.parse(consensus_file, "fasta"))[1]
        aligned_sequence.id = f"{re.findall('(?<=align_samples/)(.*?)(?=/)',consensus_file)[0]}__{aligned_sequence.id}"
        final_consensus_file.append(aligned_sequence)
    with open(output_file, "w", encoding="utf-8") as handle_fasta_out_align:
        SeqIO.write(final_consensus_file, handle_fasta_out_align, "fasta")

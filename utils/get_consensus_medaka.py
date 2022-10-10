import sys
from Bio import SeqIO

consensus_files = sys.argv[1].split(' ')
output_file = sys.argv[2]

final_consensus_file = []

for path in consensus_files:
    final_consensus_file.append(list(SeqIO.parse(path, "fasta"))[1])


with open(output_file, "w") as handle_fasta_out_align:
    SeqIO.write(final_consensus_file, handle_fasta_out_align, "fasta")
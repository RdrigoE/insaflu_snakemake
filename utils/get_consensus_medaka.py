import sys
from Bio import SeqIO
import re


consensus_files = sys.argv[1].split(' ')
output_file = sys.argv[2]

final_consensus_file = []
for consensus_file in consensus_files:
    final_consensus_file.append(list(SeqIO.parse(consensus_file, "fasta"))[1])

for record in final_consensus_file:
    record.id = f"{re.findall('(?<=align_samples/)(.*?)(?=/)',consensus_file)[0]}__{record.id}"


with open(output_file, "w") as handle_fasta_out_align:
    SeqIO.write(final_consensus_file, handle_fasta_out_align, "fasta")
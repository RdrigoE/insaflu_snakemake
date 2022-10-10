import sys
import re
from Bio import SeqIO

consensus_file = sys.argv[1]
output_file = sys.argv[2]

final_consensus_file = []

final_consensus_file.append(list(SeqIO.parse(consensus_file, "fasta")))
for record in final_consensus_file:
    for entry in record:
        entry.id = f"{re.findall('(?<=sample_)(.*?)(?=/snippy)',consensus_file)[0]}__{entry.id}"

with open(output_file, "w") as handle_fasta_out_align:
    SeqIO.write(final_consensus_file[0], handle_fasta_out_align, "fasta")
""""""
import sys
import re
from Bio import SeqIO

CONSENSUS_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2]

final_consensus_file = []

final_consensus_file.append(list(SeqIO.parse(CONSENSUS_FILE, "fasta")))
for record in final_consensus_file:
    for entry in record:
        entry.id = f"{re.findall('(?<=sample_)(.*?)(?=/snippy)',CONSENSUS_FILE)[0]}__{entry.id}"

with open(OUTPUT_FILE, "w", encoding="UTF8") as handle_fasta_out_align:
    SeqIO.write(final_consensus_file[0], handle_fasta_out_align, "fasta")

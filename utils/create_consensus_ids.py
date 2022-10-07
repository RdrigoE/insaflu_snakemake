from Bio import SeqIO
import sys
sample_id = sys.argv[1] 
file_path = sys.argv[2] 
output = sys.argv[3]
handle_fasta = open(file_path)
fasta_file = SeqIO.parse(handle_fasta, "fasta")
new_file = []
for record in fasta_file:
    record.id = f"{sample_id}__{record.id}"
    new_file.append(record)
with open(output, "w") as my_file:
    SeqIO.write(new_file, my_file, "fasta")
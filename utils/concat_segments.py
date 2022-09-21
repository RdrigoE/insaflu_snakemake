import sys
from Bio import SeqIO
from get_locus import get_locus

alignment = sys.argv[1]
reference_gb = sys.argv[2]
species = sys.argv[3]
output = sys.argv[4]
counter = 0
set_read = 0
locus = get_locus(reference_gb,species)
new_file = []
if type(locus) == type([0]):
    loop = len(locus)
    for record in SeqIO.parse(alignment, "fasta"):
        if counter == loop:
            counter = 0
            set_read += 1 
        if counter < loop:
            try:
                new_file[set_read].seq +=  record.seq
            except:
                new_file.append(record)
                #new_file[set_read].id = record.id
        counter+=1
    SeqIO.write(new_file, output, "fasta")

SeqIO.write(list(SeqIO.parse(alignment, "fasta")), output, "fasta")
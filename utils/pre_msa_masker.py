import sys
from Bio import SeqIO

alignment =sys.argv[1]
outdir = sys.argv[2]
species = sys.argv[4]
counter = 0
loop = int(sys.argv[3])
list_of_segments = [[] for x in range(loop)]
for record in SeqIO.parse(alignment, "fasta"):
    if counter == loop:
        counter = 0
    #print(counter)
    list_of_segments[counter].append(record)
    counter+=1
if loop == 1:
    SeqIO.write(list_of_segments[0], f"{outdir}/{species}/Alignment_nt_{species}.fasta", "fasta")
else:
    for i in range(loop):
        #print(i)
        SeqIO.write(list_of_segments[i], f"{outdir}/{i+1}/Alignment_nt_{i+1}.fasta", "fasta")

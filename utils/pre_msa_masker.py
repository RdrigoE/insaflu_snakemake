import sys
from Bio import SeqIO
counter = 0
loop = len(sys.argv[2])
print(loop)
for idx, record in enumerate(SeqIO.parse(sys.argv[1], "fasta")):
    print()

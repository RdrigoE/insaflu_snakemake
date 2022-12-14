import sys
from Bio import SeqIO
from manipulate_fasta import detach_reference

if __name__ == "__main__":
    INPUT_FILE = sys.argv[1]
    REF = sys.argv[2]
    detach_reference(INPUT_FILE, REF)

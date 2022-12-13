import sys
from Bio import SeqIO
from manipulate_fasta import attach_reference

if __name__ == "__main__":
    INPUT_FILE = sys.argv[1]
    REF = sys.argv[2]
    attach_reference(INPUT_FILE, REF)

import sys
from Bio import SeqIO

def prepare_msa_masker(alignment,outdir,species,loop):
    """
    The prepare_msa_masker function takes an alignment file and a species name as input.
    It then splits the alignment into multiple files, one for each loop of the species.
    The function returns nothing.
    
    :param alignment: Define the path to the alignment file
    :param outdir: Specify the directory where the files will be saved
    :param species: Name the output files
    :param loop: Split the alignment into multiple files
    :return: A list of seqio
    :doc-author: Trelent
    """
    counter = 0
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

if __name__ == '__main__':
    alignment =sys.argv[1]
    outdir = sys.argv[2]
    loop = int(sys.argv[3])
    species = sys.argv[4]
    prepare_msa_masker(alignment,outdir,species,loop)
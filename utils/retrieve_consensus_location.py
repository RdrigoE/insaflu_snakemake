from Bio import SeqIO

for record in SeqIO.parse("projects/insaflu_comp_1/main_result/mafft/mafft_masked.fasta.fasta", "fasta"):
    print(record.id)
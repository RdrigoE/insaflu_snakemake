
#with translate
from get_gene_bank import *
positions = get_positions_gb("reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.gb")

from Bio import SeqIO
from Bio.Seq import Seq
new_consensus = {}

for gene in positions:
    new_consensus[gene[0]] = {}
    for pos in gene[1]:
        #print(f"This is {gene[0]} with the pos {pos[0], pos[1]}") #This is orf1ab with the pos (265, 13468)
        for record in SeqIO.parse("projects/insaflu_comp_1/main_result/mafft/mafft_masked.fasta", "fasta"):
            # try:
            #     new_consensus[gene[0]][record.id] += record.seq[pos[0]:pos[1]].translate()
            # except:
            new_consensus[gene[0]][record.id] = record.seq[pos[0]:pos[1]].replace('-','N').translate(table=11, gap='-',to_stop=False)

def write_fasta(dictionary, filename):
    import textwrap
    with open(filename, "w") as fasta:
        for key, value in dictionary.items():
            fasta.write(f">{key}\n")
            fasta.write("\n".join(textwrap.wrap(str(value), 70)))
            fasta.write("\n")
    print(f"{filename} out!")

for gene in new_consensus:
    write_fasta(new_consensus[gene], 'utils/testfolder/'+gene+'aa.fasta')
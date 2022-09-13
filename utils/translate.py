
import sys
from Bio import SeqIO
import get_gene_bank as ggb 

def write_fasta(dictionary, filename):
        import textwrap
        with open(filename, "w") as fasta:
            for key, value in dictionary.items():
                fasta.write(f">{key}\n")
                fasta.write("\n".join(textwrap.wrap(str(value), 70)))
                fasta.write("\n")

def write_fast_aa(reference,ref_name, consensus, outdir): 
    positions = ggb.get_positions_gb(reference)
    new_consensus = {}
    print(positions)
    for gene in positions:
        new_consensus[gene[0]] = {}
        for pos in gene[1]:
            #print(f"This is {gene[0]} with the pos {pos[0], pos[1]}") #This is orf1ab with the pos (265, 13468)
            for record in SeqIO.parse(consensus, "fasta"):
                #need to find new solution to this problem
                #the final aa alignment has no gaps due to error in this segment of the pipline
                try:
                    new_consensus[gene[0]][record.id] += record.seq[pos[0]:pos[1]].replace('-','').translate(table=11, gap='-',to_stop=False)
                except:
                    new_consensus[gene[0]][record.id] = record.seq[pos[0]:pos[1]].replace('-','').translate(table=11, gap='-',to_stop=False)
    for gene in new_consensus:
        write_fasta(new_consensus[gene], f"{outdir}/Alignment_aa_{ref_name}_{gene}.fasta")
        
if __name__ == '__main__':
    #write_fast_aa("reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.gb","SARS_CoV_2","projects/insaflu_comp_1/main_result/mafft/mafft_masked.fasta","projects/insaflu_comp_1/main_result/SARS_CoV_2")
    write_fast_aa(sys.argv[1],sys.argv[2], sys.argv[3], sys.argv[4])

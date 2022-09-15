
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

def get_ref_adjusted_positions(alignment, positions, locus, gene): 
    references = []
    new_positions = []
    ref = list(SeqIO.parse(alignment, "fasta"))[0]
  
    index_list = []
    for idx in range(len(ref)):
        if ref[idx] != "-":
            index_list.append(idx)
    for gene_group in positions:
        if gene_group[0] == gene:
            for group in gene_group[1]:
                group[0] = index_list[group[0]]
                if group[1] == len(index_list):
                    group[1] = index_list[-1]
                else:
                    group[1] = index_list[group[1]]
            new_positions.append(gene_group)
    return new_positions

def write_fast_aa(reference,alignment, output, locus, gene): 
    positions = ggb.get_positions_gb(reference)
    positions = get_ref_adjusted_positions(alignment, positions, locus, gene)
    new_consensus = {}
    for gene in positions:
        new_consensus[gene[0]] = {}
        for pos in gene[1]:
            #print(f"This is {gene[0]} with the pos {pos[0], pos[1]}") #This is orf1ab with the pos (265, 13468)
            for record in SeqIO.parse(alignment, "fasta"):
                #need to find new solution to this problem
                try:
                    new_consensus[gene[0]][record.id] += record.seq[pos[0]:pos[1]].replace('-','').translate(table=11,to_stop=False)
                except:
                    new_consensus[gene[0]][record.id] = record.seq[pos[0]:pos[1]].replace('-','').translate(table=11,to_stop=False)
    for gene in new_consensus:
        write_fasta(new_consensus[gene], output)
        
if __name__ == '__main__':
    #write_fast_aa("reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.gb","SARS_CoV_2","projects/insaflu_comp_1/main_result/mafft/mafft_masked.fasta","projects/insaflu_comp_1/main_result/SARS_CoV_2")
    reference = sys.argv[1]
    alignment = sys.argv[2]
    output = sys.argv[3]
    locus = sys.argv[4]
    gene = sys.argv[5]
    
    write_fast_aa(reference, alignment, output, locus, gene)

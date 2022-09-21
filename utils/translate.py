
import sys
from Bio import SeqIO
import get_gene_bank as ggb 
import csv
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
                if group[1] >= len(index_list): #if this is == it will throw an error cuz index list é mais pequena que o numero no group[1 ]
                    group[1] = index_list[-1]
                else:
                    print(len(index_list),group[1])
                    group[1] = index_list[group[1]]
            new_positions.append(gene_group)
    return new_positions

def get_coverage_to_translate_matrix(filename):
    with open(filename, newline='') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        coverage_dic = {}
        for row in csv_reader:
            coverage_dic[row[0]] = row[1:]
    return coverage_dic

def write_fast_aa(reference,alignment, output, locus, gene, coverage): 
    if type(locus) == type(1):
        position = locus - 1 
    else:
        position = 0
    coverage_dic = get_coverage_to_translate_matrix(coverage)
    reference_id = list(SeqIO.parse(alignment, "fasta"))[0].id
    
    positions = ggb.get_positions_gb(reference)
    positions = get_ref_adjusted_positions(alignment, positions, locus, gene)
    new_consensus = {}
    for gene in positions:
        new_consensus[gene[0]] = {}
        for pos in gene[1]:
            #print(f"This is {gene[0]} with the pos {pos[0], pos[1]}") #This is orf1ab with the pos (265, 13468)
            for record in SeqIO.parse(alignment, "fasta"):
                if record.id != reference_id and float(coverage_dic[record.id[:record.id.index('__'+reference_id)]][position]) >= 90:
                    try:
                        new_consensus[gene[0]][record.id] += record.seq[pos[0]:pos[1]].replace('-','').translate(table=11,to_stop=False)
                    except:
                        new_consensus[gene[0]][record.id] = record.seq[pos[0]:pos[1]].replace('-','').translate(table=11,to_stop=False)
                elif record.id == reference_id:
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
    coverage = sys.argv[6]
    try:
        locus = int(locus)
    except:
        locus = locus

    write_fast_aa(reference, alignment, output, locus, gene, coverage)

import sys
from Bio import SeqIO


genes = []
genbank_file = "../reference/A_H1N1pdm09_A_Michigan_45_2015.gbk"
#genbank_file = "../reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.gb"
def get_locus_protein(genbank_file,possible_name):
    locus_protein = []
    handle_gb = open(genbank_file)
    for record in SeqIO.parse(handle_gb, "genbank"):
        for features in record.features:
            if (features.type == 'CDS'):
                try:
                    a = features.qualifiers["locus_tag"]
                    locus_protein.append(f"{int(features.qualifiers['locus_tag'][0][-5:])}/Alignment_aa_{int(features.qualifiers['locus_tag'][0][-5:])}_{features.qualifiers['gene'][0]}")
                except:
                    locus_protein.append(f"{possible_name}/Alignment_aa_{possible_name}_{features.qualifiers['gene'][0]}")
    handle_gb.close()
    return locus_protein

# def get_locus(genbank_file,possible_name):
#     locus = []
#     handle_gb = open(genbank_file)
#     for record in SeqIO.parse(handle_gb, "genbank"):
#         for features in record.features:
#             if (features.type == 'CDS'):
#                 try:
#                     a = features.qualifiers["locus_tag"]
#                     locus.append(f"{int(features.qualifiers['locus_tag'][0][-5:])}")
#                 except:
#                     return possible_name
#     handle_gb.close()
#     return locus

def get_locus(genbank_file,possible_name):
    count = 1
    locus = []
    handle_gb = open(genbank_file)
    for record in SeqIO.parse(handle_gb, "genbank"):
        locus.append(count)
        count+=1
    handle_gb.close()
    if len(locus) == 1:
        return possible_name
    else:
        return locus

#when locus exists:
# ['1_PB2', '2_PB1', '3_PA', '4_HA', '5_NP', '6_NA', '7_M', '8_NS']
# ['1', '2', '3', '4', '5', '6', '7', '8']
#when locus does not exist:
# ['SARS_CoV_2_orf1ab', 'SARS_CoV_2_S', 'SARS_CoV_2_ORF3a', 'SARS_CoV_2_E', 'SARS_CoV_2_M', 'SARS_CoV_2_ORF6', 'SARS_CoV_2_ORF7a', 'SARS_CoV_2_ORF8', 'SARS_CoV_2_N', 'SARS_CoV_2_ORF10']
# SARS_CoV_2
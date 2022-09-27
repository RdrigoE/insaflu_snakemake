import sys
from Bio import SeqIO
from get_locus import get_locus
import csv
import re

def get_fasta_reference_concat(fasta_ref):
    concat_reference = []
    count = 0
    for record in SeqIO.parse(fasta_ref, "fasta"):
        if count == 0:
            concat_reference.append(record)
            start = fasta_ref.index('reference')
            end = fasta_ref.index('.fasta')
            concat_reference[0].id = fasta_ref[start+10:end]
        else:
            concat_reference[0].seq += record.seq
        count += 1
    return concat_reference[0]

def get_id_sequence_from_consensus(consensus,sample_name, value_list):
    concat_consensus = []
    count = 1

    starting_segment = value_list[0]
    for record in SeqIO.parse(consensus, "fasta"):
        if count == starting_segment:
            concat_consensus.append(record)
            concat_consensus[0].id = sample_name
        elif count in value_list:
            concat_consensus[0].seq += record.seq
        count += 1
    return concat_consensus[0]

def get_id_sequence_from_consensus_strict(consensus,sample_name, value_list):
    concat_consensus = []
    count = 1
    starting_segment = value_list[0]
    for record in SeqIO.parse(consensus, "fasta"):
        if count == starting_segment:
            concat_consensus.append(record)
            concat_consensus[0].id = sample_name
        elif count in value_list:
            concat_consensus[0].seq += record.seq
        count += 1
    return concat_consensus[0]

    


consensus = sys.argv[1].split(' ')
reference_gb = sys.argv[2]
species = sys.argv[3]
output = sys.argv[4]
coverage_file = sys.argv[5]
fasta_file = sys.argv[6]
output_only_90_plus = sys.argv[7]

# consensus = '../projects/flu_testing_BVictoria/sample_SRR10885405/snippy/SRR10885405_consensus.fasta ../projects/flu_testing_BVictoria/sample_SRR10885406/snippy/SRR10885406_consensus.fasta'.split(' ')
# reference_gb = "../reference/B_Vic_B_Brisbane_60_2008.gbk"
# species = "Flu"
# output = "../projects/flu_testing_BVictoria/main_result/All_nt.fasta"
# coverage_file = "../projects/flu_testing_BVictoria/main_result/coverage_translate.csv"
# fasta_file = '../reference/B_Vic_B_Brisbane_60_2008.fasta'
# output_only_90_plus = "../projects/flu_testing_BVictoria/main_result/All_nt_only_90plus.fasta"

with open(coverage_file, newline='') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            coverage_list = list(csv_reader)

for row in coverage_list:
    for value in range(len(row)-1, 0, -1):
        # print(row[value])
        if float(row[value]) > 90:
            row[value] = value
        else:
            row.pop(value)

coverage_dic = {}
for i in coverage_list:
    coverage_dic[i[0]] = i[1:]


fasta_reference_concat = get_fasta_reference_concat(fasta_file)
new_concat_file = [fasta_reference_concat]
new_concat_file_strict = [fasta_reference_concat]
if species == 'Flu':
    for file in consensus:
        sample_name = re.findall("(?<=sample_)(.*?)(?=/snippy/)",file)[0]
        if len(coverage_dic[sample_name]) > 0:
            new_sequence = get_id_sequence_from_consensus(file,sample_name,coverage_dic[sample_name])
            new_concat_file.append(new_sequence)
        if len(coverage_dic[sample_name]) == 8:
            new_sequence = get_id_sequence_from_consensus_strict(file,sample_name,coverage_dic[sample_name])
            new_concat_file_strict.append(new_sequence)
    SeqIO.write(new_concat_file, output_only_90_plus, "fasta")
    SeqIO.write(new_concat_file_strict, output, "fasta")

else:
    for file in consensus:
        sample_name = re.findall("(?<=/snippy/)(.*?)(?=_consensus.fasta)",file)[0]
        if len(coverage_dic[sample_name]) == 1:
            new_sequence = get_id_sequence_from_consensus_strict(file,sample_name,coverage_dic[sample_name])
            new_concat_file_strict.append(new_sequence)
    SeqIO.write(new_concat_file_strict, output_only_90_plus, "fasta")
    SeqIO.write(new_concat_file_strict, output, "fasta")

import os
import sys
from Bio import SeqIO

file_path = f"align_samples/{sys.argv[2]}/snippy/"
file_list = os.listdir(file_path)


os.system(f"mkdir projects/{sys.argv[1]}/main_result/consensus/ -p")

if "snps.aligned.fa" in file_list:
    new_file = []
    for record in SeqIO.parse(file_path+"snps.aligned.fa", "fasta"):
        record.id = f"{sys.argv[2]}__{record.id}"
        new_file.append(record)
    
    SeqIO.write(new_file, f"projects/{sys.argv[1]}/main_result/consensus/"+sys.argv[2]+"_consensus.fasta", "fasta")
    
    
    

# with open(file_path+"snps.consensus.fa") as f:
#     lines = f.readlines()
#     lines[0] = f">{sys.argv[2]}__{sys.argv[3]}\n"
#     print(lines[0])
#     with open(f"projects/{sys.argv[1]}/main_result/consensus/"+sys.argv[2]+"_consensus.fasta", "w") as w:
#         w.writelines(lines)



import os
import sys
from Bio import SeqIO


project=sys.argv[1]
sample=sys.argv[2]
ref_id = sys.argv[3]
print(ref_id)
file_path = f"align_samples/{sample}/snippy/"
file_list = os.listdir(file_path)

os.system(f"mkdir projects/{project}/main_result/consensus/ -p")

if "snps.consensus.fa" in file_list:
    new_file = []
    for record in SeqIO.parse(file_path+"snps.consensus.fa", "fasta"):
        if ref_id == 'SARS_CoV_2' or ref_id == "MonkeyPox":
            record.id = ref_id
        print(f"{sample}__{record.id}")
        record.id = f"{sample}__{record.id}"
        new_file.append(record)
    
    SeqIO.write(new_file, f"projects/{project}/main_result/consensus/"+sample+"_consensus.fasta", "fasta")
    
    
    

# with open(file_path+"snps.consensus.fa") as f:
#     lines = f.readlines()
#     lines[0] = f">{sys.argv[2]}__{sys.argv[3]}\n"
#     print(lines[0])
#     with open(f"projects/{sys.argv[1]}/main_result/consensus/"+sys.argv[2]+"_consensus.fasta", "w") as w:
#         w.writelines(lines)



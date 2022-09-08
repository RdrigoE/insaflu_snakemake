import os
import sys

file_path = f"align_samples/{sys.argv[2]}/snippy/"
file_list = os.listdir(file_path)


os.system(f"mkdir projects/{sys.argv[1]}/main_result/consensus/ -p")

if "snps.consensus.fa" in file_list:
    with open(file_path+"snps.consensus.fa") as f:
        lines = f.readlines()
        lines[0] = f">{sys.argv[2]}_{sys.argv[3]}\n"
        print(lines[0])
        with open(f"projects/{sys.argv[1]}/main_result/consensus/"+lines[0][1:-1]+"_consensus.fasta", "w") as w:
            w.writelines(lines)



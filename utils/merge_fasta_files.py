import os
import sys

#sys.argv[1] -> project name
#sys.argv[2] -> species

# Get the list of all files and directories
path = f"projects/{sys.argv[1]}/samples/"
dir_list = os.listdir(path)

os.system(f"mkdir projects/{sys.argv[1]}/concat/ -p")

for file in dir_list:
    file_path = path+file
    file_dir = os.listdir(file_path)
    if "snps.consensus.fa" in file_dir:
        with open(file_path+"/"+"snps.consensus.fa") as f:
            lines = f.readlines()
            lines[0] = f">{file}_{sys.argv[2]}\n"
            with open(f"projects/{sys.argv[1]}/concat/"+lines[0][1:-1]+"_consensus.fasta", "w") as w:
                w.writelines(lines)
os.system(f"cat projects/{sys.argv[1]}/concat/* > projects/{sys.argv[1]}/concat/multifile.fasta")



import os
import sys

#sys.argv[1] -> project name
#sys.argv[2] -> sample

# Get the list of all files and directories

project = sys.argv[1]
sample = sys.argv[2]
ref_id = sys.argv[3]


path = f"projects/{project}/sample_{sample}/snippy/"
dir_list = os.listdir(path)

os.system(f"mkdir projects/{project}/consensus/ -p")

for file in dir_list:
    file_path = path+file
    file_dir = os.listdir(file_path)
    if "snps.consensus.fa" in file_dir:
        with open(file_path+"/"+"snps.consensus.fa") as f:
            lines = f.readlines()
            lines[0] = f">{file}_{ref_id}\n"
            with open(f"projects/{project}/consensus/"+lines[0][1:-1]+"_consensus.fasta", "w") as w:
                w.writelines(lines)
            w.close()
        f.close()



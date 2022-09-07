import os
import sys
# Get the list of all files and directories
path = f"results/snippy_{sys.argv[1]}/"
dir_list = os.listdir(path)
dir_list
os.system(f"mkdir results/concat_{sys.argv[1]}/ -p")
for file in dir_list:
    file_path = path+file
    file_dir = os.listdir(file_path)
    if "snps.consensus.fa" in file_dir:
        with open(file_path+"/"+"snps.consensus.fa") as f:
            lines = f.readlines()
            lines[0] = f">{file[:-3]}_SARS_COV_2\n"
            with open(f"results/concat_{sys.argv[1]}/"+lines[0][1:-1]+"_consensus.fasta", "w") as w:
                w.writelines(lines)
os.system(f"cat results/concat_{sys.argv[1]}/* > results/concat_{sys.argv[1]}/multifile.fasta")



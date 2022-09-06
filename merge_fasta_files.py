import os
# Get the list of all files and directories
path = "results/snippy_pe/"
dir_list = os.listdir(path)
dir_list

for file in dir_list:
    file_path = path+file
    file_dir = os.listdir(file_path)
    if "snps.consensus.fa" in file_dir:
        with open(file_path+"/"+"snps.consensus.fa") as f:
            lines = f.readlines()
            lines[0] = f">{file[:-3]}_SARS_COV_2\n"
            with open("results/concat/"+lines[0][1:]+"consensus.fasta", "w") as w:
                w.writelines(lines)
os.system("cat results/concat/* > results/concat/multifile.fasta")



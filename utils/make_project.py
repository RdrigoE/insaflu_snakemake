import os
import sys
import os
def makeproject(PROJECT,SAMPLES):

    path = "samples/"
    dir_list = os.listdir(path)

    os.system(f"mkdir projects/{PROJECT}/main_result/ -p")

    for file in dir_list:
        file_path = path+file
        file_dir = os.listdir(file_path)
        if "snps.consensus.fa" in file_dir:
            with open(file_path+"/"+"snps.consensus.fa") as f:
                lines = f.readlines()
                lines[0] = f">{file}_{SAMPLES}\n"
                with open(f"projects/{PROJECT}/concat/"+lines[0][1:-1]+"_consensus.fasta", "w") as w:
                    w.writelines(lines)
    os.system(f"cat projects/{PROJECT}/concat/* > projects/{PROJECT}/concat/multifile.fasta")





    os.system(f"mkdir projects/{PROJECT}/samples -p")
    os.system(f"cp samples/{{"".join(list)}} projects/{PROJECT}/samples/")


makeproject(sys.argv[1], sys.argv[2])
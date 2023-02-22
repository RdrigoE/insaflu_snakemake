import os

with open("./SRR6849735.txt") as handle:

    lines = handle.readlines()
    last_item = "replace"
    for ref in lines[1:]:
        ref = ref.strip()
        os.system(f" sed -i 's/{last_item}/{ref}/g'  config/constants.yaml  ")
        last_item = ref
        os.system("snakemake -c 12 --use-conda --ri")

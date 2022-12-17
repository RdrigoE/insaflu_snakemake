import os

files = [
    "demo_SARS_CoV_2_INSaFLU",
    "demo_RSV_INSaFLU",
    "demo_MPXV_INSaFLU",
    "demo_influenza_INSaFLU",
]

last_file = "demo_influenza_INSaFLU"

for file in files:
    print(f"Replacing string: {file}")
    os.system(f"sed -i 's/{last_file}/{file}/' config/constants.yaml")
    os.system("snakemake -c 12 --use-conda --ri")
    last_file = file

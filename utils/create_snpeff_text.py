import sys
import re
from Bio import SeqIO
import os

snpeff_path = sys.argv[1]
ref_path = sys.argv[2]
locus = sys.argv[3]
reference_name = sys.argv[4]

if locus == 'Flu':
    text = [f"{reference_name}.chromosome : 1, 2, 3, 4, 5, 6, 7, 8\n",
    f"{reference_name}.genome : flu\n",
    f"{reference_name}.1.codonTable : Bacterial_and_Plant_Plastid\n",
    f"{reference_name}.2.codonTable : Bacterial_and_Plant_Plastid\n",
    f"{reference_name}.3.codonTable : Bacterial_and_Plant_Plastid\n",
    f"{reference_name}.4.codonTable : Bacterial_and_Plant_Plastid\n",
    f"{reference_name}.5.codonTable : Bacterial_and_Plant_Plastid\n",
    f"{reference_name}.6.codonTable : Bacterial_and_Plant_Plastid\n",
    f"{reference_name}.7.codonTable : Bacterial_and_Plant_Plastid\n",
    f"{reference_name}.8.codonTable : Bacterial_and_Plant_Plastid\n"]
else:
    version = SeqIO.parse(ref_path, "genbank")

    for i in list(version.records):
        identification = i.annotations['accessions'][-1]
        version = i.annotations['sequence_version']
    text = [f"\n{reference_name}.genome: {reference_name}\n",
    f"{reference_name}.{identification}.{version}.codonTable : Bacterial_and_Plant_Plastid\n",]
    project = sys.argv[5]
    samples = sys.argv[6]

    samples = samples.split(" ")
    for sample in samples:
        os.system(f"sed -i 's/{locus}/{identification}.{version}/g' projects/{project}/main_result/freebayes/{sample}_var.vcf")
    for sample in samples:
        os.system(f"sed -i 's/{locus}/{identification}.{version}/g' projects/{project}/sample_{sample}/snippy/snps.vcf")

with open('config/snpeff.config', 'a') as snpeff:
    snpeff.writelines(text)


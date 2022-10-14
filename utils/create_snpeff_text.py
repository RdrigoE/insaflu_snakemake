import subprocess
import sys
import re
from Bio import SeqIO
import os
from get_locus import get_id_version, get_locus

def prepare_snpeff_run(ref_path_gb,locus,ref_path_fa,reference_name,output):
    """
    The prepare_snpeff_run function creates a snpEFF database for the reference genome.
    It takes as input:
        - The path to the reference genome in GenBank format, and 
        - The name of the output file that will be created by this function. 
    
       This function creates a directory called 'config/data' if it does not already exist, 
       and then creates three files inside that directory: genes.gbk, sequences.fa, 
       and snpeff_genes_prod_snpsiftdb-passed-only-geneidmap 
       (the last one is created by running SnpSift).
    
    :param snpeff_path: Specify the path to snpeff
    :param ref_path_gb: Specify the path to the genbank file of a reference genome
    :param locus: Determine the reference genome to use for snpeff
    :param ref_path_fa: Specify the path to the reference genome fasta file
    :param reference_name: Name the reference files in the snpeff config file
    :param output: Create a file that indicates the database is ready
    :return: A text file with a message &quot;database ready!&quot;
    :doc-author: Trelent
    """
    if len(locus) > 1:
        text = [f"{reference_name}.chromosome : {' '.join(str(locus))}\n",
        f"{reference_name}.genome : flu\n"]
        for i in locus:
            text.append(f"{reference_name}.{i}.codonTable : Bacterial_and_Plant_Plastid\n")
    else:
        text = [f"\n{reference_name}.genome: {reference_name}\n",
        f"{reference_name}.{get_id_version(ref_path_gb)}.codonTable : Bacterial_and_Plant_Plastid\n",]
        # project = sys.argv[5]
        # samples = sys.argv[6]

        # samples = samples.split(" ")
        # for sample in samples:
        #     os.system(f"sed -i 's/{locus}/{identification}.{version}/g' projects/{project}/main_result/freebayes/{sample}_var.vcf")
        # for sample in samples:
        #     os.system(f"sed -i 's/{locus}/{identification}.{version}/g' projects/{project}/sample_{sample}/snippy/snps.vcf")

    with open('config/snpeff.config' , 'r') as f:
        lines = f.readlines()
        joined = ''.join(lines)
        find = joined.find(''.join(text))

    if find == -1:
        subprocess.run(f"mkdir  config/data/{reference_name} -p", shell=True)
        subprocess.run( f"cat {ref_path_gb} > config/data/{reference_name}/genes.gbk", shell=True)
        subprocess.run( f"cat {ref_path_fa} > config/data/{reference_name}/sequences.fa", shell=True)
        with open('config/snpeff.config', 'a') as snpeff:
            snpeff.writelines(text)
        
        subprocess.run( f"snpEff build -genbank {reference_name} -c config/snpeff.config > {output}", shell=True)
    else:
        with open(output, 'w') as out:
            out.write("Database Ready!")



if __name__ == '__main__':
    snpeff_path = sys.argv[1]
    ref_path_gb = sys.argv[2]
    ref_path_fa = sys.argv[3]
    reference_name = sys.argv[4]
    output = sys.argv[5]
    locus = get_locus(ref_path_gb)
    prepare_snpeff_run(ref_path_gb,locus,ref_path_fa,reference_name,output)
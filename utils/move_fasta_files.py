import os
import sys
from Bio import SeqIO

def move_fasta_files(project,sample, ref_id):
    """
    The move_fasta_files function moves the snippy fasta files to a new directory.
    It also changes the name of the file to include both sample and reference id.
    
    :param project: Specify the project directory
    :param sample: Create a new folder for each sample
    :param ref_id: Rename the reference sequence in the snippy consensus file
    :return: The list of files that were created
    :doc-author: Trelent
    """
    file_path = f"align_samples/{sample}/snippy/"
    file_list = os.listdir(file_path)

    new_file = []
    for record in SeqIO.parse(file_path+"snps.consensus.fa", "fasta"):
        if ref_id != "Flu":
            record.id = ref_id
        print(f"{sample}__{record.id}")
        record.id = f"{sample}__{record.id}"
        new_file.append(record)

    SeqIO.write(new_file, f"projects/{project}/sample_{sample}/snippy/{sample}_consensus.fasta", "fasta")

if __name__ == '__main__':
    project=sys.argv[1]
    sample=sys.argv[2]
    ref_id = sys.argv[3]
    move_fasta_files(project,sample, ref_id)
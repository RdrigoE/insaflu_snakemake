import sys
from Bio import SeqIO
from extract_gb_info import get_locus


def prepare_msa_masker(alignment, outdir, reference_gb):
    """
    The prepare_msa_masker function takes a multiple sequence alignment and a reference genbank file as input.
    It then creates an output directory for each locus in the reference genome, and writes the sequences of that locus to
    a fasta file within each output directory.

    :param alignment: Specify the path to the alignment file
    :param outdir: Specify the directory where the files will be saved
    :param reference_gb: Get the locus of each segment in the alignment
    :return: A list of lists
    :doc-author: Trelent
    """
    counter = 0
    seg_names = get_locus(reference_gb)
    n_locus = len(seg_names)
    list_of_segments = [[] for x in range(n_locus)]
    for record in SeqIO.parse(alignment, "fasta"):
        if counter == len(list_of_segments):
            counter = 0
        list_of_segments[counter].append(record)
        counter += 1
    for index, seg in enumerate(list_of_segments):
        SeqIO.write(
            seg,
            f"{outdir}/{seg_names[index]}/Alignment_nt_{seg_names[index]}.fasta",
            "fasta",
        )


if __name__ == "__main__":
    alignment = sys.argv[1]
    outdir = sys.argv[2]
    reference_gb = sys.argv[3]
    prepare_msa_masker(alignment, outdir, reference_gb)

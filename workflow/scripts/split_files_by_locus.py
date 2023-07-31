import sys
from Bio import SeqIO
from extract_gb_info import get_locus


def get_locus_out_id(string_id: str):
    if string_id.find("__") == -1:
        return False
    return string_id[len(string_id) - string_id[::-1].index("__"):]


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
    locus_names = get_locus(reference_gb)

    locus_records = {x: [] for x in locus_names}

    for locus in locus_records:
        for record in SeqIO.parse(alignment, "fasta"):
            if locus == record.id or locus == get_locus_out_id(record.id):
                locus_records[locus].append(record)
    for locus, list_of_records in locus_records.items():
        SeqIO.write(
            list_of_records,
            f"{outdir}/{locus}/Alignment_nt_{locus}.fasta",
            "fasta",
        )


if __name__ == "__main__":
    alignment = sys.argv[1]
    outdir = sys.argv[2]
    reference_gb = sys.argv[3]
    prepare_msa_masker(alignment, outdir, reference_gb)

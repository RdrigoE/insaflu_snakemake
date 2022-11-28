import csv
import sys
import subprocess
from Bio import SeqIO


class Sequence:
    def __init__(self, sequence, reference) -> None:
        self.reference = str(reference)
        self.sequence = str(sequence)
        self.ref_len = len(self.reference)
        self.seq_len = len(self.sequence)

    def __str__(self) -> str:
        return self.sequence

    def get_number_gaps(self):
        count_gaps = 0
        for entry in self.sequence:
            if entry == "-":
                count_gaps += 1
        return count_gaps

    def get_number_Ns(self):
        count_gaps = 0
        for entry in self.sequence:
            if entry == "N":
                count_gaps += 1
        return count_gaps

    def get_match_nucleotides(self):
        max_length = (
            len(self.reference)
            if len(self.reference) < len(self.sequence)
            else len(self.sequence)
        )
        num_equals = 0
        missmatched = 0
        missligned = 0
        Ns = 0
        positions = []
        diff_dic = start_dic()
        count = 0
        for nucleotide_index in range(max_length):
            if self.reference[nucleotide_index] == self.sequence[nucleotide_index]:
                num_equals += 1
            else:
                count +=1
                print(nucleotide_index,self.reference[nucleotide_index], self.sequence[nucleotide_index], count, sample)
                if (
                    self.reference[nucleotide_index] in "ATCG"
                    and self.sequence[nucleotide_index] in "ATCG"
                ):
                    missmatched += 1
                if (
                    self.reference[nucleotide_index] in "ATCG"
                    and self.sequence[nucleotide_index] == "-"
                ) or (
                    self.reference[nucleotide_index] == "-"
                    and self.sequence[nucleotide_index] in "AGTC"
                ):
                    missligned += 1
                if (
                    self.reference[nucleotide_index] in "ATCG-"
                    and self.sequence[nucleotide_index] == "N"
                ) or (
                    self.reference[nucleotide_index] == "N"
                    and self.sequence[nucleotide_index] in "AGTC-"
                ):
                    Ns += 1

                diff_dic[
                    f"{self.reference[nucleotide_index]}>{self.sequence[nucleotide_index]}"
                ] += 1
                # print(f"{nucleotide_index} =====> reference = {self.reference[nucleotide_index]} > snakemake = {self.sequence[nucleotide_index]}")

        return num_equals, diff_dic, missmatched, missligned, Ns


def get_locus(genbank_file):
    """
    The get_locus function takes a genbank file and returns the locus number of that record.
    If there is only one record in the genbank file, it will return the possible_name. If there are multiple records,
    it will return a list of all locus numbers.

    :param genbank_file: Open the genbank file, and then parse it using seqio
    :param possible_name: Determine if the file is a single genbank file or not
    :return: A list of locus numbers
    :doc-author: Trelent
    """
    locus = []
    handle_gb = open(genbank_file)
    for record in SeqIO.parse(handle_gb, "genbank"):
        locus.append(record.name)
    return locus


def start_dic():
    dic = {}
    for i in "ATCG-N":
        for j in "ATCG-N":
            if i != j:
                dic[f"{i}>{j}"] = 0
    return dic


def normalize_id(seqIO_object, reference_gb):
    if len(get_locus(reference_gb)) > 1:
        new_obj = []
        for record in seqIO_object:
            new_obj.append(record)
        return new_obj

    new_obj = []
    for record in seqIO_object:
        if "_" in record.id:
            record.id = record.id[: record.id.index("_")]
        new_obj.append(record)
    return new_obj


def get_sequence_dic(context):
    dic = {}
    for record in list(context):
        # dic[record.id[:record.id.index('_')]] = record.seq
        dic[record.id] = record.seq
    return dic


def get_object_dic(context):
    dic = {}
    for record in list(context):
        # dic[record.id[:record.id.index('_')]] = record.seq
        dic[record.id] = record
    return dic


def get_output(analyse, match_nucleotides, diff, missmatched, missligned, Ns):
    return {
        "match": match_nucleotides,
        "snakemake length": analyse.seq_len,
        "insaflu length": analyse.ref_len,
        "missmatched": missmatched,
        "missligned": missligned,
        "Ns": Ns,
        "total_differences": analyse.ref_len - match_nucleotides,
        "total_differences %": round(
            (analyse.ref_len - match_nucleotides) * 100 / analyse.ref_len, 2
        ),
        "diff": diff,
    }


reference_gb = sys.argv[3]

snakemake = sys.argv[1]
snakemake_handle = open(snakemake)
snakemake_records = SeqIO.parse(snakemake_handle, "fasta")
snakemake_records = normalize_id(snakemake_records, reference_gb)
# SeqIO.write(snakemake_records, snakemake, "fasta")

insaflu = sys.argv[2]
insaflu_handle = open(insaflu)
insaflu_records = SeqIO.parse(insaflu_handle, "fasta")
insaflu_records = normalize_id(insaflu_records, reference_gb)
# SeqIO.write(insaflu_records, insaflu, "fasta")

snakemake_obj = get_object_dic(snakemake_records)

insaflu_obj = get_object_dic(insaflu_records)
snakemake_seq = get_sequence_dic(snakemake_records)
insaflu_seq = get_sequence_dic(insaflu_records)

# Check if everything is equal between the two
results = {}
to_align = []
for sample in snakemake_seq:
    analyse = Sequence(snakemake_seq[sample], insaflu_seq[sample])
    (
        match_nucleotides,
        diff,
        missmatched,
        missligned,
        Ns,
    ) = analyse.get_match_nucleotides()
    if match_nucleotides != max(analyse.seq_len, analyse.ref_len):
        print(sample,missmatched,missligned,Ns)
        to_align.append(sample)
        # print(sample)
        # print(
        #     get_output(analyse,match_nucleotides,diff,missmatched, missligned, Ns)
        # )
    else:
        results[sample] = get_output(
            analyse, match_nucleotides, diff, missmatched, missligned, Ns
        )
        # print(sample)
        # print(
        #     get_output(analyse,match_nucleotides,diff,missmatched, missligned, Ns)
        # )

for sample in snakemake_obj:
    if sample in to_align:
        new_file_to_align_name = f"{sample}_to_align.fasta"
        aligned_file = f"{sample}_aligned.fasta"
        snakemake_sample = snakemake_obj[sample]
        snakemake_sample.id = f"{sample}_snakemake"
        snakemake_sample.description = ""
        insaflu_sample = insaflu_obj[sample]
        insaflu_sample.id = f"{sample}_insaflu"
        insaflu_sample.description = ""
        new_file_to_align = [snakemake_sample, insaflu_sample]
        SeqIO.write(new_file_to_align, new_file_to_align_name, "fasta")
        subprocess.run(
            f"mafft --thread 12 --quiet --preservecase {new_file_to_align_name} > {aligned_file}",
            shell=True,
        )
        aligned_file_handle = open(aligned_file)
        aligned_file_records = list(SeqIO.parse(aligned_file_handle, "fasta"))
        analyse = Sequence(aligned_file_records[0].seq, aligned_file_records[1].seq)
        (
            match_nucleotides,
            diff,
            missmatched,
            missligned,
            Ns,
        ) = analyse.get_match_nucleotides()
        results[sample] = get_output(
            analyse, match_nucleotides, diff, missmatched, missligned, Ns
        )

final_file = []
for sample in results:
    sample_list_values = [sample]
    for field in results[sample]:
        if type(results[sample][field]) != type({}):
            sample_list_values.append(results[sample][field])
        else:
            for subfield in results[sample][field]:
                sample_list_values.append(results[sample][field][subfield])
    final_file.append(sample_list_values)

fields = [
    "sample_name",
    "match",
    "snakemake length",
    "insaflu length",
    "missmatched",
    "missligned",
    "Ns",
    "total_differences",
    "total_differences %",
    "A>T",
    "A>C",
    "A>G",
    "A>-",
    "A>N",
    "T>A",
    "T>C",
    "T>G",
    "T>-",
    "T>N",
    "C>A",
    "C>T",
    "C>G",
    "C>-",
    "C>N",
    "G>A",
    "G>T",
    "G>C",
    "G>-",
    "G>N",
    "->A",
    "->T",
    "->C",
    "->G",
    "->N",
    "N>A",
    "N>T",
    "N>C",
    "N>G",
    "N>-",
]
with open(sys.argv[4], "w") as file:
    writer = csv.writer(file)
    writer.writerow(fields)
    writer.writerows(final_file)

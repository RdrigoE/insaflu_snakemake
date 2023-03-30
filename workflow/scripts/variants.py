import csv
import re
import sys

DESCRIPTORS = [
    "CHROM",
    "POS",
    "TYPE",
    "REF",
    "ALT",
    "FREQ",
    "COVERAGE",
    "EVIDENCE",
    "FTYPE",
    "STRAND",
    "NT_POS",
    "AA_POS",
    "EFFECT",
    "NT CHANGE",
    "AA CHANGE",
    "AA CHANGE ALT",
    "LOCUS_TAG",
    "GENE",
    "PRODUCT",
    "VARIANTS IN INCOMPLETE LOCUS",
]

#
# class SnakemakeEntry(TypedDict):
#     chrom: str
#     pos: int
#     ref: str
#     alt: str
#     qual: float
#     info: dict[str, str]
#     ann: dict[str, str]
#     names: list[str]
#     values: list[str]
#
#


def remove_smallcase(string):
    word = ""
    for letter in string:
        if not letter.islower():
            word += letter
    return word


def smaller(string: str) -> str:
    return string[: 2] + remove_smallcase(string[2:])


def convert_list_to_dict(data: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    snp_annotation = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID",
                      "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos/cDNA.length", "CDS_CDSLength", "AA.pos/AA.length", "Distance", "ERRORS/WARNINGS/INFO"]
    new_dict = {}

    ann_dict = {}
    for item in data:
        if item:
            k, v = item.split("=")
            if k.upper() == "ANN":
                ann_dict = {}
                new_v = v.split("|")
                for idx, v in enumerate(new_v):
                    ann_dict[snp_annotation[idx]] = v
            else:
                new_dict[k] = v
    return new_dict, ann_dict


# class DataEntry(TypedDict, total=False):
#     "CHROM"
#     "POS"
#     "TYPE"
#     "REF"
#     "ALT"
#     "FREQ"
#     "COVERAGE"
#     "EVIDENCE"
#     "FTYPE"
#     "STRAND"
#     "NT_POS"
#     "AA_POS"
#     "EFFECT"
#     "NT CHANGE"
#     "AA CHANGE"
#     "AA CHANGE ALT"
#     "LOCUS_TAG"
#     "GENE"
#     "PRODUCT"
#     "VARIANTS IN INCOMPLETE LOCUS"
#
#
def analyse_file(input_file):
    with open(input_file) as handler:
        lines = [line for line in handler.readlines() if not line[0] == '#']
    vcf_data = []
    for line in lines:
        l = line.split()
        info, ann = convert_list_to_dict(l[7].split(";"))
        entry = dict(
            chrom=l[0],
            pos=int(l[1]),
            ref=l[3],
            alt=l[4],
            qual=float(l[5]),
            info=info,
            ann=ann,
            names=l[8].split(":"),
            values=l[9].split(":")
        )
        print(entry["info"], input_file, entry["pos"])

        get_freq = round(float(entry["info"]["AO"].replace(
            ",", "."))/float(entry["info"]["DP"].replace(",", ".")) * 100, 1)
        default = [entry["chrom"],
                   entry["pos"],
                   entry["info"]["TYPE"],
                   entry["ref"],
                   entry["alt"],
                   get_freq,
                   entry["info"]["DP"],
                   entry["alt"] + ":" + entry["info"]["AO"],
                   entry["ref"] + ":" + entry["info"]["RO"],
                   ]

        if entry.get("ann"):
            ann = entry["ann"]

            default.extend(["CDS",
                            "+",
                            ann["CDS_CDSLength"],
                            ann["AA.pos/AA.length"],
                            ann["Annotation"],
                            ann["HGVS.c"],
                            ann["HGVS.p"],
                            smaller(ann["HGVS.p"]),
                            ann["Gene_Name"],
                            ann["Gene_Name"],
                            "yes"
                            ])
        else:
            default.extend(['', '', '', '', '', '', '', '', '', '', ''])
        vcf_data.append(default)
    return vcf_data


def word_in_list(word, check_list):
    for item in check_list:
        if item in word:
            return True
    return False


def comparison(value1, value2, signal):
    if signal == "bigger":
        return value1 >= value2
    elif signal == "smaller":
        return value1 < value2


def filter_variants(data: list[list[str]], signal: str, list_of_types: list[str], identifier: str) -> list[list[str]]:
    filtered_data = []
    for entry in data:
        entry_d = dict() # DataEntry()
        for idx, value in enumerate(entry):
            entry_d[DESCRIPTORS[idx]] = value
        if not comparison(entry_d["FREQ"], 0.51, signal):
            continue
        if not word_in_list(entry_d["TYPE"], list_of_types):
            continue
        entry.insert(0, identifier)
        filtered_data.append(entry)

    return filtered_data


def validated_variants(
    files_path, output_file, signal, re_expression, list_of_types
):
    """
    The validated_variants function takes a list of vcf files and outputs a csv file with the following columns:
        ID, CHROM, POS, TYPE, REF, ALT, FREQ (allele frequency), COVERAGE (total depth at that position), EVIDENCE (number of reads supporting variant)
        FTYPE(functional type), STRAND(strand bias for reads supporting variant), NT_POS(nucleotide position in codon)
        AA_POS(amino acid position in protein sequence) EFFECT(effect on gene function as reported by snpeff software package).


    :param files_path: Pass the path to the folder containing all of your vcf files
    :param output_file: Specify the name of the file that will be created
    :return: A csv file with the following columns:
    :doc-author: Trelent
    """
    vcf_final = [
        [
            "ID",
            "CHROM",
            "POS",
            "TYPE",
            "REF",
            "ALT",
            "FREQ",
            "COVERAGE",
            "EVIDENCE",
            "FTYPE",
            "STRAND",
            "NT_POS",
            "AA_POS",
            "EFFECT",
            "NT CHANGE",
            "AA CHANGE",
            "AA CHANGE ALT",
            "LOCUS_TAG",
            "GENE",
            "PRODUCT",
            "VARIANTS IN INCOMPLETE LOCUS",
        ]
    ]
    for file in files_path:
        identifier = re.findall(re_expression, file)[0]
        clean_file = analyse_file(file)
        filtered_file = filter_variants(
            clean_file, signal, list_of_types, identifier
        )
        vcf_final.extend(filtered_file)
    print(vcf_final)
    with open(output_file, mode="w") as f:
        f_writer = csv.writer(f, delimiter=",")
        f_writer.writerows(vcf_final)
        f.close()


if __name__ == "__main__":
    files_path = sys.argv[1].split(" ")
    output_file = sys.argv[2]
    type_of_file = sys.argv[3]
    if type_of_file == "validated_variants":
        signal = "bigger"
        re_expression = "(?<=sample__)(.*?)(?=_snpeff.vcf)"
        list_of_words = ["snp", "mnp", "complex", "ins", "del"]
    elif type_of_file == "minor_iSNVs_inc_indels":
        signal = "smaller"
        re_expression = "(?<=/freebayes/)(.*?)(?=_snpeff.vcf)"
        list_of_words = ["snp", "del", "ins"]
    elif type_of_file == "minor_iSNVs":
        signal = "smaller"
        re_expression = "(?<=/freebayes/)(.*?)(?=_snpeff.vcf)"
        list_of_words = ["snp"]
    else:
        raise TypeError  # to change
    validated_variants(
        files_path,
        output_file,
        signal,
        re_expression,
        list_of_words,
    )

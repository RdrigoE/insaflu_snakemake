"""Tbd"""
import os
import sys
from Bio import SeqIO
from Bio.Seq import MutableSeq
from pysam import VariantFile, VariantRecord


def get_type_variation(ref, alt):
    """return type of variation based on change
    possible in variation type
                    snp	Single Nucleotide Polymorphism	A => T
                    mnp	Multiple Nuclotide Polymorphism	GC => AT
                    ins	Insertion	ATT => AGTT
                    del	Deletion	ACGG => ACG
                    complex	Combination of snp/mnp
    """
    if len(ref) == len(alt) and len(alt) == 1:
        return "snp"
    if len(ref) > len(alt):
        return "del"
    if len(ref) < len(alt):
        return "ins"
    if len(ref) == len(alt) and len(alt) > 1:
        return "mnp"
    return "complex"


def read_text_file(file_name):
    """
    read text file and put the result in an vector
    """
    vect_out = []
    with open(file_name) as handle:
        for line in handle:
            sz_temp = line.strip()
            if len(sz_temp) == 0:
                continue
            vect_out.append(sz_temp)
    return vect_out


def get_coverage_by_pos(
    file_coverage, chr_name, position_start, position_end, temp_file
):
    # print(" i was called get_Coveraaaage")
    if not os.path.exists(file_coverage):

        print(" failed")
        return -1

    cmd = "{} {} {}:{}-{} > {}".format(
        "tabix", file_coverage, chr_name, position_start, position_end, temp_file
    )
    os.system(cmd)
    ### get number
    vect_lines = read_text_file(temp_file)
    # print(len(vect_lines) == 1 and len(vect_lines[0].split("\t")) == 3)
    if (
        len(vect_lines) == 1
        and len(vect_lines[0].split("\t")) == 3
        and isinstance(int(vect_lines[0].split("\t")[2]), int)
    ):
        return int(vect_lines[0].split("\t")[2])

    return -1


def add_freq_ao_ad_and_type_to_vcf(
    vcf_file,
    vcf_file_out,
    vcf_file_out_removed_by_filter,
    freq_vcf_limit,
    coverage_limit,
    file_coverage,
):
    # Set constants
    FREQ: str = "FREQ"
    AO: str = "AO"
    RO: str = "RO"
    AF: str = "AF"
    TYPE: str = "TYPE"
    DP_COMPOSED: str = "DP_COMPOSED"

    # freq="Ratio of AO/(DPSP-AR)"

    vcf_hanlder = VariantFile(vcf_file, "r")

    vcf_hanlder_write = VariantFile(vcf_file_out, "w")

    vcf_handler_write_rmv_by_filter = VariantFile(vcf_file_out_removed_by_filter, "w")

    variant: VariantRecord
    for variant in vcf_hanlder:
        # print("Variant in pos", variant.__str__().split("\t")[1])
        ### DP must be replaced by DPSP. DPSP is the sum of all reads Span and Ambiguous
        if (
            "SR" in variant.info and "DPSP" in variant.info and "AR" in variant.info
        ):  ## SR=0,0,15,6
            ### don't process this VCF because has a low coverage
            total_deep = int(variant.info["DPSP"]) - sum(
                [int(_) for _ in variant.info["AR"]]
            )
            total_deep_samtools = get_coverage_by_pos(
                file_coverage, variant.chrom, variant.pos, variant.pos, sys.argv[6]
            )
            # print(
            #     f"Total Deepth: {total_deep} \t Total_depth_samtools {total_deep_samtools}"
            # )
            if (
                coverage_limit > 0
                and total_deep_samtools >= 0
                and total_deep_samtools < coverage_limit
            ):
                continue
            if ((len(variant.info["SR"]) // 2) - 1) != len(variant.alts):  # type: ignore
                # vcf_hanlder_write.write(variant)
                continue  ### different numbers of Alleles and References

            #### extra info
            vect_out_ao = []  ### AO
            out_ro = -1  ### RO
            vect_out_af: list[int] = []  ### AF
            vect_out_freq: list[float] = []  ### FREQ
            vect_out_freq_filtered = []  ### FREQ
            vect_out_type = []  ### TYPE
            for value_ in range(0, len(variant.info["SR"]), 2):
                # print(value_)
                if value_ > 0:
                    allele_count = int(variant.info["SR"][value_]) + int(
                        variant.info["SR"][value_ + 1]
                    )
                    # print(allele_count)

                    if total_deep > 0:
                        ### incongruences in Medaka,
                        ### these values are collected in different stages of the Medaka workflow, (email from support@nanoporetech.com at 23 Dec 2020)
                        # print(
                        #     f"Total deep: {total_deep} \t Allele Count: {allele_count}"
                        # )
                        if total_deep <= allele_count:
                            vect_out_freq.append(100)
                        else:
                            freq_value = allele_count / float(total_deep)
                            if freq_value >= freq_vcf_limit / 100:
                                vect_out_freq.append(
                                    float("{:.1f}".format(freq_value * 100))
                                )
                            elif vcf_file_out_removed_by_filter:
                                vect_out_freq_filtered.append(
                                    float("{:.1f}".format(freq_value * 100))
                                )
                        # print(
                        #     variant.pos,
                        #     variant.ref,
                        #     str(variant.alts),
                        #     variant.info["DP"],
                        #     vect_out_ao,
                        #     vect_out_freq,
                        # )

                    vect_out_ao.append(allele_count)
                    vect_out_type.append(
                        get_type_variation(variant.ref, variant.alts[(value_ - 2) >> 1])
                    )
                    # These are bitwise shift operators.
                    # Quoting from the docs:
                    # x << y
                    # Returns x with the bits shifted to the left by y places (and new bits on the right-hand-side are zeros). This is the same as multiplying x by 2**y.
                    # x >> y
                    # Returns x with the bits shifted to the right by y places. This is the same as dividing x by 2**y
                vect_out_af.append(
                    int(variant.info["SR"][value_])
                    + int(variant.info["SR"][value_ + 1])
                )
                if out_ro == -1:
                    out_ro = int(variant.info["SR"][value_]) + int(
                        variant.info["SR"][value_ + 1]
                    )

            ### has some variant to save
            if len(vect_out_freq) > 0:
                ### Only save the ones with FREQ
                # vcf_hanlder_write.write(variant)
                print(
                    variant.pos,
                    tuple([out_ro]),
                    tuple(vect_out_ao),
                    tuple(vect_out_af),
                    tuple(vect_out_type),
                    tuple([f"{total_deep}/{total_deep_samtools}"]),
                    tuple(vect_out_freq),
                )

            ### save the filtered
            if len(vect_out_freq_filtered) > 0:
                print(
                    variant.pos,
                    tuple([out_ro]),
                    tuple(vect_out_ao),
                    tuple(vect_out_af),
                    tuple(vect_out_type),
                    tuple([f"{total_deep}/{total_deep_samtools}"]),
                    tuple(vect_out_freq),
                )
                # if vect_out_freq_filtered[0] > 50 and vect_out_freq_filtered[0] < 80:
                ### Only save the ones with FREQ
                vcf_handler_write_rmv_by_filter.write(variant)
    # print(count)
    vcf_hanlder_write.close()
    vcf_hanlder.close()
    vcf_handler_write_rmv_by_filter.close()
    return vcf_file_out


if __name__ == "__main__":
    freq_vcf_limit = 80
    coverage_limit = 1
    element_name_old = ""
    vect_sites = []
    vect_ranges = []
    final_mask_consensus = []

    vcf_file = sys.argv[2]
    vcf_file_out = sys.argv[4]
    vcf_file_out_removed_by_filter = sys.argv[5]
    file_coverage = sys.argv[3]

    add_freq_ao_ad_and_type_to_vcf(
        vcf_file,
        vcf_file_out,
        vcf_file_out_removed_by_filter,
        freq_vcf_limit,
        coverage_limit,
        file_coverage,
    )

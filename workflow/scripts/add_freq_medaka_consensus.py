"""Tbd"""
import os
import sys
from Bio import SeqIO
from Bio.Seq import MutableSeq
import pysam


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

    if not os.path.exists(file_coverage):

        return -1
    cmd = "{} {} {}:{}-{} > {}".format(
        "tabix", file_coverage, chr_name, position_start, position_end, temp_file
    )
    os.system(cmd)
    ### get number
    vect_lines = read_text_file(temp_file)
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
    temp_file,
):
    """add FREQ, AO, AF and TYPE to VCF, FREQ=AO/DP
    This case is used in MEDAKA only
    :param vcf_file_out_removed_by_filter -> can be None, keep the variants that are filter by  freq_vcf_limit
    vcffile must be gzip and tbi included
    :param coverage_limit -> filter by this coverage (this is necessary because medaka doesn't have)
    :param cut off for VCF freq
    returns: vcf file with freq, AO and AF
    """
    FREQ = "FREQ"
    AO = "AO"
    RO = "RO"
    AF = "AF"
    TYPE = "TYPE"
    DP_COMPOSED = "DP_COMPOSED"  ### this is used to get

    # read the input file
    vcf_hanlder = pysam.VariantFile(vcf_file, "r")
    if (
        FREQ in vcf_hanlder.header.info
        and AO in vcf_hanlder.header.info
        and AF in vcf_hanlder.header.info
    ):
        vcf_hanlder.close()
        return

    vcf_hanlder_write = pysam.VariantFile(vcf_file_out, "w")
    if not vcf_file_out_removed_by_filter is None:
        vcf_hanlder_write_removed_by_filter = pysam.VariantFile(
            vcf_file_out_removed_by_filter, "w"
        )
    if not FREQ in vcf_hanlder.header.info:
        vcf_hanlder.header.info.add(
            FREQ, number="A", type="Float", description="Ratio of AO/(DPSP-AR)"
        )
    if not AO in vcf_hanlder.header.info:
        vcf_hanlder.header.info.add(
            AO,
            number="A",
            type="Integer",
            description="Alternate allele observation count, SR (alt1 fwd + alt1 rev, etc.)",
        )
    if not RO in vcf_hanlder.header.info:
        vcf_hanlder.header.info.add(
            RO,
            number="1",
            type="Integer",
            description="Reference allele observation count, SR (ref fwd + ref rev)",
        )
    if not AF in vcf_hanlder.header.info:
        vcf_hanlder.header.info.add(
            AF,
            number="R",
            type="Integer",
            description="Number of observation for each allele, SR (ref fwd + ref rev, alt1 fwd + alt1 rev, etc.)",
        )
    if not TYPE in vcf_hanlder.header.info:
        vcf_hanlder.header.info.add(
            TYPE,
            number="A",
            type="String",
            description="The type of allele, either snp, mnp, ins, del, or complex",
        )
    if not DP_COMPOSED in vcf_hanlder.header.info:
        vcf_hanlder.header.info.add(
            DP_COMPOSED,
            number="1",
            type="String",
            description="Coverage at position (DPSP-AR)/(samtools -aa). First is collected by Medaka, Second is collected by samtools.",
        )

    ## write the header
    for variant_header_records in vcf_hanlder.header.records:
        vcf_hanlder_write.header.add_record(variant_header_records)
        if not vcf_file_out_removed_by_filter is None:
            vcf_hanlder_write_removed_by_filter.header.add_record(
                variant_header_records
            )

    for variant_sample in vcf_hanlder.header.samples:
        vcf_hanlder_write.header.add_sample(variant_sample)
        if not vcf_file_out_removed_by_filter is None:
            vcf_hanlder_write_removed_by_filter.header.add_sample(variant_sample)

    for variant in vcf_hanlder:
        ### DP must be replaced by DPSP. DPSP is the sum of all reads Span and Ambiguous
        if (
            "SR" in variant.info and "DPSP" in variant.info and "AR" in variant.info
        ):  ## SR=0,0,15,6
            ### don't process this VCF because has a low coverage
            total_deep = int(variant.info["DPSP"]) - sum(
                [int(_) for _ in variant.info["AR"]]
            )
            total_deep_samtools = get_coverage_by_pos(
                file_coverage, variant.chrom, variant.pos, variant.pos, temp_file
            )
            if (
                coverage_limit > 0
                and total_deep_samtools >= 0
                and total_deep_samtools < coverage_limit
            ):
                continue
            if ((len(variant.info["SR"]) // 2) - 1) != len(variant.alts):
                # vcf_hanlder_write.write(variant)
                continue  ### different numbers of Alleles and References

            #### extra info
            vect_out_ao = []  ### AO
            out_ro = -1  ### RO
            vect_out_af = []  ### AF
            vect_out_freq = []  ### FREQ
            vect_out_freq_filtered = []  ### FREQ
            vect_out_type = []  ### TYPE

            for value_ in range(0, len(variant.info["SR"]), 2):
                if value_ > 0:
                    allele_count = int(variant.info["SR"][value_]) + int(
                        variant.info["SR"][value_ + 1]
                    )

                    if total_deep > 0:
                        ### incongruences in Medaka,
                        ### these values are collected in different stages of the Medaka workflow, (email from support@nanoporetech.com at 23 Dec 2020)
                        if total_deep <= allele_count:
                            vect_out_freq.append(100)
                        else:
                            freq_value = allele_count / float(total_deep)
                            if freq_value >= freq_vcf_limit:
                                vect_out_freq.append(
                                    float("{:.1f}".format(freq_value * 100))
                                )
                            elif not vcf_file_out_removed_by_filter is None:

                                vect_out_freq_filtered.append(
                                    float("{:.1f}".format(freq_value * 100))
                                )
                        # print(variant.pos, variant.ref, str(variant.alts), variant.info['DP'], vect_out_ao[-1], vect_out_freq[-1])

                    vect_out_ao.append(allele_count)
                    vect_out_type.append(
                        get_type_variation(variant.ref, variant.alts[(value_ - 2) >> 1])
                    )
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
                if out_ro > -1:
                    variant.info[RO] = tuple([out_ro])
                variant.info[AO] = tuple(vect_out_ao)
                variant.info[AF] = tuple(vect_out_af)
                variant.info[TYPE] = tuple(vect_out_type)
                variant.info[DP_COMPOSED] = tuple(
                    ["{}/{}".format(total_deep, total_deep_samtools)]
                )
                variant.info[FREQ] = tuple(vect_out_freq)

                ### Only save the ones with FREQ
                vcf_hanlder_write.write(variant)

            ### save the filtered
            if len(vect_out_freq_filtered) > 0:
                if out_ro > -1:
                    variant.info[RO] = tuple([out_ro])
                variant.info[AO] = tuple(vect_out_ao)
                variant.info[AF] = tuple(vect_out_af)
                variant.info[TYPE] = tuple(vect_out_type)
                variant.info[DP_COMPOSED] = tuple(
                    ["{}/{}".format(total_deep, total_deep_samtools)]
                )
                variant.info[FREQ] = tuple(vect_out_freq_filtered)

                ### Only save the ones with FREQ
                vcf_hanlder_write_removed_by_filter.write(variant)

    vcf_hanlder_write.close()
    vcf_hanlder.close()
    if not vcf_file_out_removed_by_filter is None:
        vcf_hanlder_write_removed_by_filter.close()
    return vcf_file_out


def compute_masking_sites(
    sequence, ranges=None, single_positions=None, from_beggining=None, from_end=None
) -> list[int]:
    length = len(sequence)
    masking_sites = []
    if ranges is not None:
        for pos in ranges:
            # print(ranges)
            pos[0] = int(pos[0])
            pos[1] = int(pos[1])
            if pos[1] < length and pos[1] > pos[0]:
                for position in range(pos[0] - 1, pos[1] - 1):
                    masking_sites.append(position)
    if single_positions is not None:
        for pos in single_positions:
            pos = int(pos)
            if pos < length:
                masking_sites.append(pos - 1)
    if from_beggining is not None:
        for pos in range(int(from_beggining)):
            masking_sites.append(pos)
    if from_end is not None:
        for pos in range(length - int(from_end), length):
            masking_sites.append(pos)
    return masking_sites


def main():
    freq_vcf_limit: float = float(sys.argv[10])
    coverage_limit: float = float(sys.argv[9])
    element_name_old = ""
    vect_sites = []
    vect_ranges = []
    final_mask_consensus = []

    vcf_file = sys.argv[2]
    vcf_file_out = sys.argv[4]
    vcf_file_out_removed_by_filter = sys.argv[5]
    file_coverage = sys.argv[3]
    temp_file = sys.argv[8]
    add_freq_ao_ad_and_type_to_vcf(
        vcf_file,
        vcf_file_out,
        vcf_file_out_removed_by_filter,
        freq_vcf_limit,
        coverage_limit,
        file_coverage,
        temp_file,
    )

    final_vcf_with_removed_variants = sys.argv[5]

    vcf_hanlder = pysam.VariantFile(final_vcf_with_removed_variants, "r")

    vect_sites = []
    vect_ranges = []
    for variant in vcf_hanlder:
        if element_name_old != variant.chrom:
            element_name_old = variant.chrom
        ### MEDAKA output must have "TYPE" in info
        if variant.info["TYPE"][0] == "snp":
            vect_sites.append(str(variant.pos))
        elif variant.info["TYPE"][0] == "ins":
            vect_sites.append(str(variant.pos))
        elif variant.info["TYPE"][0] == "del":
            vect_ranges.append(
                "{}-{}".format(
                    variant.pos + len(variant.alts[0]),
                    variant.pos - len(variant.alts[0]) + len(variant.ref),
                )
            )
        else:
            vect_ranges.append(
                "{}-{}".format(variant.pos, variant.pos + len(variant.ref) - 1)
            )

    sites = ""
    for idx, i in enumerate(vect_sites):
        if idx == len(vect_sites) - 1:
            sites += i
        else:
            sites += i + ","

    regions = ""
    for idx, i in enumerate(vect_ranges):
        if idx == len(vect_ranges) - 1:
            regions += i
        else:
            regions += i + ","

    # print(sites)
    # print(regions)

    consensus = sys.argv[1]
    output = sys.argv[7]
    final_mask_consensus = []

    ranges, single_positions, from_beggining, from_end = (
        [x.split("-") for x in regions.split(",")] if regions != "" else None,
        sites.split(",") if sites != "" else None,
        None,
        None,
    )

    for record in SeqIO.parse(consensus, "fasta"):
        sequence = MutableSeq(record.seq)
        masking_sites = compute_masking_sites(
            sequence, ranges, single_positions, from_beggining, from_end
        )
        masking_sites = sorted(masking_sites)
        # print(masking_sites)
        ### Taken from insaflu
        ref_pos = 0
        ref_insertions = 0
        gap = 0
        for _ in range(len(sequence)):
            if sequence[_] == "-":
                # print("Hello")
                gap += 1
                # print(gap)

            # ref_insertions += 1
            # continue
            # print(_, ref_pos)
            if ref_pos in masking_sites:
                # print(
                #     f"placing N in position {ref_pos + ref_insertions}, it was a {sequence[ref_pos + ref_insertions]}"
                # )
                sequence[ref_pos + ref_insertions] = "N"
            ref_pos += 1
            if (ref_pos + ref_insertions) >= len(record.seq):
                break
        ### End of insaflu code
        record.seq = sequence
        final_mask_consensus.append(record)

    SeqIO.write(final_mask_consensus, output, "fasta")


if __name__ == "__main__":
    main()

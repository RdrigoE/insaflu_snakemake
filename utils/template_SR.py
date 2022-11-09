def add_freq_ao_ad_and_type_to_vcf(
    self,
    vcf_file,
    file_coverage,
    vcf_file_out,
    vcf_file_out_removed_by_filter,
    coverage_limit,
    freq_vcf_limit,
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
            total_deep_samtools = self.get_coverage_by_pos(
                file_coverage, variant.chrom, variant.pos, variant.pos
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
                        self._get_type_variation(
                            variant.ref, variant.alts[(value_ - 2) >> 1]
                        )
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

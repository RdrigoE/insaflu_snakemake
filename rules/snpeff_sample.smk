configfile: "config/parameters.yaml"


replace = (
    f"sed -i 's/{SEGMENTS[0]}/{get_id_version(REFERENCE_GB)}/g' "
    if len(SEGMENTS) == 1
    else "true "
)


rule snpeff_sample:
    input:
        i1="projects/{project}/sample_{sample}/snps.vcf",
        i2="projects/{project}/main_result/snp_ready.txt",
    output:
        "projects/{project}/sample_{sample}/{sample}_snpeff.vcf",
    conda:
        "../envs/snpeff.yaml"
    threads: config["snpeff_threads"]
    params:
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats -c config/snpeff.config",
    shell:
        "{replace} {input.i1} && "
        "snpEff {params} -v {REFERENCE_NAME} {input.i1}  > {output}"
        #-t {threads}

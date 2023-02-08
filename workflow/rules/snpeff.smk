configfile: "../config/threads.yaml"


replace = (
    f"sed -i 's/{SEGMENTS[0]}/{get_identification_version_string(SEGMENTS,REFERENCE_GB)}/g' "
    if len(SEGMENTS) == 1
    else "true "
)


rule prepare_snpeff:
    input:
        ref_gb=REFERENCE_GB,
        ref_fa=REFERENCE_FASTA,
    output:
        "projects/{project}/main_result/snp_ready.txt",
    conda:
        "../envs/snpeff.yaml"
    log:
        "logs/projects/{project}/main_result/snpeff_prepare.log",
    benchmark:
        "benchmark/projects/{project}/main_result/snpeff_prepare.tsv",
    shell:
        "python {scripts_directory}create_snpeff_text.py $CONDA_PREFIX {input.ref_gb} {input.ref_fa} {REFERENCE_NAME} {output} "


rule snpeff:
    input:
        snp_file="projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
        dependency_file="projects/{project}/main_result/snp_ready.txt",
    output:
        "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
    conda:
        "../envs/snpeff.yaml"
    params:
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats -c ../workflow/db/snpeff.config",
    log:
        "logs/projects/{project}/{sample}/snpeff.log",
    benchmark:
        "benchmark/projects/{project}/{sample}/snpeff.tsv",
    shell:
        "{replace} {input.snp_file} &&"
        "snpEff {params} -v {REFERENCE_NAME} {input.snp_file} > {output}"


rule snpeff_sample:
    input:
        i1="projects/{project}/sample_{sample}/snps.vcf",
        i2="projects/{project}/main_result/snp_ready.txt",
    output:
        "projects/{project}/sample_{sample}/sample__{sample}_snpeff.vcf",
    conda:
        "../envs/snpeff.yaml"
    log:
        "logs/projects/{project}/{sample}/{sample}_snpeff.log",
    benchmark:
        "benchmark/projects/{project}/{sample}/{sample}_snpeff.tsv",
    params:
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats -c ../workflow/db/snpeff.config",
    shell:
        "{replace} {input.i1} && "
        "snpEff {params} -v {REFERENCE_NAME} {input.i1}  > {output}"

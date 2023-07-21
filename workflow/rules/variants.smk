localrules: variant_validated,snpeff_concat,snpeff_concat_indels,proportions_iSNVs_graph,
rule variant_validated:
    input:
        expand(
            "projects/{project}/sample_{sample}/sample__{sample}_snpeff.vcf",
            sample=config_user["samples"],
            project=config_user["project"],
        ),
    conda:
        "../envs/base.yaml"
    output:
        "projects/{project}/main_result/validated_variants.csv",
    resources:
        mem_mb=memory["variant_validated"],
    log:
        "logs/projects/{project}/main_result/validated_variants.log",
    benchmark:
        "benchmark/projects/{project}/main_result/validated_variants.tsv"
    shell:
        "python {scripts_directory}variants.py '{input}' '{output}' validated_variants"


rule snpeff_concat:
    input:
        expand(
            "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
            sample=illumina_samples,
            project=config_user["project"],
        ),
    conda:
        "../envs/base.yaml"
    output:
        "projects/{project}/main_result/validated_minor_iSNVs.csv",
    resources:
        mem_mb=memory["snpeff_concat"],
    log:
        "logs/projects/{project}/main_result/validated_minor_iSNVs.log",
    benchmark:
        "benchmark/projects/{project}/main_result/validated_minor_iSNVs.tsv"
    shell:
        "python {scripts_directory}variants.py '{input}' {output} minor_iSNVs"


rule snpeff_concat_indels:
    input:
        expand(
            "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
            sample=illumina_samples,
            project=config_user["project"],
        ),
    conda:
        "../envs/base.yaml"
    output:
        "projects/{project}/main_result/validated_minor_iSNVs_inc_indels.csv",
    resources:
        mem_mb=memory["snpeff_concat_indels"],
    log:
        "logs/projects/{project}/main_result/validated_minor_iSNVs_inc_indels.log",
    benchmark:
        "benchmark/projects/{project}/main_result/validated_minor_iSNVs_inc_indels.tsv"
    shell:
        "python {scripts_directory}variants.py '{input}' {output} minor_iSNVs_inc_indels"


rule proportions_iSNVs_graph:
    input:
        expand(
            "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
            sample=config_user["samples"],
            project=config_user["project"],
        ),
    output:
        out_file="projects/{project}/main_result/proportions_iSNVs_graph.csv",
    conda:
        "../envs/base.yaml"
    resources:
        mem_mb=memory["proportions_iSNVs_graph"],
    log:
        "logs/projects/{project}/main_result/proportions_iSNVs_graph.log",
    benchmark:
        "benchmark/projects/{project}/main_result/proportions_iSNVs_graph.tsv"
    shell:
        "python {scripts_directory}proportions_iSNVs_graph.py '{input}' {output.out_file}"

rule variant_validated:
    input:
        expand(
            "projects/{project}/sample_{sample}/{sample}_snpeff.vcf",
            sample=config_user["samples"],
            project=config_user["project"],
        ),
    conda:
        "../envs/base.yaml"
    output:
        "projects/{project}/main_result/validated_variants.csv",
    shell:
        "python {scripts_directory}variants.py '{input}' '{output}' validated_variants"


rule snpeff_concat:
    input:
        expand(
            "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
            sample=config_user["samples"],
            project=config_user["project"],
        ),
    conda:
        "../envs/base.yaml"
    output:
        "projects/{project}/main_result/validated_minor_iSNVs.csv",
    shell:
        "python {scripts_directory}variants.py '{input}' {output} minor_iSNVs"


rule snpeff_concat_indels:
    input:
        expand(
            "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
            sample=config_user["samples"],
            project=config_user["project"],
        ),
    conda:
        "../envs/base.yaml"
    output:
        "projects/{project}/main_result/validated_minor_iSNVs_inc_indels.csv",
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
    log:
        "projects/{project}/main_result/proportions_iSNVs_graph.log",
    shell:
        "python {scripts_directory}proportions_iSNVs_graph.py '{input}' {output.out_file}"

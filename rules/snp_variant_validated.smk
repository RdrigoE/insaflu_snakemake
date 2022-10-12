rule variant_validated:
    input:
        expand("projects/{project}/sample_{sample}/{sample}_snpeff.vcf",sample=config_user['samples'], project=config_user['project']),
    output:
        "projects/{project}/main_result/validated_variants.csv"
    shell:
        "python utils/validated_variants.py '{input}' '{output}'"
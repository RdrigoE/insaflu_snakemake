rule snpeff_sample:
    input:
        "projects/{project}/sample_{sample}/snippy/snps.vcf"
    output:
        "projects/{project}/main_result/{sample}_snpeff.vcf"
    conda:
        "../envs/snpeff.yaml"
    params:
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats "#-c config/snpeff.config"
    shell:
        "snpEff {params} -v MN908947.3 {input}  > {output}"
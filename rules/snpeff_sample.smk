rule snpeff_sample:
    input:
        i1 = "projects/{project}/sample_{sample}/snippy/snps.vcf",
        i2 = "projects/{project}/main_result/snpeff/ready.txt"
    output:
        "projects/{project}/main_result/snpeff_samples/{sample}_snpeff.vcf"
    conda:
        "../envs/snpeff.yaml"
    params:
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats -c config/snpeff.config"
    shell:
        #"snpEff {params} -v MN908947.3 {input.1}  > {output}"
        "snpEff {params} -v {REFERENCE_NAME} {input.i1}  > {output}"

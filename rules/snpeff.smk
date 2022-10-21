configfile: "config/parameters.yaml"

replace = f"sed -i 's/{SEGMENTS[0]}/{get_id_version(REFERENCE_GB)}/g' " if len(SEGMENTS) == 1 else 'true '

rule prepare_snpeff:
    input:
        ref_gb = REFERENCE_GB,
        ref_fa = REFERENCE,
    output:
        temp("projects/{project}/main_result/snp_ready.txt"),
    conda:
        "../envs/snpeff.yaml"
    shell:
        "python utils/create_snpeff_text.py $CONDA_PREFIX {input.ref_gb} {input.ref_fa} {REFERENCE_NAME} {output} "
        

rule snpeff:
    input:
        snp_file = "projects/{project}/sample_{sample}/freebayes/{sample}_var.vcf",
        dependency_file = "projects/{project}/main_result/snp_ready.txt"
    output:
        "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
    conda:
        "../envs/snpeff.yaml"
    threads: 
        config['snpeff_threads']
    params:
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats -c config/snpeff.config"
    shell:
        #not ready for multithreading >< -t {threads}
        "echo '{replace}' && "
        "{replace} {input.snp_file} &&"
        "snpEff {params} -v {REFERENCE_NAME} {input.snp_file} > {output}"


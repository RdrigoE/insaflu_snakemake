checkpoint abricate_pangolin:
    output:
        csv=temp(expand("projects/{project}/ref.csv", project=config_user["project"])),
        yaml=temp(expand("projects/{project}/ref.yaml", project=config_user["project"])),
    conda:
        "../envs/abricate.yaml"
    params:
        "--db insaflu --minid 70 --mincov 30",
    shell:
        "abricate {params} {REFERENCE_FASTA} > {output.csv}"
        " && python3 {scripts_directory}get_abricate_info.py {output.csv} {output.yaml}"


rule not_pangolin:
    input:
        yaml=expand("projects/{project}/ref.yaml", project=config_user["project"]),
    output:
        temp("projects/{project}/main_result/not_pangolin.csv"),
    shell:
        "echo 'This is not SARS-CoV-2' > {output}"


rule pangolin:
    input:
        yaml=expand("projects/{project}/ref.yaml", project=config_user["project"]),
        consensus="projects/{project}/main_result/AllConsensus_no_ref.fasta",
    output:
        "projects/{project}/main_result/lineage_report.csv",
    conda:
        "../envs/pangolin.yaml"
    params:
        "--analysis-mode fast",  #pangolearn
    shell:
        "pangolin {input.consensus} --outfile {output} -t 2"

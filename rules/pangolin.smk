rule pangolin:
    input:
        "projects/{project}/main_result/AllConsensus_no_ref.fasta"
    output:
        "projects/{project}/main_result/lineage_report.csv"
    conda:
        "../envs/pangolin.yaml"
    params:
        "--analysis-mode fast"
    shell:
        "pangolin {input} --outfile {output} -t 2"

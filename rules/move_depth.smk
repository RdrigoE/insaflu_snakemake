rule move_depth:
    input:
        "projects/{project}/sample_{sample}/snippy/snps.depth.gz"
    output:
        "projects/{project}/main_result/depth/{sample}__SARS_COV_2.depth"
    shell:
        "gzip -d -k -c {input} > {output}"
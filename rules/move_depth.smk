rule move_depth:
    input:
        "projects/{project}/sample_{sample}/snps.depth.gz",
    output:
        temp("projects/{project}/main_result/depth/{sample}__{ref}.depth"),
    shell:
        "gzip -d -k -c {input} > {output} && python utils/split_depth_file.py {output} {REFERENCE_GB}"

rule move_depth:
    input:
        "projects/{project}/sample_{sample}/snippy/snps.depth.gz"
    output:
        "projects/{project}/main_result/depth/{sample}.depth"
    shell:
        "gzip -d -k -c {input} > {output} | python utils/split_depth_file.py {output}"
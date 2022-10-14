rule move_depth:
    input:
        i = "projects/{project}/sample_{sample}/snps.depth.gz",
        ref_gb = REFERENCE_GB
    output:
        o = temp("projects/{project}/main_result/depth/{sample}__{ref}.depth"),
    shell:
        "gzip -d -k -c {input.i} > {output.o} && python utils/split_depth_file.py {output.o} {input.ref_gb}"

with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

rule move_depth:
    input:
        "projects/{project}/sample_{sample}/snippy/snps.depth.gz"
    output:
        o = "projects/{project}/main_result/depth/{sample}__{ref}.depth",
    shell:
        "gzip -d -k -c {input} > {output.o} | python utils/split_depth_file.py {output.o}"

rule move_depth_keep:
    input:
        "projects/{project}/sample_{sample}/snippy/snps.depth.gz"
    output:
        o = "projects/{project}/main_result/depth/{sample}.depth"
    shell:
        "gzip -d -k -c {input} > {output.o}"
#not being used at the time of realease 1.0
rule rabbitqc:
    input:
        i1 = "user_data/{sample}.fastq.gz",
        i2 = "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
    output:
        dir=directory("samples/{sample}/rabbitqc/"),
        a = "samples/{sample}/rabbitqc/rabbit.html"
    conda:
        "../envs/rabbitqc.yaml"
    threads:
        6
    params:
        "-w 3 -D"
    shell:
        "rabbit_qc -w {threads} -i {input.i1} {input.i2} -h {output.dir} {params}"

configfile: "../config/threads.yaml"


rule spades_se:
    input:
        "samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz",
    output:
        "samples/{sample}/spades_se/contigs.fasta",
        dir=directory("samples/{sample}/spades_se/"),
    conda:
        "../envs/spades.yaml"
    threads: config["spades_threads"]
    params:
        "--only-assembler",
    log:
        "spades/{sample}/spades.log",
    shell:
        "spades.py -t {threads} {params} -s {input} -o {output.dir}"


rule spades_pe:
    input:
        read_1="samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        read_2="samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",
    output:
        o="samples/{sample}/spades_pe/contigs.fasta",
        dir=directory("samples/{sample}/spades_pe/"),
    threads: config["spades_threads"]
    conda:
        "../envs/spades.yaml"
    params:
        "--only-assembler",
    log:
        "spades/{sample}/spades.log",
    shell:
        "spades.py -t {threads}  {params} -1 {input.read_1} -2 {input.read_2} -o {output.dir}"

rule spades_se:
    input:
        "samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz"
    output:
        "samples/{sample}/spades/contigs.fasta",
        dir = directory('samples/{sample}/spades/')
    conda:
        "../envs/spades.yaml"
    params:
        "--only-assembler"
    shell:
        'spades.py {params} -s {input} -o {output.dir}' #isolate is just to keep this working

rule spades_pe:
    input:
        i1 = "samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        i2 = "samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",
    output:
        o = "samples/{sample}/spades/contigs.fasta",
        dir = directory('samples/{sample}/spades/')
    conda:
        "../envs/spades.yaml"
    params:
        "--only-assembler --isolate"
    shell:
        'spades.py {params} -1 {input.i1} -2 {input.i2} -o {output.dir}' #isolate is just to keep this working
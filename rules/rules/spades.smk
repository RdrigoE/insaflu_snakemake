rule spades_se:
    input:
        "results/trimmed_reads_se/{sample}.trimmed.fastq.gz"
    output:
        "results/spades_se/{sample}/contigs.fasta",
        dir = directory('results/spades_se/{sample}')
    conda:
        "../envs/spades.yaml"
    shell:
        'spades.py --only-assembler -s {input} -o {output.dir}' #isolate is just to keep this working

rule spades_pe:
    input:
        i1 = "results/trimmed_reads_pe/{sample}_1.trimmed.fastq.gz",
        i2 = "results/trimmed_reads_pe/{sample}_2.trimmed.fastq.gz",
    output:
        o = "results/spades_pe/{sample}/contigs.fasta",
        dir = directory('results/spades_pe/{sample}')
    conda:
        "../envs/spades.yaml"
    shell:
        'spades.py --only-assembler -1 {input.i1} -2 {input.i2} -o {output.dir}' #isolate is just to keep this working
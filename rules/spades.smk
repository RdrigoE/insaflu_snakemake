rule spades_se:
    input:
        "results/trimmed_reads_se/{sample}.trimmed.fastq.gz"
    output:
        "results/spades/{sample}_spades/contigs.fasta",
        dir = directory('results/spades/{sample}_spades')
    conda:
        "../envs/spades.yaml"
    shell:
        'spades.py --only-assembler -s {input} -o {output.dir}' #isolate is just to keep this working
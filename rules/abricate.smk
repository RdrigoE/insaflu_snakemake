rule abricate:
    input:
        "results/spades/{sample}_spades/contigs.fasta"
    output:
        "results/abricate/{sample}/abricate_{sample}.csv"
    conda:
        "../envs/abricate.yaml"
    shell:
        "abricate --minid 70 --mincov 60 {input} > {output}"
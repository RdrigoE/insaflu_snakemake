rule abricate_se:
    input:
        "results/spades_se/{sample}/contigs.fasta"
    output:
        "results/abricate_se/{sample}/abricate_{sample}.csv"
    conda:
        "../envs/abricate.yaml"
    shell:
        "abricate --db insaflu --minid 70 --mincov 60 {input} > {output}"

rule abricate_pe:
    input:
        "results/spades_pe/{sample}/contigs.fasta"
    output:
        "results/abricate_pe/{sample}/abricate_{sample}.csv"
    conda:
        "../envs/abricate.yaml"
    shell:
        "abricate --db insaflu --minid 70 --mincov 60 {input} > {output}"


rule abricate:
    input:
        "samples/{sample}/spades/contigs.fasta"
    output:
        "samples/{sample}/abricate/abricate_{sample}.csv"
    conda:
        "../envs/abricate.yaml"
    params:
        "--db insaflu --minid 70 --mincov 60"
    shell:
        "abricate {params} {input} > {output}"

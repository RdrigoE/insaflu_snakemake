rule abricate_se:
    input:
        "samples/{sample}/spades_se/contigs.fasta",
    output:
        csv="samples/{sample}/abricate_se/abricate_{sample}.csv",
        yaml="samples/{sample}/abricate_se/abricate_{sample}.yaml",
    conda:
        "../envs/abricate.yaml"
    params:
        "--db insaflu --minid 70 --mincov 30",
    shell:
        "abricate {params} {input} > {output.csv}"
        " && python3 utils/get_abricate_type_notypes.py {output.csv} {output.yaml}"


rule abricate_pe:
    input:
        "samples/{sample}/spades_pe/contigs.fasta",
    output:
        csv="samples/{sample}/abricate_pe/abricate_{sample}.csv",
        yaml="samples/{sample}/abricate_pe/abricate_{sample}.yaml",
    conda:
        "../envs/abricate.yaml"
    params:
        "--db insaflu --minid 70 --mincov 30",
    shell:
        "abricate {params} {input} > {output.csv}"
        " && python3 utils/get_abricate_type_notypes.py {output.csv} {output.yaml}"


rule abricate_ont:
    input:
        "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
    output:
        csv="samples/{sample}/abricate_ont/abricate_{sample}.csv",
        yaml="samples/{sample}/abricate_ont/abricate_{sample}.yaml",
    conda:
        "../envs/abricate.yaml"
    params:
        "--db insaflu --minid 70 --mincov 30",
    shell:
        "abricate {params} {input} > {output.csv}"
        " && python3 utils/get_abricate_type_notypes.py {output.csv} {output.yaml}"

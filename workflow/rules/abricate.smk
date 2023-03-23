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
    resources:
        mem_mb=memory["abricate_se"],
    log:
        "logs/samples/{sample}/abricate.log",
    benchmark:
        "benchmark/samples/{sample}/abricate_se.tsv"
    shell:
        "abricate {params} {input} > {output.csv}"
        " && python3 {scripts_directory}get_abricate_info_list.py {output.csv} {output.yaml}"


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
    resources:
        mem_mb=memory["abricate_pe"],
    log:
        "logs/samples/{sample}/abricate.log",
    benchmark:
        "benchmark/samples/{sample}/abricate_pe.tsv"
    shell:
        "abricate {params} {input} > {output.csv}"
        " && python3 {scripts_directory}get_abricate_info_list.py {output.csv} {output.yaml}"


rule abricate_ont:
    input:
        "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
    output:
        csv="samples/{sample}/abricate_ont/abricate_{sample}.csv",
        yaml="samples/{sample}/abricate_ont/abricate_{sample}.yaml",
    conda:
        "../envs/abricate.yaml"
    resources:
        mem_mb=memory["abricate_ont"],
    log:
        "logs/samples/{sample}/abricate.log",
    benchmark:
        "benchmark/samples/{sample}/abricate_ont.tsv"
    params:
        "--db insaflu --minid 70 --mincov 30",
    shell:
        "abricate {params} {input} > {output.csv}"
        " && python3 {scripts_directory}get_abricate_info_list.py {output.csv} {output.yaml}"

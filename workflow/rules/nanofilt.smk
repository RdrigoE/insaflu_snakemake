rule nanofilt_SE:
    input:
        get_raw_input_ont,
    output:
        "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
    conda:
        "../envs/nanofilt.yaml"
    params:
        get_nanofilt_parameters(software_parameters),
    resources:
        mem_mb=memory["nanofilt_SE"],
    log:
        "logs/samples/{sample}/nanofilt.log",
    benchmark:
        "benchmark/samples/{sample}/nanofilt.tsv"
    shell:
        "gunzip -cd {input} | NanoFilt {params} | gzip > {output}"

configfile: "config/parameters.yaml"

rule nanofilt_SE:
    input:
        "user_data/{sample}.fastq.gz"
    output:
        o = "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz"
    conda:
        "../envs/nanofilt.yaml"
    params:
        "--quality 10 --length 50 --headcrop 70 --tailcrop 70 "
    shell:
        "gunzip -c {input} | NanoFilt {params} | gzip > {output.o}"
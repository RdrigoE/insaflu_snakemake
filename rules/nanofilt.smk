configfile: "config/parameters.yaml"

rule nanofilt_SE:
    input:
        "user_data/{sample}.fastq.gz"
    output:
        "samples/{sample}/nano_trimmed_reads/{sample}.trimmed.fastq.gz"
    conda:
        "../envs/nanofilt.yaml"
    threads: 
        config['trimmomatic_threads']
    params:
        "-q 10 -l 50 --headcrop 70 --tailcrop 70 --maxlength 0"
    shell:
        "gunzip -c {input} | NanoFilt {params} | gzip > {output}"
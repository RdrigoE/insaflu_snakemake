rule trimme_reads_SE:
    input:
        "user_data/{sample}."+ext
    output:
        "samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz"
    conda:
        "../envs/trimmomatic.yaml"
    params:
        "HEADCROP:30 SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33"
    shell:
        "trimmomatic SE "
        "{input} "
        "{output} " 
        "{params}"

rule trimme_reads_PE:
    input:
        i1 = "user_data/{sample}_1."+ext,
        i2 = "user_data/{sample}_2."+ext
    output:
        o1="samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        o2="samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",
        
        o_un1="samples/{sample}/trimmed_reads/{sample}_1.untrimmed.fastq.gz",
        o_un2="samples/{sample}/trimmed_reads/{sample}_2.untrimmed.fastq.gz"
    params:
        "HEADCROP:30 SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33"
    shell:
        "trimmomatic PE "
        "{input.i1} "
        "{input.i2} "
        "{output.o1} "
        "{output.o_un1} "
        "{output.o2} "
        "{output.o_un2} "
        "{params}"
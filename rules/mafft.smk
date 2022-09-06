rule mafft_pe:
    input:
        "results/concat_pe/multifile.fasta"
    output:
        "results/mafft_pe/mafft.fasta" 
    conda:
        "../envs/mafft.yaml"
    shell:
        "mafft --preservecase {input} > {output}"

rule mafft_se:
    input:
        "results/concat_se/multifile.fasta"
    output:
        "results/mafft_se/mafft.fasta" 
    conda:
        "../envs/mafft.yaml"
    shell:
        "mafft --preservecase {input} > {output}"
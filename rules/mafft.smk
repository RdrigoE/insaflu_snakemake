rule mafft_pe:
    input:
        "results/concat/multifile.fasta"
    output:
        "results/mafft_pe/mafft.fasta" 
    conda:
        "../envs/mafft.yaml"
    shell:
        "mafft --preservecase {input} > {output}"
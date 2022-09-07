rule mafft:
    input:
        "projects/{project}/concat/multifile.fasta"
    output:
        "projects/{project}/mafft/mafft.fasta" 
    conda:
        "../envs/mafft.yaml"
    params:
        "--preservecase"
    shell:
        "mafft {params} {input} > {output}"


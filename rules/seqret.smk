rule seqret:
    input:
        "projects/{project}/mafft/mafft.fasta" 

    output:
        "projects/{project}/seqret/seqret.nex"

    conda:
        "../envs/seqret.yaml"

    params:
        "-sformat fasta -osformat2 nexusnon"
    shell:
        "seqret {params} -sequence {input} -outseq {output}"
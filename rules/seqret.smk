rule seqret:
    input:
        "projects/{project}/main_result/mafft/mafft.fasta" 

    output:
        "projects/{project}/main_result/seqret/seqret.nex"

    conda:
        "../envs/seqret.yaml"

    params:
        "-sformat fasta -osformat2 nexusnon"
    shell:
        "seqret {params} -sequence {input} -outseq {output}"
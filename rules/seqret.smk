rule seqret:
    input:
        "results/mafft_pe/mafft.fasta" 

    output:
        "results/seqret/seqret.nex"

    conda:
        "../envs/seqret.yaml"
    shell:
        "seqret -sformat fasta -osformat2 nexusnon -sequence {input} -outseq {output}"
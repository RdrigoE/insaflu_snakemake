rule fasttree:
    input:
        "results/mafft_pe/mafft.fasta"    
    output:
        "results/fasttre/tree"    

    conda:
        "../envs/fasttree.yaml"
    shell:
        "fasttree -gtr -boot 1000 -nt {input} > {output}"
rule fasttree_pe:
    input:
        "results/mafft_pe/mafft.fasta"    
    output:
        "results/fasttre_pe/tree"    

    conda:
        "../envs/fasttree.yaml"
    shell:
        "fasttree -gtr -boot 1000 -nt {input} > {output}"

rule fasttree_se:
    input:
        "results/mafft_se/mafft.fasta"    
    output:
        "results/fasttre_se/tree"    

    conda:
        "../envs/fasttree.yaml"
    shell:
        "fasttree -gtr -boot 1000 -nt {input} > {output}"     
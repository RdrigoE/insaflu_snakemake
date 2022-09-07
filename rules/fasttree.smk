rule fasttree:
    input:
        "projects/{project}/mafft/mafft.fasta"    
    output:
        "projects/{project}/fasttre/tree"    

    conda:
        "../envs/fasttree.yaml"
    params:
        "-gtr -boot 1000 -nt"
    shell:
        "fasttree {params} {input} > {output}"

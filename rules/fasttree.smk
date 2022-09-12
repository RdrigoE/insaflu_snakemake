rule fasttree:
    input:
        "projects/{project}/main_result/mafft/mafft_masked.fasta"    
    output:
        "projects/{project}/main_result/fasttre/tree"    

    conda:
        "../envs/fasttree.yaml"
    params:
        "-gtr -boot 1000 -nt"
    shell:
        "fasttree {params} {input} > {output}"

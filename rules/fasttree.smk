rule fasttree:
    input:
        "projects/{project}/main_result/Alignment_nt_All.fasta"    
    output:
        "projects/{project}/main_result/Tree_ML_All.tree"    
    conda:
        "../envs/fasttree.yaml"
    params:
        "-gtr -boot 1000 -nt"
    shell:
        "fasttree {params} {input} > {output}"
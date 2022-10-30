rule fasttree:
    input:
        "projects/{project}/main_result/Alignment_nt_All.fasta"    
    output:
        tree = "projects/{project}/main_result/Tree_ML_All.tree",
        nwk = "projects/{project}/main_result/Tree_ML_All.nwk"
    conda:
        "../envs/fasttree.yaml"
    params:
        "-gtr -boot 1000 -nt"
    shell:
        #FastTreeDbl
        "fasttree {params} {input} > {output.tree} && cp {output.tree} {output.nwk}"
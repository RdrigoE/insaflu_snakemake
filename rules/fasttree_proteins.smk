rule fasttree_proteins:
    input:
        "projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_mafft.fasta"    
    output:
        tree = "projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_tree.tree", 
        nwk = "projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_tree.nwk" 
    conda:
        "../envs/fasttree.yaml"
    params:
        "-gtr -boot 1000"
    shell:
        #FastTreeDbls
        "fasttree {params} {input} > {output.tree} && cp {output.tree} {output.nwk}"
    

rule fasttree_proteins:
    input:
        "projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_mafft.fasta"    
    output:
        "projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_tree.tree" 
    conda:
        "../envs/fasttree.yaml"
    params:
        "-gtr -boot 1000 -nt"
    shell:
        "fasttree {params} {input} > {output}"
    

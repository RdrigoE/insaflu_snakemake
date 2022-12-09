# rule fasttree:
#     input:
#         "projects/{project}/main_result/Alignment_nt_All.fasta",
#     output:
#         tree="projects/{project}/main_result/Tree_ML_All.tree",
#         nwk="projects/{project}/main_result/Tree_ML_All.nwk",
#     conda:
#         "../envs/fasttree.yaml"
#     params:
#         "-gtr -boot 1000 -nt",
#     shell:
#         #FastTreeDbl
#         "fasttree {params} {input} > {output.tree} && cp {output.tree} {output.nwk}"


rule fasttree_proteins:
    input:
        "projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_mafft.fasta",
    output:
        tree="projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_tree.tree",
        nwk="projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_tree.nwk",
    conda:
        "../envs/fasttree.yaml"
    params:
        "-gtr -boot 1000",
    shell:
        #FastTreeDbls
        "fasttree {params} {input} > {output.tree} && cp {output.tree} {output.nwk}"


rule fasttree_nt:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft.fasta",
    output:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}_tree.tree",
    conda:
        "../envs/fasttree.yaml"
    params:
        "-gtr -boot 1000",
    shell:
        "fasttree {params} {input} > {output}"


# rule cp_Alignment_nt_tree:
#     input:
#         tree="projects/{project}/main_result/{seg}/Alignment_nt_{seg}_tree.tree",
#         # nwk="projects/{project}/main_result/Tree_ML_All.nwk",
#     output:
#         tree="projects/{project}/main_result/Alignment_nt_{seg}.tree",
#         nwk="projects/{project}/main_result/Alignment_nt_{seg}.nwk",
#     shell:
#         "cp {input.tree} {output.tree} &&"
#         "cp {input.tree} {output.nwk}"

rule concat_:
    input:
        expand("projects/{project}/samples/{sample}/snippy/snps.consensus.fa")
    output:
        o="projects/{project}/concat/multifile.fasta",
        dir = directory("projects/{project}/concat")   
    shell:
        "python merge_fasta_files.py {project} SARS_COV_2"
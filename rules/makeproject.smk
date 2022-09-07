rule makeproject:
    input:
        "samples/{sample}/snippy/snps.consensus.fa"
    output:
        directory("projects/{project}/sample_{sample}/")
    shell:
        "python3 utils/make_project.py {wildcards.project} {wildcards.sample}"
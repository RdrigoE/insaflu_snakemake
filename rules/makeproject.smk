rule makeproject:
    input:
        o1="samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        o2="samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",

    output:
        "projects/{project}/sample_{sample}/main_result/project.txt"
    shell:
        "mkdir projects/{wildcards.project}/main_result/ -p | "
        "mkdir projects/{wildcards.project}/sample_{wildcards.sample}/main_result -p |"
        " echo 'project created' > projects/{wildcards.project}/sample_{wildcards.sample}/main_result/project.txt"
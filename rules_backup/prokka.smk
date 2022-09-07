rule prokka:
    input:

    output:
        dir = directory("results/prokka")
    conda:
        "../envs/prokka.yaml"
    params: 
        "--kingdom Viruses --locustag locus --genus Influenzavirus --species Influenzavirus --strain "
        "ref_PREFIX_FILES_OUT --gcode 11" #https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

    shell:
        "prokka {params} --outdir {output.dir}"
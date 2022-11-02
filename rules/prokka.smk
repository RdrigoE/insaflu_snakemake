# not being used at the time of realease 1.0
# rule prokka:
#     input:
#     output:
#         dir = directory("results/prokka")
#     conda:
#         "../envs/prokka.yaml"
#     params:
#         "--kingdom Viruses --locustag locus --genus Influenzavirus --species Influenzavirus --strain "
#         "ref_PREFIX_FILES_OUT --gcode 11"
#     shell:
#         "prokka {params} --outdir {output.dir}"

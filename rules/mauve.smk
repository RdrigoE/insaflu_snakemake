# rule mauve_se:
#     input:
#         i = "results/snippy_se/{sample}_se/snps.consensus.fa",
#     output:
#         o = "results/mauve_se/{sample}.xmfa",
#         dir = directory("results/mauve_se")
#     conda:
#         "../envs/mauve.yaml"
#     shell:
#         "progressiveMauve --output={output.o} {input.i} "

SAMPLES, = glob_wildcards("results/snippy_pe/{sample}_pe/snps.consensus.fa")

rule mauve_pe:
    input:
        expand("results/snippy_pe/{sample}_pe/snps.consensus.fa", sample=SAMPLES)
    output:
        o = "results/mauve_pe/mauve.xmfa",
        dir = directory("results/mauve_pe")    
    conda:
        "../envs/mauve.yaml"
    shell:
        "progressiveMauve {input} --output={output.o}"
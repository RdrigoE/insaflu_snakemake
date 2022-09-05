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
#alterar o consensus de lugar para depois ler apenas todos os ficheiros da pasta ->Rodrigo LE isto
#falta fazer single end 
rule concat:
    input:
        expand("results/snippy_pe/{sample}_pe/snps.consensus.fa", sample=SAMPLES)
    output:
        o = "results/concat/multifile.fasta",
        dir = directory("results/concat")   
    shell:
        "cat {input} > {output.o}"
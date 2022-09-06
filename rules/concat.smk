SAMPLES_PE, = glob_wildcards("results/snippy_pe/{sample}_pe/snps.consensus.fa")

#alterar o consensus de lugar para depois ler apenas todos os ficheiros da pasta ->Rodrigo LE isto
#falta fazer single end 
rule concat_pe:
    input:
        expand("results/snippy_pe/{sample}_pe/snps.consensus.fa", sample=SAMPLES_PE)
    output:
        o="results/concat_pe/multifile.fasta",
        dir = directory("results/concat_pe")   
    shell:
        "python merge_fasta_files.py pe"

SAMPLES_SE, = glob_wildcards("results/snippy_pe/{sample}_se/snps.consensus.fa")

rule concat_se:
    input:
        expand("results/snippy_se/{sample}_se/snps.consensus.fa", sample=SAMPLES_SE)
    output:
        o="results/concat_se/multifile.fasta",
        dir = directory("results/concat_se")   
    shell:
        "python merge_fasta_files.py se"
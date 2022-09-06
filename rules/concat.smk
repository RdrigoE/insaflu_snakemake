SAMPLES, = glob_wildcards("results/snippy_pe/{sample}_pe/snps.consensus.fa")
#alterar o consensus de lugar para depois ler apenas todos os ficheiros da pasta ->Rodrigo LE isto
#falta fazer single end 
rule concat:
    input:
        expand("results/snippy_pe/{sample}_pe/snps.consensus.fa", sample=SAMPLES)
    output:
        o="results/concat/multifile.fasta",
        dir = directory("results/concat")   
    shell:
        "python merge_fasta_files.py"
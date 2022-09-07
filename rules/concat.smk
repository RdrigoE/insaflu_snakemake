checkpoint cp_directory:
    input:
        "projects/{project}/sample_{sample}/snippy/snps.consensus.fa"
    output:
        "projects/{project}/main_result/consensus/{sample}_SARS_COV_2_consensus.fasta"
    shell:
        "python utils/move_fasta_files.py {wildcards.project} {wildcards.sample} SARS_COV_2"

rule all_consensus:
    input:
        aggregate_input
    output:
        o="projects/{project}/main_result/AllConsensus.fasta"    
    shell:
        "cat projects/{wildcards.project}/main_result/consensus/* > projects/{wildcards.project}/main_result/AllConsensus.fasta"

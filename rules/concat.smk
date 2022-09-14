with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)


fr = REFERENCE
rule cp_directory:
    input:
        i1 = expand("projects/{project}/sample_{sample}/snippy/snps.consensus.fa",project=config_user['project'], sample=config_user['samples']),
    output:
        o1 = "projects/{project}/main_result/consensus/{sample}_consensus.fasta"
    shell:
        "python utils/move_fasta_files.py {wildcards.project} {wildcards.sample}" 

rule all_consensus:
    input:
        expand("projects/{project}/main_result/consensus/{sample}_consensus.fasta",project=config_user['project'], sample=config_user['samples'])
    output:
        o="projects/{project}/main_result/AllConsensus.fasta"    
    shell:
        "cat {fr} {input} > {output.o}"

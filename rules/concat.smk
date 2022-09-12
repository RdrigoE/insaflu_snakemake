with open('config_user/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

rule cp_directory:
    input:
        i1 = "projects/{project}/sample_{sample}/snippy/snps.consensus.fa"
    output:
        o1 = "projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta"
    shell:
        "python utils/move_fasta_files.py {wildcards.project} {wildcards.sample} SARS_COV_2" 

rule all_consensus:
    input:
        expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta",project=config_user['project'], sample=config_user['samples'])
    output:
        o="projects/{project}/main_result/AllConsensus.fasta"    
    shell:
        "cat reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.fasta {input} > {output.o}"

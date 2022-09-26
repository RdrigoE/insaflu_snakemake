with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

species = config_user['locus']

fr = REFERENCE
rule cp_directory:
    input:
        i1 = expand("projects/{project}/sample_{sample}/snippy/snps.consensus.fa",project=config_user['project'], sample=config_user['samples']),
    output:
        o1 = "projects/{project}/main_result/consensus/{sample}_consensus.fasta"
    shell:
        "python utils/move_fasta_files.py {wildcards.project} {wildcards.sample} {species}" 

rule all_consensus:
    input:
        i = expand("projects/{project}/main_result/consensus/{sample}_consensus.fasta",project=config_user['project'], sample=config_user['samples']),
        ref = REFERENCE_GB
    output:
        o="projects/{project}/main_result/AllConsensus.fasta",
        o_2 = "projects/{project}/main_result/All_nt.fasta"
    shell:
        "cat {fr} {input.i} > {output.o} | python utils/concat_segments.py {output.o} {input.ref} {species} {output.o_2}"

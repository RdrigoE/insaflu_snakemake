with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)


fr = REFERENCE
rule cp_directory:
    input:
        i1 = "projects/{project}/sample_{sample}/snippy/snps.consensus.fa",
    output:
        o1 = "projects/{project}/sample_{sample}/snippy/{sample}_consensus.fasta"
    shell:
        "python utils/get_consensus.py {input} {output}" 

rule all_consensus:
    input:
        i = expand("projects/{project}/sample_{sample}/snippy/{sample}_consensus.fasta",project=config_user['project'], sample=config_user['samples']),
        ref = REFERENCE_GB,
        coverage = "projects/{project}/main_result/coverage_translate.csv",
        fasta = REFERENCE
    output:
        o="projects/{project}/main_result/AllConsensus.fasta",
        o_2 = "projects/{project}/main_result/All_nt.fasta",
        o_3 = "projects/{project}/main_result/All_nt_only_90plus.fasta"
    shell:
        "cat {fr} {input.i} > {output.o} | python utils/concat_segments.py '{input.i}' {input.ref} {output.o_2} {input.coverage} {input.fasta} {output.o_3}"

rule all_consensus:
    input:
        i = expand("projects/{project}/sample_{sample}/{sample}_consensus.fasta",project = PROJECT_NAME, sample=config_user['samples']),
        coverage = "projects/{project}/main_result/coverage_translate.csv",
    output:
        o="projects/{project}/main_result/AllConsensus.fasta",
        all_consensus_no_ref = "projects/{project}/main_result/AllConsensus_no_ref.fasta",
        o_2 = "projects/{project}/main_result/All_nt.fasta",
        o_3 = "projects/{project}/main_result/All_nt_only_90plus.fasta"
    shell:
        "python utils/generate_AllConsensus.py {input.coverage} {REFERENCE_GB} '{input.i}' {REFERENCE} {output.o} {output.all_consensus_no_ref} && python utils/concat_segments.py '{input.i}' {REFERENCE_GB} {output.o_2} {input.coverage} {REFERENCE} {output.o_3}"
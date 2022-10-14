rule all_consensus:
    input:
        i = expand("projects/{project}/sample_{sample}/{sample}_consensus.fasta",project = PROJECT_NAME, sample=config_user['samples']),
        ref = REFERENCE_GB,
        coverage = "projects/{project}/main_result/coverage_translate.csv",
        fasta = REFERENCE
    output:
        o="projects/{project}/main_result/AllConsensus.fasta",
        all_consensus_no_ref = "projects/{project}/main_result/AllConsensus_no_ref.fasta",
        o_2 = "projects/{project}/main_result/All_nt.fasta",
        o_3 = "projects/{project}/main_result/All_nt_only_90plus.fasta"
    shell:
        "cat {REFERENCE} {input.i} > {output.o} && cat {input.i} > {output.all_consensus_no_ref} && python utils/concat_segments.py '{input.i}' {input.ref} {output.o_2} {input.coverage} {input.fasta} {output.o_3}"
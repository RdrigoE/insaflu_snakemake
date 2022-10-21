rule all_consensus:
    input:
        every_consensus = expand("projects/{project}/sample_{sample}/{sample}_consensus.fasta",project = PROJECT_NAME, sample=config_user['samples']),
        coverage = "projects/{project}/main_result/coverage_translate.csv",
    output:
        AllConsensus = "projects/{project}/main_result/AllConsensus.fasta",
        all_consensus_no_ref = "projects/{project}/main_result/AllConsensus_no_ref.fasta",
        All_nt = "projects/{project}/main_result/All_nt.fasta",
        All_nt_only_90plus = "projects/{project}/main_result/All_nt_only_90plus.fasta"
    shell:
        "python utils/generate_AllConsensus.py {input.coverage} {REFERENCE_GB} '{input.every_consensus}' {REFERENCE} {output.AllConsensus} {output.all_consensus_no_ref} && python utils/concat_segments.py '{input.every_consensus}' {REFERENCE_GB} {output.All_nt} {input.coverage} {REFERENCE} {output.All_nt_only_90plus}"
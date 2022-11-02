rule create_segments:
    input:
        expand(
            "projects/{project}/main_result/Alignment_nt_All_sep.fasta",
            project=PROJECT_NAME,
        ),
    output:
        expand(
            "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
            project=PROJECT_NAME,
            seg=SEGMENTS,
        ),
    shell:
        "python utils/pre_msa_masker.py {input} projects/{PROJECT_NAME}/main_result {REFERENCE_GB}"

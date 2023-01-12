rule initiate_folder:
    output:
        dir("projects/{project}/main_result/"),
    shell:
        "mkdir {output} && cp {user_metadata_directort}parameters.yaml projects/{wildcards.project}/"


rule makeproject:
    input:
        get_consensus,
    output:
        consensus="projects/{project}/sample_{sample}/{sample}_consensus.fasta",
        depth="projects/{project}/sample_{sample}/snps.depth.gz",
        vcf="projects/{project}/sample_{sample}/snps.vcf",
    params:
        get_directory,
    shell:
        "mkdir projects/{wildcards.project}/sample_{wildcards.sample}/ -p && "
        " cp -r {params} projects/{wildcards.project}/sample_{wildcards.sample}/ "


rule assemble_consensus:
    input:
        every_consensus=expand(
            "projects/{project}/sample_{sample}/{sample}_consensus.fasta",
            project=PROJECT_NAME,
            sample=config_user["samples"],
        ),
        coverage="projects/{project}/main_result/coverage_translate.csv",
    output:
        AllConsensus="projects/{project}/main_result/AllConsensus.fasta",
        all_consensus_no_ref="projects/{project}/main_result/AllConsensus_no_ref.fasta",
        All_nt="projects/{project}/main_result/All_nt.fasta",
        All_nt_only_90plus="projects/{project}/main_result/All_nt_only_90plus.fasta",
    shell:
        "python {scripts_directory}generate_AllConsensus.py {input.coverage} {REFERENCE_GB} '{input.every_consensus}' {REFERENCE_FASTA} {output.AllConsensus} {output.all_consensus_no_ref} "
        "&& python {scripts_directory}concat_segments.py '{input.every_consensus}' {REFERENCE_GB} {output.All_nt} {input.coverage} {REFERENCE_FASTA} {output.All_nt_only_90plus}"


rule create_segments:
    input:
        actual_input=expand(
            "projects/{project}/main_result/Alignment_nt_All_sep.fasta",
            project=PROJECT_NAME,
        ),
    output:
        temp(
            expand(
                "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
                project=PROJECT_NAME,
                seg=SEGMENTS,
            )
        ),
    shell:
        "python {scripts_directory}split_files_by_locus.py {input} projects/{PROJECT_NAME}/main_result {REFERENCE_GB}"

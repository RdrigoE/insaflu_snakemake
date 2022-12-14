configfile: "../config/threads.yaml"


#########################################################################################
# Assemblers


rule align_w_nextalign_medaka:
    input:
        "align_samples/{sample}/medaka/medaka_align_{seg}.fasta",
    output:
        output=temp("align_samples/{sample}/medaka/medaka_aligned_{seg}.fasta"),
        temp_ref=temp("align_samples/{sample}/medaka/temp_ref_file_{seg}"),
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    params:
        "-j {threads}",
    shell:
        "python {scripts_directory}remove_ref.py {input} {output.temp_ref} && "
        "nextalign run --include-reference -r {output.temp_ref} {params} --output-fasta {output.output} {input} "
        # "python {scripts_directory}attach_ref.py {input} {output.temp_ref}"


rule align_nextalign_snippy:
    input:
        align_file="align_samples/{sample}/snippy/snippy_align_{seg}.fasta",
    output:
        output="align_samples/{sample}/snippy/snippy_aligned_{seg}.fasta",
        temp_ref=temp("align_samples/{sample}/snippy/temp_ref_file_{seg}.fasta"),
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    params:
        "-j 12",
    shell:
        "python {scripts_directory}remove_ref.py {input} {output.temp_ref} && "
        "nextalign run --include-reference -r {output.temp_ref} {params} --output-fasta {output.output} {input} "
        # "python {scripts_directory}attach_ref.py {output} {output.temp_ref}"


#######################################################################


rule nextalign_pre_aa:
    input:
        "projects/{project}/main_result/AllConsensus.fasta",
    output:
        temp("projects/{project}/main_result/Alignment_nt_All_sep.fasta"),
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    # params:
    #     "-j {threads}",
    shell:
        "nextalign run -r {REFERENCE_FASTA} -j {threads} --output-fasta {output} {input}"


rule nextalign:
    input:
        "projects/{project}/main_result/All_nt.fasta",
    output:
        output="projects/{project}/main_result/Alignment_nt_All.fasta",
        temp_ref=temp("projects/{project}/main_result/temp_ref_file"),
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    # params:
    #     "-j {threads}",
    shell:
        "python {scripts_directory}remove_ref.py {input} {output.temp_ref} && "
        "nextalign run -r {output.temp_ref} -j {threads} --output-fasta {output.output} {input} "
        # "python {scripts_directory}attach_ref.py {input} {output.temp_ref}"


rule nextalign_nt:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
    output:
        output="projects/{project}/main_result/{seg}/Alignment_nt_{seg}_mafft.fasta",
        temp_ref=temp("projects/{project}/main_result/{seg}/temp_ref_file"),
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    # params:
    #     "-j {threads}",
    shell:
        "python {scripts_directory}remove_ref.py {input} {output.temp_ref} && "
        "nextalign run --include-reference -r {output.temp_ref} -j {threads} --output-fasta {output.output} {input} "
        # "python {scripts_directory}attach_ref.py {input} {output.temp_ref}"


rule nextalign_proteins:
    input:
        "projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_trans.fasta",
    output:
        output="projects/{project}/main_result/{locus}/Alignment_aa_{locus}_{gene}_mafft.fasta",
        temp_ref=temp("projects/{project}/main_result/{locus}/temp_ref_file_{gene}"),
    conda:
        "../envs/nextalign.yaml"
    threads: config["nextalign_threads"]
    # params:
    #     "-j {threads}",
    shell:
        "python {scripts_directory}remove_ref.py {input} {output.temp_ref} && "
        "nextalign run --include-reference -r {output.temp_ref} -j {threads} --output-fasta {output.output} {input} "
        # "python {scripts_directory}attach_ref.py {input} {output.temp_ref}"


rule cp_Alignment_nt_nextalign:
    input:
        "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
    output:
        "projects/{project}/main_result/Alignment_nt_{seg}.fasta",
    shell:
        "cp {input} {output}"

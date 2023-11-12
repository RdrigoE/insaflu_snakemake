def get_output_sample(
    paired_illumina: list[str],
    single_illumina: list[str],
    ont_samples: list[str],
    classification: bool,
    **kwargs,
) -> list[str]:
    files = [
        expand(
            "samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html",
            sample=paired_illumina,
            direction=["1", "2"],
        ),
        expand(
            "samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html",
            sample=paired_illumina,
            direction=["1", "2"],
        ),
        expand(
            "samples/{sample}/raw_fastqc/{sample}_fastqc.html",
            sample=single_illumina,
        ),
        expand(
            "samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html",
            sample=single_illumina,
        ),
        expand(
            "samples/{sample}/raw_nanostat/{sample}_stats.txt",
            sample=ont_samples,
        ),
        expand(
            "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
            sample=ont_samples,
        ),
        expand(
            "samples/{sample}/nano_trimmed_fastqc/{sample}_stats.txt",
            sample=ont_samples,
        ),
        expand(
            "samples/{sample}/raw_rabbit_qc/{sample}_stats.html",
            sample=ont_samples,
        ),
        expand(
            "samples/{sample}/processed_rabbit_qc/{sample}_stats.html",
            sample=ont_samples,
        ),
        expand(
            "samples/{sample}/raw_rabbit_qc/{sample}_stats.html",
            sample=ont_samples,
        ),
        expand(
            "samples/{sample}/processed_rabbit_qc/{sample}_stats.html",
            sample=ont_samples,
        ),
    ]
    if classification:
        files.extend(
            [
                expand(
                    "samples/{sample}/spades_se/contigs.fasta",
                    sample=single_illumina,
                ),
                expand(
                    "samples/{sample}/abricate_se/abricate_{sample}.csv",
                    sample=single_illumina,
                ),
                expand(
                    "samples/{sample}/abricate_se/abricate_{sample}.yaml",
                    sample=single_illumina,
                ),
                expand(
                    "samples/{sample}/spades_pe/contigs.fasta",
                    sample=paired_illumina,
                ),
                expand(
                    "samples/{sample}/abricate_pe/abricate_{sample}.csv",
                    sample=paired_illumina,
                ),
                expand(
                    "samples/{sample}/abricate_pe/abricate_{sample}.yaml",
                    sample=paired_illumina,
                ),
                expand(
                    "samples/{sample}/abricate_ont/abricate_{sample}.csv",
                    sample=ont_samples,
                ),
                expand(
                    "samples/{sample}/abricate_ont/abricate_{sample}.yaml",
                    sample=ont_samples,
                ),
            ]
        )
    return files


def get_output_project(
    paired_illumina: list[str],
    single_illumina: list[str],
    ont_samples: list[str],
    classification: bool,
    illumina_vc_software: str,
    project_name: str,
    min_coverage: int,
    reference_gb: str,
) -> list[str]:
    files = [
        get_output_sample(
            paired_illumina, single_illumina, ont_samples, classification
        ),
        expand(
            "align_samples/{sample}/{illumina_genome_assembly_software}/{sample}_consensus.fasta",
            illumina_genome_assembly_software=illumina_assembler,
            sample=single_illumina,
        ),
        expand(
            "align_samples/{sample}/{illumina_genome_assembly_software}/{sample}_consensus.fasta",
            illumina_genome_assembly_software=illumina_assembler,
            sample=paired_illumina,
        ),
        expand(
            "align_samples/{sample}/{sample}_coverage.csv",
            sample=single_illumina,
        ),
        expand(
            "align_samples/{sample}/{sample}_coverage.csv",
            sample=paired_illumina,
        ),
        expand(
            "projects/{project}/main_result/coverage_translate.csv",
            project=config_user["project"],
        ),
        expand(
            "samples/{sample}/raw_nanostat/{sample}_stats.txt",
            sample=ont_samples,
        ),
        expand(
            "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
            sample=ont_samples,
        ),
        expand(
            "samples/{sample}/nano_trimmed_fastqc/{sample}_stats.txt",
            sample=ont_samples,
        ),
        expand(
            "samples/{sample}/raw_rabbit_qc/{sample}_stats.html",
            sample=ont_samples,
        ),
        expand(
            "samples/{sample}/processed_rabbit_qc/{sample}_stats.html",
            sample=ont_samples,
        ),
        Checkpoint_Main(
            [
                (f'projects/{config_user["project"]}/main_result/', "_trans.fasta"),
                (f'projects/{config_user["project"]}/main_result/', "_mafft.fasta"),
                (f'projects/{config_user["project"]}/main_result/', "_mafft.nex"),
                (f'projects/{config_user["project"]}/main_result/', "_tree.tree"),
            ],
            reference_gb,
            f"projects/{config_user['project']}/main_result/coverage_translate.csv",
            min_coverage,
            project_name,
            samples,
            SEGMENTS
        ),
        Abricate_Pangolin(
            f"projects/{config_user['project']}/ref.yaml",
            f"projects/{config_user['project']}/main_result/lineage_report.csv",
            project_name,
        ),
    ]
    return files

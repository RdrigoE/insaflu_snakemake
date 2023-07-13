rule iVar_align_pe:
    input:
        reads_1="samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        reads_2="samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz",
    output:
        bam="align_samples/{sample}/iVar/pre_snps.bam",
    conda:
        "../envs/ivar.yaml"
    resources:
        mem_mb=memory["iVar_align_pe"],
    log:
        "logs/align_samples/{sample}/iVar/mapping.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/mapping_pe.tsv"
    params:
        minqual=software_parameters["mapqual"],
    shell:
        "mkdir align_samples/{wildcards.sample}/reference -p && "
        "cp {REFERENCE_FASTA} align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta && "
        "bwa index align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta && "
        "bwa mem align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta {input.reads_1} {input.reads_2} | "
        "samtools view -u -T align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta -q {params.minqual} | "
        "samtools sort --reference align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta > {output}"


rule iVar_align_se:
    input:
        reads_1="samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz",
    output:
        "align_samples/{sample}/iVar/pre_snps.bam",
    conda:
        "../envs/ivar.yaml"
    resources:
        mem_mb=memory["iVar_align_se"],
    log:
        "logs/align_samples/{sample}/iVar/mapping.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/mapping_se.tsv"
    params:
        minqual=software_parameters["mapqual"],
    shell:
        "mkdir align_samples/{wildcards.sample}/reference -p && "
        "cp {REFERENCE_FASTA} align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta && "
        "bwa index align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta && "
        "bwa mem align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta {input.reads_1} | "
        "samtools view -u -T align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta -q {params.minqual} | "
        "samtools sort --reference align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta > {output}"


ruleorder: iVar_align_pe > iVar_align_se


rule primers_bam:
    input:
        pre_snps="align_samples/{sample}/iVar/pre_snps.bam",
    output:
        snps="align_samples/{sample}/iVar/second_snps.bam",
    conda:
        "../envs/ivar.yaml"
    resources:
        mem_mb=memory["primers_bam"],
    log:
        "logs/align_samples/{sample}/iVar/primers.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/primers.tsv"
    params:
        folder="align_samples/{sample}/iVar/",
    shell:
        "bwa mem -k 5 -T 16 align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta {PRIMER_FASTA} | samtools view -b -F 4 > {params.folder}primers.bam"
        " && bedtools bamtobed -i {params.folder}primers.bam > {params.folder}primers.bed &&"
        "samtools sort -o {input.pre_snps}.sorted {input.pre_snps} "
        "&& ivar trim -b {params.folder}primers.bed -i {input.pre_snps}.sorted -m 0 -q 0 -p {output} -e "
        "&& samtools sort -o {output} {output}"


rule injection:
    input:
        bam="align_samples/{sample}/iVar/second_snps.bam",
    output:
        snps="align_samples/{sample}/iVar/snps.bam",
    conda:
        "../envs/ivar.yaml"
    params:
        prefix="align_samples/{sample}/iVar/temp_snps",
        folder="align_samples/{sample}/iVar/",
    resources:
        mem_mb=memory["primers_bam"],
    log:
        "logs/align_samples/{sample}/iVar/primers.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/primers.tsv"
    shell:
        " ivar trim -m 0 -q 0 -e -b {params.folder}primers.bed -p {params.prefix}.trimmed -i {input.bam}"
        " && samtools sort -o {params.prefix}.trimmed.sorted.bam {params.prefix}.trimmed.bam"
        " && samtools index {params.prefix}.trimmed.sorted.bam"
        " && samtools mpileup -A -d 0 -Q 0 {input.bam} | ivar consensus -m 0 -n N -p {params.prefix}.ivar_consensus"
        " && ./../workflow/scripts/run_check_consensus {params.prefix}.ivar_consensus.fa align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta"
        " && bwa index -p {params.prefix}.ivar_consensus {params.prefix}.ivar_consensus.fa"
        " && bwa mem -k 5 -T 16 {params.prefix}.ivar_consensus {PRIMER_FASTA} | samtools view -bS -F 4 | samtools sort -o {params.folder}primers_consensus.bam"
        " && samtools mpileup -A -d 0 --reference {params.prefix}.ivar_consensus.fa -Q 0 {params.folder}primers_consensus.bam | ivar variants -p {params.folder}primers_consensus -t 0.03"
        " && bedtools bamtobed -i {params.folder}primers_consensus.bam > {params.folder}primers_consensus.bed"
        " && ivar getmasked -i {params.folder}primers_consensus.tsv -b {params.folder}primers_consensus.bed -f {PRIMER_FASTA}.pair_information.tsv -p {params.folder}primer_mismatchers_indices"
        " && ivar removereads -i {params.prefix}.trimmed.sorted.bam -p {params.prefix}.masked.bam -t {params.folder}primer_mismatchers_indices.txt -b {params.folder}primers.bed"
        " && samtools sort -o {params.prefix}.masked.sorted.bam {params.prefix}.masked.bam"
        " && cp {params.prefix}.masked.sorted.bam {output}"
        " && samtools index {output}"


rule call_variant:
    input:
        bam="align_samples/{sample}/iVar/snps.bam",
    output:
        tsv_file="align_samples/{sample}/iVar/snps.tsv",
    conda:
        "../envs/ivar.yaml"
    resources:
        mem_mb=memory["call_variant"],
    log:
        "logs/align_samples/{sample}/iVar/call_variants.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/call_variants.tsv"
    params:
        mapqual=software_parameters["mapqual"],
        mincov=software_parameters["mincov"],
        minfrac=software_parameters["minfrac"],
    shell:
        "samtools mpileup -A -d 600000 -B -Q 0 {input.bam} --reference align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta | "
        "ivar variants -p align_samples/{wildcards.sample}/iVar/snps -q {params.mapqual} -t {params.minfrac} -m {params.mincov}  align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta "


rule generate_consensus:
    input:
        bam="align_samples/{sample}/iVar/snps.bam",
    output:
        consensus="align_samples/{sample}/iVar/snps.consensus.fa",
    conda:
        "../envs/ivar.yaml"
    params:
        consensus="snps.consensus",
        mapqual=software_parameters["mapqual"],
        mincov=software_parameters["mincov"],
        minfrac=software_parameters["minfrac"],
    resources:
        mem_mb=memory["generate_consensus"],
    log:
        "logs/align_samples/{sample}/iVar/generate_consensus.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/generate_consensus.tsv"
    shell:
        "samtools mpileup -d 0 -A -Q {params.mapqual} {input.bam} --reference align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta | "
        "ivar consensus -p align_samples/{wildcards.sample}/iVar/{params.consensus} -q {params.mapqual} -t {params.minfrac} -n N -m {params.mincov} "


rule get_depth:
    input:
        i="align_samples/{sample}/iVar/snps.bam",
        consensus="align_samples/{sample}/iVar/snps.consensus.fa",
    output:
        only_depth="align_samples/{sample}/iVar/snps.depth.gz",
        depth="align_samples/{sample}/iVar/snps.depth",
        ind="align_samples/{sample}/iVar/snps.depth.gz.tbi",
    conda:
        "../envs/snippy.yaml"
    params:
        "-aa -q 20",
    resources:
        mem_mb=memory["get_depth"],
    log:
        "logs/align_samples/{sample}/iVar/get_depth.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/get_depth.tsv"
    shell:
        "samtools depth {params} {input.i} | bgzip -c > {output.only_depth} "
        "&& tabix -p vcf {output.only_depth} && gunzip -c {output.only_depth}  > {output.depth} "


rule iVar_depth_1_2:
    input:
        consensus="align_samples/{sample}/iVar/snps.consensus.fa",
        depth="align_samples/{sample}/iVar/snps.depth",
    output:
        consensus="align_samples/{sample}/iVar/new_snps.consensus.fa",
    conda:
        "../envs/base.yaml"
    resources:
        mem_mb=memory["iVar_depth_1_2"],
    log:
        "logs/align_samples/{sample}/iVar/depth_1_2.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/depth_1_2.tsv"
    shell:
        "python ../workflow/scripts/split_consensus.py {input.consensus} {input.depth} {REFERENCE_GB} {output.consensus}"


rule iVar_depth_step_2:
    input:
        zipped="align_samples/{sample}/iVar/snps.depth",
    output:
        unzipped="align_samples/{sample}/iVar/{seg}.depth",
    conda:
        "../envs/base.yaml"
    resources:
        mem_mb=memory["iVar_depth_step_2"],
    log:
        "logs/align_samples/{sample}/iVar/depth_step_2/{seg}.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/depth_step_2/{seg}.tsv"
    shell:
        "python ../workflow/scripts/split_depth_file.py {input.zipped} {REFERENCE_GB}"


rule create_align_file_iVar:
    input:
        first_consensus="align_samples/{sample}/iVar/new_snps.consensus.fa",
    output:
        align_file=temp("align_samples/{sample}/iVar/iVar_align_{seg}.fasta"),
    conda:
        "../envs/base.yaml"
    resources:
        mem_mb=memory["create_align_file_iVar"],
    log:
        "logs/align_samples/{sample}/iVar/create_align_file/{seg}.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/create_align_file/{seg}.tsv"
    shell:
        "python ../workflow/scripts/mask_consensus_by_deep.py align_samples/{wildcards.sample}/reference/{REFERENCE_NAME}.fasta {input.first_consensus} {output.align_file} {wildcards.seg}"


rule align_mafft_iVar:
    input:
        align_file="align_samples/{sample}/iVar/iVar_align_{seg}.fasta",
    output:
        aligned_file=temp("align_samples/{sample}/iVar/iVar_aligned_{seg}.fasta"),
    conda:
        "../envs/mafft.yaml"
    threads: config["mafft_threads"]
    params:
        "--preservecase",
    resources:
        mem_mb=memory["align_mafft_iVar"],
    log:
        "logs/align_samples/{sample}/iVar/align_mafft/{seg}.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/align_mafft/{seg}.tsv"
    shell:
        "mafft --thread {threads} {params} {input.align_file} > {output.aligned_file}"


rule msa_masker_iVar:
    input:
        align_file="align_samples/{sample}/iVar/iVar_aligned_{seg}.fasta",
        depth="align_samples/{sample}/iVar/{seg}.depth",
    output:
        temp("align_samples/{sample}/iVar/consensus_aligned_{seg}.fasta"),
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c " + str(software_parameters["mincov"] - 1),
    resources:
        mem_mb=memory["msa_masker_iVar"],
    log:
        "logs/align_samples/{sample}/iVar/msa_masker/{seg}.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/msa_masker/{seg}.tsv"
    shell:
        "python ../workflow/scripts/msa_masker.py -i {input.align_file} -df {input.depth} -o {output} {params}"


rule get_masked_consensus_iVar:
    input:
        lambda wildcards: expand(
            "align_samples/{sample}/iVar/consensus_aligned_{seg}.fasta",
            sample=wildcards.sample,
            seg=get_locus(REFERENCE_GB),
        ),
    output:
        final_consensus="align_samples/{sample}/iVar/pre_{sample}_consensus.fasta",
    conda:
        "../envs/base.yaml"
    resources:
        mem_mb=memory["get_masked_consensus_iVar"],
    log:
        "logs/align_samples/{sample}/iVar/get_masked_consensus.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/get_masked_consensus.tsv"
    shell:
        "python ../workflow/scripts/get_consensus_medaka.py '{input}' {output}"


rule mask_regions_consensus_iVar:
    input:
        consensus="align_samples/{sample}/iVar/pre_{sample}_consensus.fasta",
        snps="align_samples/{sample}/iVar/snps.vcf",
        depth="align_samples/{sample}/iVar/snps.depth.gz",
    output:
        final_consensus="align_samples/{sample}/iVar/{sample}_consensus.fasta",
    conda:
        "../envs/base.yaml"
    params:
        mask_regions_parameters(software_parameters),
    resources:
        mem_mb=memory["mask_regions_consensus_iVar"],
    log:
        "logs/align_samples/{sample}/iVar/mask_regions_consensus.log",
    benchmark:
        "benchmark/align_samples/{sample}/iVar/mask_regions_consensus.tsv"
    shell:
        "python {scripts_directory}mask_regions.py {input.consensus} {output.final_consensus} {params} "


rule get_vcf:
    input:
        tsv_file="align_samples/{sample}/iVar/snps.tsv",
    output:
        vcf_file="align_samples/{sample}/iVar/snps.vcf",
    resources:
        mem_mb=memory["get_vcf"],
    conda:
        "../envs/base.yaml"
    shell:
        "python {scripts_directory}convert_vcf.py {input.tsv_file} {output.vcf_file} "

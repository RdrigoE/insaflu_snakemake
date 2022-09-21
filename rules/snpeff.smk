with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)
locus = config_user['locus']
rule prepare_snpeff:
    input:
        ref = REFERENCE_GB,
        i=expand("projects/{project}/sample_{sample}/snippy/snps.vcf",sample=config_user['samples'], project=config_user['project']),
    output:
        "projects/{project}/main_result/snpeff/ready.txt"
    conda:
        "../envs/snpeff.yaml"
    shell:
        "python utils/create_snpeff_text.py $CONDA_PREFIX {input.ref} {locus} {REFERENCE_NAME} |"
        "conda activate $CONDA_PREFIX | "
        "mkdir  config/data/{REFERENCE_NAME} -p |"
        "cat {input.ref} > config/data/{REFERENCE_NAME}/genes.gbk | "
        "snpEff build -genbank {REFERENCE_NAME} -c config/snpeff.config > {output}"

rule snpeff:
    input:
        i = "projects/{project}/main_result/freebayes/{sample}_var.vcf",
        i2 = "projects/{project}/main_result/snpeff/ready.txt"
    output:
        o = "projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",
    conda:
        "../envs/snpeff.yaml"
    params:
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats -c config/snpeff.config"
    shell:
        #"snpEff {params} -v MN908947.3 {input.i}  > {output.o}"
        "snpEff {params} -v {REFERENCE_NAME} {input.i}  > {output.o}"

#ir ao snpEff config e colocar:❯ nano /home/reusebio/tese/insaflu_snakemake/.snakemake/conda/42a5a38a1b1e62c38dd21415a318b1ee_/share/snpeff-4.3.1t-5/snpEff.config 
    # # SARS CoV 2
    #     sarscov2.genome :  Severe_acute_respiratory_syndrome_coronavirus_2
    #             sarscov2.chromosomes : SC2.1
    #             sarscov2.MN908947.codonTable : Bacterial_and_Plant_Plastid
    #
    #tentar novamente a apontar para a ref
#snpEff build -c $refdir/snpeff.config -dataDir . -gff3 ref"


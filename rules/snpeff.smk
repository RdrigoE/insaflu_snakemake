with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)


rule snpeff:
    input:
        i = "projects/{project}/main_result/freebayes/{sample}_var.vcf"
    output:
        o = "projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",
    conda:
        "../envs/snpeff.yaml"
    params:
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats -c config/snpeff.config"
    shell:
        "snpEff {params} -v sarscov2 {input.i}  > {output.o}"


#ir ao snpEff config e colocar:❯ nano /home/reusebio/tese/insaflu_snakemake/.snakemake/conda/42a5a38a1b1e62c38dd21415a318b1ee_/share/snpeff-4.3.1t-5/snpEff.config 
    # # SARS CoV 2
    #     sarscov2.genome :  Severe acute respiratory syndrome coronavirus 2
    #             sarscov2.chromosomes : SC2.1
    #             sarscov2.MN908947.codonTable : Bacterial_and_Plant_Plastid
    #
    #tentar novamente a apontar para a ref
#snpEff build -c $refdir/snpeff.config -dataDir . -gff3 ref"


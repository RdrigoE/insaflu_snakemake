rule getCoverage_se:
#get all coverage files here...
    input: 
        i  =  "results/snippy_se/{sample}_se/snps.depth.gz",
        _1 = "results/snippy_se/{sample}_se/snps.bam",
        _2 = "results/snippy_se/{sample}_se/snps.tab",
        _3 = "results/snippy_se/{sample}_se/snps.consensus.fa",
        ref = REFERENCE
    output:
        o = "results/coverage_se/{sample}_coverage.tab"
    conda:
        "../envs/coverage.yaml"
    shell:
        "python software/getCoverage/getCoverage.py -i {input.i} -r {input.ref} -o {output.o}"

rule getCoverage_pe:
    input: 
        i  =  "results/snippy_pe/{sample}_pe/snps.depth.gz",
        _1 = "results/snippy_pe/{sample}_pe/snps.bam",
        _2 = "results/snippy_pe/{sample}_pe/snps.tab",
        _3 = "results/snippy_pe/{sample}_pe/snps.consensus.fa",
        ref = REFERENCE
    output:
        o = "results/coverage_pe/{sample}_coverage.tab"
    conda:
        "../envs/coverage.yaml"
    shell:
        "python software/getCoverage/getCoverage.py -i {input.i} -r {input.ref} -o {output.o}"
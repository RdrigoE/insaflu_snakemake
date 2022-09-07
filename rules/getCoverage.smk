rule getCoverage:
#get all coverage files here...
    input: 
        i  = "projects/{project}/sample_{sample}/snippy/snps.depth.gz",
        _0 = "projects/{project}/sample_{sample}/snippy/snps.bam",
        _2 = "projects/{project}/sample_{sample}/snippy/snps.tab",
        _3 = "projects/{project}/sample_{sample}/snippy/snps.consensus.fa",
        ref = REFERENCE
    output:
        o = "projects/{project}/main_result/coverage/{sample}_coverage.tab"
    conda:
        "../envs/coverage.yaml"
    shell:
        "python software/getCoverage/getCoverage.py -i {input.i} -r {input.ref} -o {output.o}"

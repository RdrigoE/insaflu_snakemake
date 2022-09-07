rule getCoverage:
#get all coverage files here...
    input: 
        i  = "projects/{project}/samples/snippy/snps.depth.gz",
        _1 = "projects/{project}/samples/snippy/snps.bam",
        _2 = "projects/{project}/samples/snippy/snps.tab",
        _3 = "projects/{project}/samples/snippy/snps.consensus.fa",
        ref = REFERENCE
    output:
        o = "projects/{project}/coverage/{sample}_coverage.tab"
    conda:
        "../envs/coverage.yaml"
    shell:
        "python software/getCoverage/getCoverage.py -i {input.i} -r {input.ref} -o {output.o}"

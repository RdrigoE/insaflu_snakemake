rule getCoverage:
#get all coverage files here...
    input: 
        i  = "align_samples/{sample}/snippy/snps.depth.gz",
        ref = REFERENCE
    output:
        o = "projects/{project}/main_result/coverage/{sample}_coverage.csv"
    conda:
        "../envs/coverage.yaml"
    shell:
        "python software/getCoverage/getCoverage.py -i {input.i} -r {input.ref} -o {output.o}"

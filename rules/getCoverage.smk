
def get_consensus(wildcards):
    return f"align_samples/{wildcards.sample}/{config_user['sample_type'][wildcards.sample]}/snps.depth.gz"


rule getCoverage:
#get all coverage files here...
    input: 
        i  = get_consensus,
        ref = REFERENCE
    output:
        o = temp("projects/{project}/main_result/{sample}_coverage.csv")
    conda:
        "../envs/coverage.yaml"
    shell:
        "python software/getCoverage/getCoverage.py -i {input.i} -r {input.ref} -o {output.o}"

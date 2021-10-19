rule freebayes:
    input:
        ref=config["reference"]["ref"],
        samples="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        index="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
    output:
        temp("variantCalls/callers/freebayes/{sample}_{seqID}.freebayes.unsort.vcf"),  # either .vcf or .bcf
    log:
        "logs/variantCalling/freebayes/{sample}_{seqID}.log",
    singularity:
        config["singularitys"]["freebayes"]
    params:
        extra=(
            " --min-alternate-fraction 0.01 --allele-balance-priors-off --pooled-discrete --pooled-continuous --report-genotype-likelihood-max -t "
            + config["bed"]["bedfile"]
        ),
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
    threads: 1
    shell:
        "(freebayes {params.extra} -f {input.ref} {input.samples} > {output}) &> {log}"


rule sortFreebayes:
    input:
        "variantCalls/callers/freebayes/{sample}_{seqID}.freebayes.unsort.vcf",
    output:
        temp("variantCalls/callers/freebayes/{sample}_{seqID}.freebayes.vcf"),
    singularity:
        config["singularitys"]["bcftools"]
    log:
        "logs/variantCalling/freebayes/{sample}_{seqID}.sort.log",
    shell:
        "(bcftools sort -o {output} -O v {input}) &> {log}"

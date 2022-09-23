rule bgzip_fixAF:
    input:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.weirdAF.vcf",
    output:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.weirdAF.vcf.gz",
        tabix="variantCalls/callers/{method}/{sample}_{seqID}.{method}.weirdAF.vcf.gz.tbi",
    log:
        "logs/variantCalling/bgzip/{method}/{sample}_{seqID}.weirdAF.log",
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input.vcf} && tabix {output.vcf}) &> {log}"

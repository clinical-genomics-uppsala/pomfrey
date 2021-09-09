localrules:
    bgzipCallers,


rule bgzipCallers:
    input:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf",
    output:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz",
        tabix="variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz.tbi",
    log:
        "logs/variantCalling/bgzip/{method}/{sample}_{seqID}.log",
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input.vcf} && tabix {output.vcf}) &> {log}"

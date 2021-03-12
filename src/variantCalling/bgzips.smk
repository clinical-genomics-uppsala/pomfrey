localrules: bgzipCallers
rule bgzipCallers:
    input:
        vcf = "variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf" #[m+"/{sample}."+m+".vcf" for m in config["methods"]]
    output:
        vcf = temp("variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf.gz"),
        tabix = temp("variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf.gz.tbi")
    log:
        "logs/variantCalling/bgzip/{method}/{sample}_{seqID}.log"
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input.vcf} && tabix {output.vcf}) 2> {log}"

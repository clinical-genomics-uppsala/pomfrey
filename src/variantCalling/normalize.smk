localrules:
    indexDecomp,


rule decompose:
    input:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf.gz",
        tbi="variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf.gz.tbi",
    output:
        temp("variantCalls/callers/{method}/{sample}_{seqID}.{method}.decomposed.vcf.gz"),
    log:
        "logs/variantCalling/vt/{sample}_{seqID}.{method}.decomposed.log",
    singularity:
        config["singularitys"]["vt"]
    shell:
        "(vt decompose -s {input.vcf} | vt decompose_blocksub -o {output} -) &> {log}"


rule normalizeAll:
    input:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.decomposed.vcf.gz",
        fasta=config["reference"]["ref"],  #,
        # tbi = "variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf.gz.tbi"
    output:
        "variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz",
    log:
        "logs/variantCalling/vt/{sample}_{seqID}.{method}.normalized.log",
    singularity:
        config["singularitys"]["vt"]
    shell:
        "(vt normalize -n -r {input.fasta} -o {output} {input.vcf} ) &> {log}"


rule indexDecomp:
    input:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz",
    output:
        tbi="variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz.tbi",
    log:
        "logs/variantCalling/vt/{sample}_{seqID}.{method}.index.log",
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(tabix {input.vcf}) 2> {log}"

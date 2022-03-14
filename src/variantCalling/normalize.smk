localrules:
    indexNormalized,


rule decompose:
    input:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf",
    output:
        temp("variantCalls/callers/{method}/{sample}_{seqID}.{method}.decomposed.vcf.gz"),
    log:
        "logs/variantCalling/vt/{sample}_{seqID}.{method}.decomposed.log",
    container:
        config["singularitys"]["vt"]
    shell:
        "(vt decompose -s {input.vcf} | vt decompose_blocksub -o {output} -) &> {log}"


rule normalizeAll:
    input:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.decomposed.vcf.gz",
        fasta=config["reference"]["ref"],
    output:
        "variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz",
    log:
        "logs/variantCalling/vt/{sample}_{seqID}.{method}.normalized.log",
    container:
        config["singularitys"]["vt"]
    shell:
        "(vt normalize -n -r {input.fasta} -o {output} {input.vcf} ) &> {log}"


rule indexNormalized:
    input:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz",
    output:
        tbi="variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz.tbi",
    log:
        "logs/variantCalling/vt/{sample}_{seqID}.{method}.index.log",
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(tabix {input.vcf}) 2> {log}"

localrules: indexDecomp

rule decompose: #Do we need decompose as well, maybe for all but vardict??
    input:
        vcf = "variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf.gz",  #[m+"/{sample}_{seqID}."+m+".normalized.vcf.gz" for m in config["methods"]] ##inte normalized.vcf filer! Hur?!
        tbi = "variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf.gz.tbi"
    output:
        temp("variantCalls/callers/{method}/{sample}_{seqID}.{method}.decomposed.vcf.gz")
    log:
        "logs/variantCalling/vt/{sample}_{seqID}.{method}.decomposed.log"
    container:
        config["singularitys"]["vt"]
    shell:
        "(vt decompose -s {input.vcf} | vt decompose_blocksub -o {output} -) &> {log}"

rule normalizeAll:
    input:
        vcf = "variantCalls/callers/{method}/{sample}_{seqID}.{method}.decomposed.vcf.gz", #"variantCalls/callers/{method}/{sample}.{method}.vcf.gz", #[m+"/{sample}."+m+".vcf.gz" for m in config["methods"]], ##inte normalized.vcf filer! Hur?!
        fasta = config["reference"]["ref"] #,
        # tbi = "variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf.gz.tbi"
    output:
        "variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz"
    log:
        "logs/variantCalling/vt/{sample}_{seqID}.{method}.normalized.log"
    container:
        config["singularitys"]["vt"]
    shell:
        "(vt normalize -n -r {input.fasta} -o {output} {input.vcf} ) &> {log}"

rule indexDecomp:
    input:
        vcf = "variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz" #"variantCalls/callers/{method}/{sample}.{method}.decomposed.vcf.gz" #[m+"/{sample}."+m+".vcf" for m in config["methods"]]
    output:
        tbi = "variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf.gz.tbi" #"variantCalls/callers/{method}/{sample}.{method}.decomposed.vcf.gz.tbi"
    log:
        "logs/variantCalling/vt/{sample}_{seqID}.{method}.index.log"
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(tabix {input.vcf}) 2> {log}"

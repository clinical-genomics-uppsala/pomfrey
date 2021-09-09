localrules:
    indexSnver,
    concatSnver,
    indexSnverIndel,


rule snver:
    input:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        index="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
        ref=config["reference"]["ref"],
        bed=config["bed"]["bedfile"],
    output:
        temp("variantCalls/callers/snver/{sample}_{seqID}.snver.raw.vcf"),
        temp("variantCalls/callers/snver/{sample}_{seqID}.snver.filter.vcf"),
        temp("variantCalls/callers/snver/{sample}_{seqID}.snver.indel.raw.vcf"),
        temp("variantCalls/callers/snver/{sample}_{seqID}.snver.indel.filter.vcf"),
    params:
        outfolder="variantCalls/callers/snver/{sample}_{seqID}.snver",
    log:
        "logs/variantCalling/snver/{sample}_{seqID}.log",
    singularity:
        config["singularitys"]["snver"]
    shell:
        "(snver -i {input.bam} -r {input.ref} -l {input.bed} -o {params.outfolder}) &> {log}"


rule indexSnver:
    input:
        "variantCalls/callers/snver/{sample}_{seqID}.snver.filter.vcf",
    output:
        temp("variantCalls/callers/snver/{sample}_{seqID}.snver.filter.vcf.gz.tbi"),
        "variantCalls/callers/snver/{sample}_{seqID}.snver.filter.vcf.gz",
    log:
        "logs/variantCalling/snver/{sample}_{seqID}.index.log",
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input} && tabix {input}.gz) &> {log}"


rule indexSnverIndel:
    input:
        "variantCalls/callers/snver/{sample}_{seqID}.snver.indel.filter.vcf",
    output:
        temp("variantCalls/callers/snver/{sample}_{seqID}.snver.indel.filter.vcf.gz.tbi"),
        "variantCalls/callers/snver/{sample}_{seqID}.snver.indel.filter.vcf.gz",
    log:
        "logs/variantCalling/snver/{sample}_{seqID}.index.indel.log",
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input} && tabix {input}.gz) &> {log}"


rule concatSnver:
    input:
        snver="variantCalls/callers/snver/{sample}_{seqID}.snver.filter.vcf.gz",
        indel="variantCalls/callers/snver/{sample}_{seqID}.snver.indel.filter.vcf.gz",
        index="variantCalls/callers/snver/{sample}_{seqID}.snver.indel.filter.vcf.gz.tbi",
        index2="variantCalls/callers/snver/{sample}_{seqID}.snver.filter.vcf.gz.tbi",
    output:
        temp("variantCalls/callers/snver/{sample}_{seqID}.snver.weirdAF.vcf"),
    log:
        "logs/variantCalling/snver/concat_{sample}_{seqID}.log",
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools concat -a -Ou {input.snver} {input.indel} | bcftools sort -Ov -o {output} -) &> {log}"

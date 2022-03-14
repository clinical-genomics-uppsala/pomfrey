localrules:
    fixAF,


include: "freebayes.smk"
include: "vardict_T.smk"
include: "pisces.smk"
include: "mutect2.smk"


rule fixAF:
    input:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.weirdAF.vcf",
    output:
        temp("variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf"),
    params:
        config["programdir"]["dir"],
    log:
        "logs/variantCalling/fixAF/{method}/{sample}_{seqID}.log",
    container:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params}/src/variantCalling/fix_af.py {input.vcf} {output}) &> {log}"


include: "normalize.smk"
include: "recall.smk"
include: "vep.smk"

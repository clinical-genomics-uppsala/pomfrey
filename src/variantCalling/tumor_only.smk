localrules:
    fixAF,


include: "freebayes.smk"
include: "vardict_T.smk"
include: "pisces.smk"
include: "mutect2.smk"


rule fixAF:
    input:
        "variantCalls/callers/{method}/{sample}_{seqID}.{method}.weirdAF.vcf",
    output:
        temp("variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf"),
    params:
        config["programdir"]["dir"],
    log:
        "logs/variantCalling/fixAF/{method}/{sample}_{seqID}.log",
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params}/src/variantCalling/fix_af.py {input} {output}) &> {log}"


include: "bgzips.smk"
include: "normalize.smk"
include: "recall.smk"
include: "vep.smk"

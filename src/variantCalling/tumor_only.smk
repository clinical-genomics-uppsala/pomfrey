localrules:
    fixAF,


include: "freebayes.smk"
include: "vardict_T.smk"
include: "pisces.smk"
include: "mutect2.smk"
include: "normalize.smk"

rule fixAF:
    input:
        vcf = "variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.weirdAF.vcf.gz",
        tbi = "variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.weirdAF.vcf.gz.tbi",
    output:
        temp("variantCalls/callers/{method}/{sample}_{seqID}.{method}.normalized.vcf"),
    params:
        config["programdir"]["dir"],
    log:
        "logs/variantCalling/fixAF/{method}/{sample}_{seqID}.log",
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params}/src/variantCalling/fix_af.py {input.vcf} {output}) &> {log}"

include: "bgzips.smk"
include: "recall.smk"
include: "vep.smk"

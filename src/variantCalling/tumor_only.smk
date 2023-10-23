localrules:
    fixAF,


include: "freebayes.smk"
include: "vardict_T.smk"
include: "pisces.smk"
include: "mutect2.smk"
include: "bgzips.smk"


rule fixAF:
    input:
        vcf="variantCalls/callers/{method}/{sample}_{seqID}.{method}.weirdAF.vcf.gz",
        tbi="variantCalls/callers/{method}/{sample}_{seqID}.{method}.weirdAF.vcf.gz.tbi",
    output:
        vcf=temp("variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf"),
    log:
        "logs/variantCalling/fixAF/{method}/{sample}_{seqID}.log",
    container:
        config["singularitys"]["python"]
    script:
        "fix_af.py"


include: "normalize.smk"
include: "recall.smk"
include: "vep.smk"

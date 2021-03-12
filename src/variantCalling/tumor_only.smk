localrules: fixAF
include:    "freebayes.smk"

# include:    "lofreq_T.smk"

# include:    "snver.smk" #since based on samtools cannot handle high base qualities, everything becomes N

include:    "vardict_T.smk"

include:    "pisces.smk"

include:    "mutect2.smk"
#include:    "manta_T.smk"
rule fixAF:
    input:
        "variantCalls/callers/{method}/{sample}_{seqID}.{method}.weirdAF.vcf"
    output:
        temp("variantCalls/callers/{method}/{sample}_{seqID}.{method}.vcf")
    params:
        config["programdir"]["dir"]
    log:
        "logs/variantCalling/fixAF/{method}/{sample}_{seqID}.log"
    container:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params}/src/variantCalling/fix_af.py {input} {output}) &> {log}"


include:    "bgzips.smk"

include:    "normalize.smk"


include:    "recall.smk"

include:    "vep.smk"

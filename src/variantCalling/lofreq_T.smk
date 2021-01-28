rule lofreq:
    input:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        bai="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
        ref=config["reference"]["ref"],
    output:
        temp("variantCalls/callers/lofreq/{sample}_{seqID}.lofreq.vcf"),
    log:
        "logs/variantCalling/lofreq_call/{sample}_{seqID}.log",
    params:
        extra="-l " + config["bed"]["bedfile"],
    singularity:
        config["singularitys"]["lofreq"]
    threads: 8
    shell:
        "(lofreq call-parallel --pp-threads {threads} -f {input.ref} {input.bam} -o {output} {params.extra} ) &> {log}"

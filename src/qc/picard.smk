rule picardHsMetrics:
    input:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        intervals=config["bed"]["intervals"],  ##Create with gatk ..
    output:
        "qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt",
    params:
        covCap=5000,
    log:
        "logs/qc/picardHsMetrics/{sample}_{seqID}.log",
    container:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk CollectHsMetrics -COVERAGE_CAP {params.covCap} -BAIT_INTERVALS {input.intervals} -TARGET_INTERVALS {input.intervals} "
        "-INPUT {input.bam} -OUTPUT {output}) &> {log}"

localrules:
    touchBatch,


rule samtools_stats:
    input:
        "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
    output:
        "qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt",
    params:
        extra="-t " + config["bed"]["bedfile"],
    log:
        "logs/qc/samtools_stats/{sample}_{seqID}.log",
    container:
        config["singularitys"]["bwa"]
    shell:
        "(samtools stats {params.extra} {input} > {output} ) &> {log}"


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


rule touchBatch:
    input:
        expand(
            "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
    output:
        temp("Results/batchQC_{seqID}/{seqID}_stats_unsorted.csv"),
    log:
        "logs/touch_{seqID}.log",
    shell:
        "(touch {output}) &> {log}"


rule getStatsforMqc:
    input:
        picardDup="qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt",
        picardMet="qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt",
        samtools="qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt",
        cartool="qc/{sample}_{seqID}/{sample}_{seqID}_Log.csv",
        batch="Results/batchQC_{seqID}/{seqID}_stats_unsorted.csv",
    output:
        batchTmp=temp("qc/{sample}_{seqID}/{sample}_batchStats.done"),
    params:
        dir=config["programdir"]["dir"],
    log:
        "logs/qc/{sample}_{seqID}_stats.log",
    container:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params.dir}/src/qc/get_stats.py {input.picardDup} {input.picardMet} {input.samtools} "
        "{input.cartool} {wildcards.sample} {input.batch} && touch {output.batchTmp}) &> {log}"


rule sortBatchStats:
    input:
        SampleSheet="SampleSheet.csv",
        batchUnsorted="Results/batchQC_{seqID}/{seqID}_stats_unsorted.csv",
        batchDone=expand(
            "qc/{sample}_{seqID}/{sample}_batchStats.done", sample=config["samples"], seqID=config["seqID"]["sequencerun"]
        ),
    output:
        batch="Results/batchQC_{seqID}/{seqID}_stats_mqc.json",
    params:
        dir=config["programdir"]["dir"],
        cov=config["cartool"]["cov"],
    log:
        "logs/qc/sortBatchStats_{seqID}.log",
    container:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params.dir}/src/qc/sortBatchStats.py {input.batchUnsorted} {input.SampleSheet} {output.batch} "
        "{params.cov}) &> {log}"

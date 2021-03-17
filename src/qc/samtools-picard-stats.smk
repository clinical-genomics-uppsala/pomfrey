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
    singularity:
        config["singularitys"]["bwa"]
    shell:
        "(samtools stats {params.extra} {input} > {output} ) &> {log}"


rule picardHsMetrics:
    input:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        intervals=config["bed"]["intervals"],  ##Create with gatk ..
    output:
        "qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt",
    log:
        "logs/qc/picardHsMetrics/{sample}_{seqID}.log",
    singularity:
        config["singularitys"]["bwa"]
    shell:
        "(java -Xmx4g -jar /opt/conda/share/picard-2.20.1-0/picard.jar CollectHsMetrics BAIT_INTERVALS={input.intervals} \
                TARGET_INTERVALS={input.intervals} INPUT={input.bam} OUTPUT={output}) &> {log}"


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
        # multiQCheader=config["programdir"]["dir"] + "src/qc/multiqc-header.txt",
        cartool="qc/{sample}_{seqID}/{sample}_{seqID}_Log.csv",
        batch="Results/batchQC_{seqID}/{seqID}_stats_unsorted.csv",
    output:
        batchTmp=temp("qc/{sample}_{seqID}/{sample}_batchStats.done"),
        # batch = "qc/{seqID}_stats_mqc.tsv",
        # sample="qc/{sample}_{seqID}/{sample}_{seqID}_stats_mqc.csv",
    params:
        dir=config["programdir"]["dir"],
    log:
        "logs/qc/{sample}_{seqID}_stats.log",
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params.dir}/src/qc/get_stats.py {input.picardDup} {input.picardMet} {input.samtools} \
                    {input.cartool} {wildcards.sample} {input.batch} && touch {output.batchTmp}) &> {log}" # # {input.multiQCheader} {output.sample}


rule sortBatchStats:
    input:
        SampleSheetUsed="SampleSheet.csv",
        batchUnsorted="Results/batchQC_{seqID}/{seqID}_stats_unsorted.csv",
        batchDone=expand(
            "qc/{sample}_{seqID}/{sample}_batchStats.done", sample=config["samples"], seqID=config["seqID"]["sequencerun"]
        ),
    output:
        batch="Results/batchQC_{seqID}/{seqID}_stats_mqc.json",
    params:
        dir=config["programdir"]["dir"],
    log:
        "logs/qc/sortBatchStats_{seqID}.log",
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params.dir}/src/qc/sortBatchStats.py {input.batchUnsorted} {input.SampleSheet} {output.batch}) &> {log}"

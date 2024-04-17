rule getStatsforMqc:
    input:
        picard_dup=expand(
            "qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        picard_metrics=expand(
            "qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        samtools_stats=expand(
            "qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        mosdepth_summary=expand(
            "qc/mosdepth/{sample}_{seqID}.mosdepth.summary.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        mosdepth_thresh_summary=expand(
            "qc/mosdepth/{sample}_{seqID}.thresholds_summary.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        sample_order="{seqID}_order.tsv",
    output:
        json="multiqc/batchQC_{seqID}/{seqID}_stats_mqc.json",
    params:
        thresholds=config["misc"]["cov_thresholds"],
    log:
        "logs/qc/{seqID}_stats.log",
    container:
        config["singularitys"]["python"]
    script:
        "get_stats.py"

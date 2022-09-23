rule getStatsforMqc:
    input:
        picard_dup=expand("qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        picard_metrics=expand("qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        samtools_stats=expand("qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        mosdepth_summary=expand("qc/mosdepth/{sample}_{seqID}.mosdepth.summary.txt", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        mosdepth_thresh_summary=expand("qc/mosdepth/{sample}_{seqID}.thresholds_summary.txt", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        samplesheet="SampleSheet.csv",
    output:
        json="Results/batchQC_{seqID}/{seqID}_stats_mqc.json",
    params:
        thresholds=config["misc"]["cov_thresholds"],
    log:
        "logs/qc/{seqID}_stats.log",
    container:
        config["singularitys"]["python"]
    script:
        "get_stats.py"

#
# rule sortBatchStats:
#     input:
#         SampleSheet=
#         batchDone=expand(
#             "qc/{sample}_{seqID}/{sample}_batchStats.done", sample=config["samples"], seqID=config["seqID"]["sequencerun"]
#         ),
#     output:
#         batch="Results/batchQC_{seqID}/{seqID}_stats_mqc.json",
#     params:
#         dir=config["programdir"]["dir"],
#         cov=config["misc"]["cov_thresholds"],
#         batchUnsorted="Results/batchQC_{seqID}/{seqID}_stats_unsorted.csv",
#     log:
#         "logs/qc/sortBatchStats_{seqID}.log",
#     container:
#         config["singularitys"]["python"]
#     shell:
#         "(python3.6 {params.dir}/src/qc/sortBatchStats.py {params.batchUnsorted} {input.SampleSheet} {output.batch} "
#         "{params.cov}) &> {log}"

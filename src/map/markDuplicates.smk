rule markDuplicates:
    input:
        bam="data_processing/{sample}_{seqID}/{sample}_{seqID}.bam",
        bai="data_processing/{sample}_{seqID}/{sample}_{seqID}.bam.bai",
    output:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        metric="qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt",
    log:
        "logs/map/{sample}_{seqID}-dedup.log",
    threads: 5
    container:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk MarkDuplicates -INPUT {input.bam} -OUTPUT {output.bam} -METRICS_FILE {output.metric}) &> {log} "


rule samtools_index_dedup:
    input:
        "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
    output:
        "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
    log:
        "logs/map/samtools_index/{sample}_{seqID}-dedup.log",
    container:
        config["singularitys"]["bwa"]
    shell:
        "(samtools index {input} {output}) &> {log}"

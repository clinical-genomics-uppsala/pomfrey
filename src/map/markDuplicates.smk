rule markDuplicates:
    input:
        bam="data_processing/{sample}_{seqID}/{sample}_{seqID}.bam",
        bai="data_processing/{sample}_{seqID}/{sample}_{seqID}.bam.bai",
    output:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        metric="qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt",
    log:
        "logs/map/{sample}_{seqID}-dedup.log",
    threads: 4  ##2??
    singularity:
        config["singularitys"]["bwa"]
    shell:
        "(java -Xmx4g -jar /opt/conda/share/picard-2.20.1-0/picard.jar MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} "
        "METRICS_FILE={output.metric}) &> {log}"


rule samtools_index_dedup:
    input:
        "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
    output:
        "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
    log:
        "logs/map/samtools_index/{sample}_{seqID}-dedup.log",  # optional params string
    singularity:
        config["singularitys"]["bwa"]
    shell:
        "(samtools index {input} {output}) &> {log}"

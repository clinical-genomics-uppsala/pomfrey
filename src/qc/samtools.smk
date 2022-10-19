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

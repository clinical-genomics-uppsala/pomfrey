rule mosdepth:
    input:
        bam = "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        bai = "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
        bed = config["bed"]["bedfile"],
    output:
        "qc/mosdepth/{sample}_{seqID}/{sample}_{seqID}.mosdepth.global.dist.txt",
        "qc/mosdepth/{sample}_{seqID}/{sample}_{seqID}.mosdepth.region.dist.txt",
        "qc/mosdepth/{sample}_{seqID}/{sample}_{seqID}.per-base.bed.gz",
        "qc/mosdepth/{sample}_{seqID}/{sample}_{seqID}.per-base.bed.gz.csi",
        "qc/mosdepth/{sample}_{seqID}/{sample}_{seqID}.regions.bed.gz",
        "qc/mosdepth/{sample}_{seqID}/{sample}_{seqID}.regions.bed.gz.csi",
        summary="qc/mosdepth/{sample}_{seqID}/{sample}_{seqID}.mosdepth.summary.txt",
    log:
        "logs/qc/mosdepth_{sample}_{seqID}.log",
    threads:
        4
    singularity:
        "docker://hydragenetics/mosdepth:0.3.2"
    wrapper:
        "v1.12.0/bio/mosdepth"

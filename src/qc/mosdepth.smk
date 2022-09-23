rule mosdepth:
    input:
        bam = "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        bai = "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
        bed = config["bed"]["bedfile"],
    output:
        glob="qc/mosdepth/{sample}_{seqID}.mosdepth.global.dist.txt",
        region="qc/mosdepth/{sample}_{seqID}.mosdepth.region.dist.txt",
        coverage="qc/mosdepth/{sample}_{seqID}.per-base.bed.gz",
        coverage_csi="qc/mosdepth/{sample}_{seqID}.per-base.bed.gz.csi",
        bed="qc/mosdepth/{sample}_{seqID}.regions.bed.gz",
        bed_csi="qc/mosdepth/{sample}_{seqID}.regions.bed.gz.csi",
        summary="qc/mosdepth/{sample}_{seqID}.mosdepth.summary.txt",
        thresholds="qc/mosdepth/{sample}_{seqID}.thresholds.bed.gz",
    params:
        thresholds=config["misc"]["cov_thresholds"]
    log:
        "logs/qc/mosdepth_{sample}_{seqID}.log",
    threads:
        4
    singularity:
        config["singularitys"]["mosdepth"]
    wrapper:
        "v1.12.0/bio/mosdepth"


rule fix_thresholds_file:
    input:
        mosdepth_thresholds="qc/mosdepth/{sample}_{seqID}.thresholds.bed.gz",
    output:
        fixed="qc/mosdepth/{sample}_{seqID}.thresholds_fixed.bed",
        summary="qc/mosdepth/{sample}_{seqID}.thresholds_summary.txt",
    log:
        "qc/mosdepth/{sample}_{seqID}.thresholds_fixed.bed.log",
    singularity:
        config["singularitys"]["python"]
    script:
        "fix_mosdepth_threshold.py"

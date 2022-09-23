rule bedtools_intersect_lowcov_regions:
    input:
        left="qc/mosdepth/{sample}_{seqID}.per-base.bed.gz",
        right=config["bed"]["bedfile"],
    output:
        vcf="qc/mosdepth/{sample}_{seqID}.mosdepth.lowCov.regions.txt",
    log:
        "logs/bedtools_intersect/{sample}_{seqID}_lowcov.log"
    singularity:
        config["singularitys"]["bedtools"]
    wrapper:
        "v1.3.1/bio/bedtools/intersect"


rule bedtools_intersect_hotspot:
    input:
        left="qc/mosdepth/{sample}_{seqID}.per-base.bed.gz",
        right=config["bed"]["hotspot"],
    output:
        vcf="qc/mosdepth/{sample}_{seqID}.mosdepth.hotspots.txt",
    log:
        "logs/bedtools_intersect/{sample}_{seqID}_hotspots.log"
    singularity:
        config["singularitys"]["bedtools"]
    wrapper:
        "v1.3.1/bio/bedtools/intersect"

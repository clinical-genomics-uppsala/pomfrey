rule all:
    input:
        expand(
            "Results/{sample}_{seqID}/Data/{sample}_{seqID}.indel.bam",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "Results/{sample}_{seqID}/Reports/{sample}_{seqID}.xlsx",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "Results/{sample}_{seqID}/Data/{sample}_{seqID}.SNV-pindel.vcf",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "Results/{sample}_{seqID}/Data/{sample}_{seqID}.normalized.genome.vcf.gz",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "Results/{sample}_{seqID}/Data/{sample}_{seqID}.normalized.genome.vcf.gz.tbi",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand("Results/multiqc_{seqID}.html", seqID=config["seqID"]["sequencerun"]),
        expand(
            "CNV/{sample}_{seqID}_clean.calledCNVs.modeled.png",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),


wildcard_constraints:
    seqID=config["seqID"]["sequencerun"],


### QC modules
include: "qc/fastqc.smk"  #fastq in html/text out
include: "qc/samtools.smk"  #bam in txt out
include: "qc/picard.smk"  #bam in txt out
include: "qc/batch_stats.smk"  #bam in txt out
include: "qc/mosdepth.smk"
include: "qc/bedtools.smk"
## Trimming in runfolder/{sample}_S[0-9]_R[12]_001.fastq.gz out trimming/{sample}_R[12]_trimmed.fastq.gz
include: "trimming/cutadapt.smk"
## Map in trimming/{sample}_R[12]_trimmed.fastq.gz out bam/{sample}.bam
include: "map/bwa-mem.smk"  #fastq R1 R2 from trimming in, bam out.
include: "map/markDuplicates.smk"
## Variant callers
include: "variantCalling/tumor_only.smk"
include: "variantCalling/pindel.smk"
## CNV?
include: "CNV/run_GATK_CNV.smk"
include: "CNV/cnvkit.smk"
## Rapportering
include: "report/multiqc.smk"  # per sample, add per batch as well but only certain results?
include: "report/vcf2excel.smk"
include: "report/igv-images.smk"  #per sample

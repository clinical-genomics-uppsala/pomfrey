##Specify configfile and singularity folder in snakemake command.
from snakemake.utils import min_version
min_version("6.6")

## PGX samples based on sample sheet
with open(config["samplesheet"], "r") as f:
    lines = f.readlines()
    samples = []
    head_found = False

    for start, line in enumerate(lines):
        if line.strip() == "[Data]":
            break

    for i, head in lines[start + 1].strip().split():
        if head == "PGX_Analysis":
            break
            head_found = True

    if not head_found:
        raise KeyError("No pgx_analysis head found in samplesheet!")
    else:
        for line in lines[start + 2:]:
            words = line.strip().split()
            if words[i] == "1":
                samples.append(words[0])

config["pgx_sample"] = [s for s in samples if s in config["sample"]]

rule All:
    input:
        expand(
            "Results/{sample}_{seqID}/Data/{sample}_{seqID}.indel.bam",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        # expand(
        #     "Results/{sample}_{seqID}/Reports/{sample}_{seqID}_MultiQC.html",
        #     sample=config["samples"],
        #     seqID=config["seqID"]["sequencerun"],
        # ),
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
            "variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz", sample=config["samples"], seqID=config["seqID"]["sequencerun"]
        ),
        expand("qc/{sample}_{seqID}/{sample}_batchStats.done", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        expand("Results/batchQC_{seqID}/{seqID}_MultiQC.html", seqID=config["seqID"]["sequencerun"]),
        expand(
            "CNV/{sample}_{seqID}_clean.calledCNVs.modeled.png", sample=config["samples"], seqID=config["seqID"]["sequencerun"]
        ),expand(
             "work/{seqID}/Results/Report/{sample}_{seqID}_pgx.html",
                sample=config["pgx_samples"],
                seqID=config["seqID"]["sequencerun"])


wildcard_constraints:
    # sample = "[a-zA-Z0-9-_\.]+",
    # support = "3", #"\.[0-9]+\."
    seqID=config["seqID"]["sequencerun"],


### QC modules
include: "qc/fastqc.smk"  #fastq in html/text out
include: "qc/samtools-picard-stats.smk"  #bam in txt out
include: "qc/cartool.smk"  #bam in tables out
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
## Rapportering
include: "report/multiqc.smk"  # per sample, add per batch as well but only certain results?
include: "report/vcf2excel.smk"
include: "report/igv-images.smk"  #per sample

## PGX module import
module pgx_workflow:
    snakefile: "../pgx_module/Snakefile_pgx.smk"
    config: config

use rule * from pgx_workflow as pgx_*

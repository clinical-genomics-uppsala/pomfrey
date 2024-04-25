rule multiqcBatch:
    input:
        expand(
            "qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "qc/mosdepth/{sample}_{seqID}.mosdepth.summary.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "qc/mosdepth/{sample}_{seqID}.mosdepth.region.dist.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "qc/mosdepth/{sample}_{seqID}.mosdepth.global.dist.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.zip",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "data_processing/{sample}_{seqID}/{sample}_{seqID}.qc.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        "multiqc/batchQC_{seqID}/{seqID}_stats_mqc.json",
    output:
        "multiqc/batchQC_{seqID}/multiqc_DNA.html",
    params:
        extra="-c " + config["configCache"]["multiqc"] + " --ignore *_{seqID}_stats_mqc.csv",
        output_dir="multiqc/batchQC_{seqID}",
        output_name="multiqc_DNA.html",
    log:
        "logs/report/multiqc/{seqID}.log",
    container:
        config["singularitys"]["multiqc"]
    shell:
        "( multiqc {params.extra} --force -o {params.output_dir} -n {params.output_name} {input} ) &> {log}"


rule multiqc_cp:
    input:
        html="multiqc/batchQC_{seqID}/multiqc_DNA.html",
    output:
        html="Results/multiqc_DNA.html",
    log:
        "logs/report/multiqc_cp_{seqID}.log",
    shell:
        "cp -r {input.html} {output.html}"

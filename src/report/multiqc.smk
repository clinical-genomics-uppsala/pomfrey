# rule multiqc:
#     input:
#         [
#             "qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt",
#             "qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.zip",
#             "data_processing/{sample}_{seqID}/{sample}_{seqID}.qc.txt",
#             "qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt",
#             "qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt",
#         ],
#     output:
#         "Results/{sample}_{seqID}/Reports/{sample}_{seqID}_MultiQC.html",
#     params:
#         extra="-c " + config["configCache"]["multiqc"],
#         output_dir="Results/{sample}_{seqID}/Reports/",
#         output_name="{sample}_{seqID}_MultiQC.html",
#     log:
#         "logs/report/multiqc/{sample}_{seqID}.log",
#     singularity:
#         config["singularitys"]["multiqc"]
#     shell:
#         "( multiqc {params.extra} --force -o {params.output_dir} -n {params.output_name} {input} ) &> {log}"


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
            "qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.zip",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        expand(
            "data_processing/{sample}_{seqID}/{sample}_{seqID}.qc.txt",
            sample=config["samples"],
            seqID=config["seqID"]["sequencerun"],
        ),
        "Results/batchQC_{seqID}/{seqID}_stats_mqc.json",
        expand("qc/{sample}_{seqID}/{sample}_batchStats.done", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),  #Wait until all in table
    output:
        "Results/batchQC_{seqID}/{seqID}_MultiQC.html",
    params:
        extra="-c " + config["configCache"]["multiqc"] + " --ignore *_{seqID}_stats_mqc.csv",
        output_dir="Results/batchQC_{seqID}",
        output_name="{seqID}_MultiQC.html",
    log:
        "logs/report/multiqc/{seqID}.log",
    singularity:
        config["singularitys"]["multiqc"]
    shell:
        "( multiqc {params.extra} --force -o {params.output_dir} -n {params.output_name} {input} ) &> {log}"

rule fastqcR1:
    input:
        "data_processing/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed.fastq.gz",
    output:
        html="qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.html",
        zip="qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.zip",
    params:
        outdir="qc/{sample}_{seqID}/",
    log:
        "logs/qc/fastqc/{sample}_{seqID}_R1_trimmed.log",
    singularity:
        config["singularitys"]["fastqc"]
    shell:
        "(fastqc --quiet --outdir {params.outdir} {input}) &> {log}"

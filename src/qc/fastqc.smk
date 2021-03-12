rule fastqcR1:
    input:
        "data_processing/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed.fastq.gz",
    output:
        html="qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.html",
        zip="qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        outdir = "qc/{sample}_{seqID}/"
    log:
        "logs/qc/fastqc/{sample}_{seqID}_R1_trimmed.log",
    container:
        config["singularitys"]["fastqc"]
    shell:
        "(fastqc --quiet --outdir {params.outdir} {input}) &> {log}"
    # wrapper:
    #     "0.38.0/bio/fastqc"

# rule fastqcR2:
#     input:
#         "data_processing/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed.fastq" ##one for each R1 and one for R2 should be from a samples.yaml file
#     output:
#         html="qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.html",
#         zip="qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#     params: ""
#     log:
#         "logs/qc/fastqc/{sample}_{seqID}_R2_trimmed.log"
#     singularity:
#         config["singularitys"]["fastqc"]
#     wrapper:
#         "0.38.0/bio/fastqc"


rule bwa_mem:
    input:
        reads=[
            "data_processing/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed.fastq.gz",
            "data_processing/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed.fastq.gz",
        ],
    output:
        "data_processing/{sample}_{seqID}/{sample}_{seqID}.bam",
    log:
        "logs/map/bwa_mem/{sample}_{seqID}.log",
    params:  #-M
        index=config["reference"]["bwa"],
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
    threads: 10
    container:
        config["singularitys"]["bwa"]  #bwa 0.7.17, samtools 1.9, picard 2.20.11
    shell:
        "(bwa mem -t {threads} {params.extra} {params.index} {input.reads} | samtools sort -o {output} - ) &> {log}"


rule samtools_index:
    input:
        "data_processing/{sample}_{seqID}/{sample}_{seqID}.bam",
    output:
        "data_processing/{sample}_{seqID}/{sample}_{seqID}.bam.bai",
    log:
        "logs/map/samtools_index/{sample}_{seqID}.log",  # optional params string
    container:
        config["singularitys"]["bwa"]
    shell:
        "(samtools index {input} {output}) &> {log}"

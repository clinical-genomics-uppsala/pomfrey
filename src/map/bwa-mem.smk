
rule bwa_mem:
    input:
        reads = ["data_processing/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed.fastq.gz", "data_processing/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed.fastq.gz"]
    output:
        "data_processing/{sample}_{seqID}/{sample}_{seqID}.bam"
    log:
        "logs/map/bwa_mem/{sample}_{seqID}.log"
    params:   #-M
        index = config["reference"]["bwa"],
        extra = r"-R '@RG\tID:{sample}\tSM:{sample}'",
        # sort = "samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order = "coordinate"  # Can be 'queryname' or 'coordinate'.
        # sort_extra = ""            # Extra args for samtools/picard.
    threads: 8
    container:
        config["singularitys"]["bwa"]  #bwa 0.7.17, samtools 1.9, picard 2.20.11
    shell:
        "(bwa mem -t {threads} {params.extra} {params.index} {input.reads} | samtools sort -o {output} - ) &> {log}"
    # wrapper:
    #     "0.38.0/bio/bwa/mem"

rule samtools_index:
    input:
        "data_processing/{sample}_{seqID}/{sample}_{seqID}.bam"
    output:
        "data_processing/{sample}_{seqID}/{sample}_{seqID}.bam.bai"
    params:
        "" # optional params string
    log:
        "logs/map/samtools_index/{sample}_{seqID}.log",  # optional params string
    container:
        config["singularitys"]["bwa"]
    shell:
        "(samtools index {input} {output}) &> {log}"
    # wrapper:
    #     "0.38.0/bio/samtools/index"

def getFastqs(wildcards):
    fastq1 = config["samples"][wildcards.sample]
    fastq2 = fastq1.replace("_R1_", "_R2_")
    return [fastq1, fastq2]


rule cutadapt:
    input:
        getFastqs,
    output:
        fastq1="data_processing/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed.fastq.gz",
        fastq2="data_processing/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed.fastq.gz",
        qc="data_processing/{sample}_{seqID}/{sample}_{seqID}.qc.txt",
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters_r1="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",  #"-a AGAGCACACGTCTGAACTCCAGTCAC -g AGATCGGAAGAGCACACGT",
        adapters_r2="-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",  #"-A AGAGCACACGTCTGAACTCCAGTCAC -G AGATCGGAAGAGCACACGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        others="--minimum-length 2 -q 20",
    log:
        "logs/trim/cutadapt/{sample}_{seqID}.log",
    threads: 10
    singularity:
        config["singularitys"]["cutadapt"]
    shell:
        "(cutadapt {params.adapters_r1} {params.adapters_r2} {params.others} -o {output.fastq1} -p {output.fastq2} "
        "-j {threads} {input} > {output.qc} ) &> {log}"

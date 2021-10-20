localrules:
    createReadCountPanelOfNormals,


rule all:
    input:
        "Normals/GATK4/readCountPoN.hdf5",
        expand("Normals/GATK4/{normal}.counts.hdf5", normal=config["normal"]),


## Create interval list and preprocess interval list, only needed when updating bedfile
rule bedToIntervalList:
    input:
        bed=config["bed"]["bedfile"],  ## Annotated clostest bedfile? Really needed?
        refDict=config["reference"]["ref"],  ##Have to be a .dict in same folder as .fasta
    output:
        "bedFiles/TM_TE-annotated_closest-noduplicates.interval_list",  #Should be based on bedfile...
    log:
        "logs/Normals/TM_TE-annotated_closest-noduplicates.log",
    container:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk BedToIntervalList  -I {input.bed} -O {output} "
        "-SD {input.refDict} ) &> {log} "


rule preprocessIntervals:
    input:
        ref=config["reference"]["ref"],
        intervalList="bedFiles/TM_TE-annotated_closest-noduplicates.interval_list",  #targets_C.interval_list interval list picard style
    output:
        "bedFiles/TM_TE-annotated_closest-noduplicates.preprocessed.interval_list",
    params:
        binLength=0,  #WGS 1000
        mergingRule="OVERLAPPING_ONLY",
    log:
        "logs/Normals/GATK/preprocessIntervals.log",
    container:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' PreprocessIntervals -L {input.intervalList} -R {input.ref} "
        "--bin-length {params.binLength} --interval-merging-rule {params.mergingRule} "
        "-O {output} ) &> {log} "


# From here need to be redone when added new samples
rule collectReadCounts:
    input:
        bam=lambda wildcards: config["normal"][wildcards.normal],
        interval="bedFiles/TM_TE-annotated_closest-noduplicates.preprocessed.interval_list",
    output:
        "Normals/GATK4/{normal}.counts.hdf5",  #Should have date in it?
    params:
        mergingRule="OVERLAPPING_ONLY",
    log:
        "logs/Normals/GATK/{normal}.collectReadCounts.log",
    container:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' CollectReadCounts -I {input.bam} -L {input.interval} "
        "--interval-merging-rule {params.mergingRule} -O {output} ) &> {log}"


rule createReadCountPanelOfNormals:
    input:
        expand("Normals/GATK4/{normal}.counts.hdf5", normal=config["normal"]),
    output:
        "Normals/GATK4/readCountPoN.hdf5",
    params:
        minIntervalMedianPerc=5.0,
        input=lambda wildcards, input: " -I ".join(input),
    log:
        "logs/Normals/GATK/readCountPoN.log",
    container:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' CreateReadCountPanelOfNormals -I {params.input} "
        "--minimum-interval-median-percentile {params.minIntervalMedianPerc} "
        "-O {output} ) &> {log}"

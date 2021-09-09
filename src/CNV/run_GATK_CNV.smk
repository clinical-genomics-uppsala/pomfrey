# wildcard_constraints:
#     seqID = config["seqID"]["sequencerun"]


rule collectReadCounts:
    input:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",  #lambda wildcards: config["samples"][wildcards.sample],
        bai="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
        interval=config["CNV"]["interval"],  #"/projects/wp4/nobackup/workspace/arielle_test/CNV/bedFiles/TM_TE-annotated_closest-noduplicates.preprocessed.interval_list" #Better version? Should be same as other intervallist.
    output:
        "CNV/{sample}_{seqID}/{sample}_{seqID}.counts.hdf5",
    params:
        mergingRule="OVERLAPPING_ONLY",
    log:
        "logs/CNV/{sample}_{seqID}.collectReadCounts.log",
    singularity:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' CollectReadCounts -I {input.bam} -L {input.interval} "
        "--interval-merging-rule {params.mergingRule} -O {output} ) &> {log}"


rule denoiseReadCounts:
    input:
        hdf5PoN=config["CNV"]["PoN"],
        hdf5Tumor="CNV/{sample}_{seqID}/{sample}_{seqID}.counts.hdf5",
    output:
        stdCopyRatio="CNV/{sample}_{seqID}/{sample}_{seqID}_clean.standardizedCR.tsv",
        denoisedCopyRatio="CNV/{sample}_{seqID}/{sample}_{seqID}_clean.denoisedCR.tsv",
    log:
        "logs/CNV/{sample}_{seqID}-denoise.log",
    singularity:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' DenoiseReadCounts -I {input.hdf5Tumor} --count-panel-of-normals {input.hdf5PoN} "
        "--standardized-copy-ratios {output.stdCopyRatio} --denoised-copy-ratios {output.denoisedCopyRatio} ) &> {log}"


rule collectAllelicCounts:
    input:
        intervalList=config["CNV"]["interval"],  #Better version? Should be same as other intervallist. Also should be in config.
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",  #lambda wildcards: config["samples"][wildcards.sample],
        bai="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
        ref=config["reference"]["ref"],
    output:
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.allelicCounts.tsv",
    log:
        "logs/CNV/{sample}_{seqID}_allelicCounts.log",
    singularity:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' CollectAllelicCounts -L {input.intervalList} -I {input.bam} -R {input.ref} "
        "-O {output} ) &> {log}"


rule modelSegments:
    input:
        denoisedCopyRatio="CNV/{sample}_{seqID}/{sample}_{seqID}_clean.denoisedCR.tsv",
        allelicCounts="CNV/{sample}_{seqID}/{sample}_{seqID}_clean.allelicCounts.tsv",
    output:
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.modelBegin.seg",
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.modelFinal.seg",
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.cr.seg",
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.modelBegin.af.param",
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.modelBegin.cr.param",
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.modelFinal.af.param",
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.modelFinal.cr.param",
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.hets.tsv",
    params:
        outDir="CNV/{sample}_{seqID}/",
        outPrefix="{sample}_{seqID}_clean",
    log:
        "logs/CNV/{sample}_{seqID}_modelSegments.log",
    singularity:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' ModelSegments --denoised-copy-ratios {input.denoisedCopyRatio} "
        "--allelic-counts {input.allelicCounts} --output {params.outDir} --output-prefix {params.outPrefix} ) &> {log}"


rule callCopyRatioSegments:
    input:
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.cr.seg",
    output:
        "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.calledCNVs.seg",
    log:
        "logs/CNV/{sample}_{seqID}_calledCRSegments.log",
    singularity:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk CallCopyRatioSegments --input {input} --output {output} ) &> {log}"


rule plotModeledSegments:
    input:
        denoisedCopyRatio="CNV/{sample}_{seqID}/{sample}_{seqID}_clean.denoisedCR.tsv",
        allelicCounts="CNV/{sample}_{seqID}/{sample}_{seqID}_clean.hets.tsv",
        segments="CNV/{sample}_{seqID}/{sample}_{seqID}_clean.modelFinal.seg",  #Vad hander om man anv'nder Results/{sample}/Reports/{sample}_clean.calledCNVs.seg ist'llet?
        refDict=config["reference"]["ref"][:-5] + "dict",
    output:
        "CNV/{sample}_{seqID}_clean.calledCNVs.modeled.png",
    params:
        outDir="CNV/",
        outPrefix="{sample}_{seqID}_clean.calledCNVs",  #--minimum-contig-length 46709983
        pointSize=2.0,
    log:
        "logs/CNV/{sample}_{seqID}_plotSegments.log",
    singularity:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk PlotModeledSegments --denoised-copy-ratios {input.denoisedCopyRatio} --allelic-counts {input.allelicCounts} "
        "--segments {input.segments} --sequence-dictionary {input.refDict} --point-size-allele-fraction {params.pointSize} "
        "--point-size-copy-ratio {params.pointSize} --output {params.outDir} --output-prefix {params.outPrefix} ) &> {log} "


#
# rule Filter_cnv:
#     input:
#         gatkSeg = "CNV/{sample}_{seqID}/{sample}_{seqID}_clean.calledCNVs.seg",
#         bedfile = config["CNV"]["bedPoN"],
#         png = "CNV/{sample}_{seqID}_clean.calledCNVs.modeled.png", #"Results/{sample}_{seqID}/Reports/{sample}_{seqID}_clean.calledCNVs.modeled.png",
#         cytoCoord = config["CNV"]["cyto"]
#     output:
#         relevant_cnvs = "CNV/{sample}_{seqID}_clean.calledCNV-relevant_cnv-GATK4.xlsx" #"Results/{sample}_{seqID}/Reports/{sample}_{seqID}_clean.calledCNV-relevant_cnv-GATK4.xlsx"
#         # cnv_done = "Tumor/{sample}_cnv_done.txt"
#     params:
#         outdir = "CNV/"#"Results/{sample}_{seqID}/Reports",
#     container:
#         config["singularitys"]["python"]
#     log:
#         "logs/Tumor/{sample}_{seqID}_relevant_cnvs-gatk4.log"
#     shell:
#         "( python3 /projects/wp4/nobackup/workspace/arielle_test/CNV/bin/filter_gatkCalls.py {input.gatkSeg} {input.bedfile} {input.png} {params.outdir} {input.cytoCoord} ) &> {log}"

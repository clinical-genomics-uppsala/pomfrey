rule vcf2excel:
    input:
        vcf_snv="variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz",
        vcf_pindel="variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz",
        gatk_seg="CNV/{sample}_{seqID}/{sample}_{seqID}_clean.calledCNVs.seg",
        gatk_png="CNV/{sample}_{seqID}_clean.calledCNVs.modeled.png",
        picard_dup="qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt",
        mosdepth_summary="qc/mosdepth/{sample}_{seqID}.mosdepth.summary.txt",
        mosdepth_lowcov="qc/mosdepth/{sample}_{seqID}.mosdepth.lowCov.regions.txt",
        mosdepth_regions="qc/mosdepth/{sample}_{seqID}.regions.bed.gz",
        mosdepth_thresh_summary="qc/mosdepth/{sample}_{seqID}.thresholds_summary.txt",
        mosdepth_hotspot="qc/mosdepth/{sample}_{seqID}.mosdepth.hotspots.txt",
        # In configfile
        bedfile=config["bed"]["bedfile"],
        bedfile_pindel=config["bed"]["pindel"],
        bedfile_cnv=config["CNV"]["bedPoN"],
        cyto_coord_convert=config["CNV"]["cyto"],
        hotspot=config["bed"]["hotspot"],
        artefact_snv=config["bed"]["artefact"],
        artefact_pindel=config["bed"]["pindelArtefact"],
        germline=config["bed"]["germline"],
        hemato_count=config["configCache"]["hemato"],
        variantslog=config["configCache"]["variantlist"],
        igv_wait="Results/{sample}_{seqID}/Reports/IGV/done-igv.txt",
    output:
        "Results/{sample}_{seqID}/Reports/{sample}_{seqID}.xlsx",
    params:
        seqid=config["seqID"]["sequencerun"],
        thresholds=config["misc"]["cov_thresholds"],
        singularitys=config["singularitys"],
    log:
        "logs/report/{sample}_{seqID}.vcf2excel.log",
    wildcard_constraints:
        sample="(?!HD829).*",
    container:
        config["singularitys"]["python"]
    script:
        "vcf2excel.py"


rule vcf2excelHD829:
    input:
        snv="variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz",
        indel="variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz",
        cart="qc/{sample}_{seqID}/{sample}_{seqID}_MeanCoverageShortList.csv",
        # In configfile
        bed=config["bed"]["pindel"],
        hotspot=config["bed"]["hotspot"],
        artefact=config["bed"]["artefact"],
        germline=config["bed"]["germline"],
        hematoCount=config["configCache"]["hemato"],
        variantsLog=config["configCache"]["variantlist"],
        shortCov="qc/{sample}_{seqID}/{sample}_{seqID}_coverageShortHotspot.tsv",
    output:
        "Results/{sample}_{seqID}/Reports/{sample}_{seqID}.xlsx",
    params:
        configfile=config["seqID"]["sequencerun"] + "_config.yaml",
        dir=config["programdir"]["dir"],
    wildcard_constraints:
        sample="(HD829).*",
    log:
        "logs/report/{sample}_{seqID}.vcf2excel.log",
    container:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params.dir}/src/report/vcf2excelHD829.py {input.snv} {input.indel} {input.cart} {output} "
        "{params.configfile}) &> {log}"

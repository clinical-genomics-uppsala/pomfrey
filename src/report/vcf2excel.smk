localrules:
    fixCoverageHotspot,


rule fixCoverageHotspot:
    input:
        tsv="qc/{sample}_{seqID}/{sample}_{seqID}_coverage.tsv",
        bed=config["bed"]["hotspot"],
    output:
        "qc/{sample}_{seqID}/{sample}_{seqID}_coverageShortHotspot.tsv",
    log:
        "logs/report/{sample}_{seqID}.covShortHotspot.log",
    shell:
        """ ( while read line; do chr=$(echo $line | awk '{{print $1}}'); pos=$(echo $line | awk '{{print $2}}');
        cat {input.tsv} | grep ${{chr}} | grep ${{pos}} >>{output} ; done < {input.bed} ) &> {log} """


rule vcf2excel:
    input:
        snv="variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz",
        indel="variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz",
        gatkSeg="CNV/{sample}_{seqID}/{sample}_{seqID}_clean.calledCNVs.seg",
        png="CNV/{sample}_{seqID}_clean.calledCNVs.modeled.png",
        cart="qc/{sample}_{seqID}/{sample}_{seqID}_MeanCoverageShortList.csv",
        # In configfile
        bed=config["bed"]["pindel"],
        cnvbed=config["CNV"]["bedPoN"],
        cytoCoord=config["CNV"]["cyto"],
        hotspot=config["bed"]["hotspot"],
        artefact=config["bed"]["artefact"],
        pindelArtefact=config["bed"]["pindelArtefact"],
        germline=config["bed"]["germline"],
        hematoCount=config["configCache"]["hemato"],
        variantsLog=config["configCache"]["variantlist"],
        shortCov="qc/{sample}_{seqID}/{sample}_{seqID}_coverageShortHotspot.tsv",
        igv="Results/{sample}_{seqID}/Reports/IGV/done-igv.txt",
    output:
        "Results/{sample}_{seqID}/Reports/{sample}_{seqID}.xlsx",
    params:
        configfile=config["seqID"]["sequencerun"] + "_config.yaml",
        dir=config["programdir"]["dir"],
    log:
        "logs/report/{sample}_{seqID}.vcf2excel.log",
    wildcard_constraints:
        sample="(?!HD829).*",
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params.dir}/src/report/vcf2excel.py {input.snv} {input.indel} {input.gatkSeg} {input.png} {input.cart} "
        "{output} {params.configfile}) &> {log}"


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
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params.dir}/src/report/vcf2excelHD829.py {input.snv} {input.indel} {input.cart} {output} "
        "{params.configfile}) &> {log}"

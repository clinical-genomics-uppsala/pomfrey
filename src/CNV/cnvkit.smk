rule cnvkit_batch:
    input:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        bai="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
        ref=config["CNV"]["cnvkitpool"],
    output:
        regions="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.cnr",
        segments="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.cns",
        segments_called="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.call.cns",
        bins="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.bintest.cns",
        target_coverage="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.targetcoverage.cnn",
        antitarget_coverage="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.antitargetcoverage.cnn",
    container:
        config["singularitys"]["cnvkit"]
    shell:
        "cnvkit.py batch {input.bam} -r {input.ref} -d CNV/{wildcards.sample}_{wildcards.seqID}/cnvkit/"


rule cnvkit_diagram:
    input:
        cns="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.cns",
        cnr="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.cnr",
    output:
        pdf="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.diagram.pdf",
    container:
        config["singularitys"]["cnvkit"]
    shell:
        "cnvkit.py diagram -s {input.cns} {input.cnr} -o {output.pdf}"


rule cnvkit_call:
    input:
        cns="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.cns",
        vcf="Results/{sample}_{seqID}/Data/{sample}_{seqID}.normalized.genome.vcf.gz",
    output:
        "CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.loh.cns",
    container:
        config["singularitys"]["cnvkit"]
    shell:
        "cnvkit.py call {input.cns} -v {input.vcf} -o {output}"


rule cnvkit_scatter_loh:
    input:
        cns="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.loh.cns",
        cnr="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.cnr",
        vcf="Results/{sample}_{seqID}/Data/{sample}_{seqID}.normalized.genome.vcf.gz",
    output:
        "CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.loh.scatter.png",
    container:
        config["singularitys"]["cnvkit"]
    shell:
        "cnvkit.py scatter -s {input.cns} {input.cnr} -v {input.vcf} -o {output}  "


rule cnvkit_scatter_loh_perchr:
    input:
        cns="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.loh.cns",
        cnr="CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.cnr",
        vcf="Results/{sample}_{seqID}/Data/{sample}_{seqID}.normalized.genome.vcf.gz",
    output:
        "CNV/{sample}_{seqID}/cnvkit/{sample}_{seqID}-dedup.loh.scatter_{chr}.png",
    container:
        config["singularitys"]["cnvkit"]
    params:
        gene=lambda wildcards: config["CNA"][wildcards.chr],
    shell:
        """
        if [ -z {params.gene} ]
        then
        cnvkit.py scatter -s {input.cns} {input.cnr} -v {input.vcf} -c {wildcards.chr} -o {output}
        else
        cnvkit.py scatter -s {input.cns} {input.cnr} -v {input.vcf} -c {wildcards.chr} -o {output} -g {params.gene}
        fi
        """

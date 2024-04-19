rule mergeSNVPindel:
    input:
        snv="variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz",
        pindel="variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz",
    output:
        "Results/{sample}_{seqID}/Data/{sample}_{seqID}.SNV-pindel.vcf",
    params:
        "-a -D -O v",  #Allow overlap, Remove duplicates, output format vcf
    log:
        "logs/report/{sample}_{seqID}.mergeVcf.log",
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools concat {params} -o {output} {input.snv} {input.pindel} ) &> {log}"


rule makePassVCF:
    input:
        vcf_snv="variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz",
        vcf_indel="variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz",
        artefact=config["bed"]["artefact"],
        germline=config["bed"]["germline"],
        hemato=config["configCache"]["hemato"],
    output:
        vcf=temp("Results/{sample}_{seqID}/Reports/{sample}_{seqID}.PASS.vcf"),
    log:
        "logs/report/{sample}_{seqID}.PASS.vcf.log",
    wildcard_constraints:
        sample="(?!HD829).*",
    container:
        config["singularitys"]["python"]
    script:
        "makePASSvcf.py"


rule appendPindeltoPASS:
    input:
        pindel="variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz",
        PASS="Results/{sample}_{seqID}/Reports/{sample}_{seqID}.PASS.vcf",
    output:
        temp("Results/{sample}_{seqID}/Reports/{sample}_{seqID}.pindel.done"),
    log:
        "logs/report/{sample}_{seqID}.pindel.log",
    wildcard_constraints:
        sample="(?!HD829).*",
    container:
        config["singularitys"]["python"]
    shell:
        "(zcat {input.pindel} | grep -v '^#' | grep PASS >> {input.PASS} || true && touch {output} ) &> {log}"


rule createBatFile:
    input:
        vcf="Results/{sample}_{seqID}/Reports/{sample}_{seqID}.PASS.vcf",
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        pindelWait="Results/{sample}_{seqID}/Reports/{sample}_{seqID}.pindel.done",
        bed=config["bed"]["bedfile"],
        ref=config["reference"]["ref"],  #until build own .genome file.
    output:
        bat="Results/{sample}_{seqID}/Reports/IGV/{sample}_{seqID}-igv.bat",
    params:
        outfolder="Results/{sample}_{seqID}/Reports/IGV/",
        padding="40",
        sort="base",  #Type of sorting: base, position, strand, quality, sample or readgroup.
        view="squish",  #Type of view, collaps, squished...
        format="svg",  #svg, jpg
    log:
        "logs/report/{sample}_{seqID}-makeBat.log",
    wildcard_constraints:
        sample="(?!HD829).*",
    container:
        config["singularitys"]["python"]
    script:
        "makeBatfile.py"


rule igv:
    input:
        bat="Results/{sample}_{seqID}/Reports/IGV/{sample}_{seqID}-igv.bat",
    output:
        touch("Results/{sample}_{seqID}/Reports/IGV/done-igv.txt"),
    log:
        "logs/report/{sample}_{seqID}.igv.log",
    wildcard_constraints:
        sample="(?!HD829).*",
    threads: 2
    container:
        config["singularitys"]["igv"]
    shell:
        "(xvfb-run --server-args='-screen 0 3200x2400x24' --auto-servernum --server-num=1 igv.sh -b {input.bat} ) &> {log}"

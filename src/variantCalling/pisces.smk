localrules:
    piscesFix,
    sortPisces,
    gVCFfinalIndex,
    renameSample,


from os.path import dirname


rule pisces:
    input:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        reffolder=dirname(config["reference"]["ref"]),
        index="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
    output:
        vcf=temp("variantCalls/callers/pisces/{sample}_{seqID}/{sample}_{seqID}-dedup.genome.vcf"),
    params:
        outfolder="variantCalls/callers/pisces/{sample}_{seqID}",
        bed=config["bed"]["bedfile"],
    threads: 1
    log:
        "logs/variantCalling/pisces/{sample}_{seqID}.log",
    singularity:
        config["singularitys"]["pisces"]
    shell:  #Remove gVCF False for genome vcf and save for db, and artifacts? -gVCF FALSE
        "(dotnet /app/Pisces/Pisces.dll -b {input.bam} -g {input.reffolder} -i {params.bed} -t {threads} "
        "--filterduplicates TRUE --outfolder {params.outfolder} ) &> {log}"


rule piscesFix:  ## use bcftools view --minalleles 2 {input} instead?
    input:
        "variantCalls/callers/pisces/{sample}_{seqID}/{sample}_{seqID}-dedup.genome.vcf",
    output:
        temp("variantCalls/callers/pisces/{sample}_{seqID}/{sample}_{seqID}.pisces.unsorted.vcf"),
    log:
        "logs/variantCalling/pisces/{sample}_{seqID}.2.log",
    shell:
        """( awk '{{if($5 != "." || $1 ~ /^"#"/)print $0}}' {input} >{output} ) 2> {log}"""


rule renameSample:
    input:
        "variantCalls/callers/pisces/{sample}_{seqID}/{sample}_{seqID}-dedup.genome.vcf",
    output:
        temp("variantCalls/callers/pisces/{sample}_{seqID}-name.txt"),
    log:
        "logs/variantCalling/pisces/{sample}_{seqID}-name.log",
    shell:
        "echo {wildcards.sample} > {output}"


rule sortPisces:
    input:
        vcf="variantCalls/callers/pisces/{sample}_{seqID}/{sample}_{seqID}.pisces.unsorted.vcf",
        name="variantCalls/callers/pisces/{sample}_{seqID}-name.txt",
    output:
        temp("variantCalls/callers/pisces/{sample}_{seqID}.pisces.vcf"),
    singularity:
        config["singularitys"]["bcftools"]
    log:
        "logs/variantCalling/pisces/{sample}_{seqID}.sort.log",
    shell:
        "(bcftools reheader -s {input.name} {input.vcf} | bcftools sort -o {output} -O v - ) &> {log}"


rule gVCFdecompose:
    input:
        vcf="variantCalls/callers/pisces/{sample}_{seqID}/{sample}_{seqID}-dedup.genome.vcf",
    output:
        temp("variantCalls/callers/pisces/{sample}_{seqID}/{sample}_{seqID}.decomp.genome.vcf"),
    log:
        "logs/variantCalling/pisces/{sample}_{seqID}.genome.decomp.log",
    singularity:
        config["singularitys"]["vt"]
    shell:
        "(vt decompose -s {input.vcf} | vt decompose_blocksub -o {output} -) &> {log}"


rule gVCFnormalize:
    input:
        vcf="variantCalls/callers/pisces/{sample}_{seqID}/{sample}_{seqID}.decomp.genome.vcf",
        fasta=config["reference"]["ref"],
    output:
        "Results/{sample}_{seqID}/Data/{sample}_{seqID}.normalized.genome.vcf.gz",
    log:
        "logs/variantCalling/pisces/{sample}_{seqID}.normalized.gVCF.log",
    singularity:
        config["singularitys"]["vt"]
    shell:
        "(vt normalize -n -r {input.fasta} -o {output} {input.vcf} ) &> {log}"


rule gVCFfinalIndex:
    input:
        vcf="Results/{sample}_{seqID}/Data/{sample}_{seqID}.normalized.genome.vcf.gz",
    output:
        "Results/{sample}_{seqID}/Data/{sample}_{seqID}.normalized.genome.vcf.gz.tbi",
    singularity:
        config["singularitys"]["bcftools"]
    log:
        "logs/variantCalling/pisces/{sample}_{seqID}.gz.log",
    shell:
        "(tabix {input.vcf}) 2> {log}"

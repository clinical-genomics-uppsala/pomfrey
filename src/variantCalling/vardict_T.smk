
rule vardict:
    input:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        index="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
        ref=config["reference"]["ref"],
        bed=config["bed"]["bedfile"],
    output:
        temp("variantCalls/callers/vardict/{sample}_{seqID}.vardict.vcf"),
    params:
        af="0.01",
    log:
        "logs/variantCalling/vardict/{sample}_{seqID}.log",
    threads: 4
    singularity:
        config["singularitys"]["vardict"]
    shell:
        "(vardict-java -G {input.ref} -f {params.af} -I 200 -th {threads} -N '{wildcards.sample}' -z -c 1 -S 2 -E 3 "
        "-b {input.bam} {input.bed} | teststrandbias.R | var2vcf_valid.pl -N '{wildcards.sample}' -E "
        "-f {params.af} > {output}) &> {log}"

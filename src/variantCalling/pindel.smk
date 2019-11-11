rule pindelConf: ##Add in excel file what genes were used.
    input:
        bam = "mapped/{sample}.bam",
        bai = "mapped/{sample}.bam.bai"
    output:
        "variantCalls/pindel/{sample}/{sample}-config.txt"
    log:
        "logs/pindel/{sample}.config.log"
    shell:
        "( echo -e '{input.bam}\t250\t{wildcards.sample}'>{output} ) &> {log}"

rule pindel:
    input:
        bed = lambda wildcards: config["bed"]["pindel"],
        ref = "/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta",
        bamconfig = "variantCalls/pindel/{sample}/{sample}-config.txt" #path to bam \t insert size \t sample name
    output:
        "variantCalls/pindel/{sample}/{sample}_BP",
        "variantCalls/pindel/{sample}/{sample}_CloseEndMapped",
        "variantCalls/pindel/{sample}/{sample}_D",
        "variantCalls/pindel/{sample}/{sample}_INT_final",
        "variantCalls/pindel/{sample}/{sample}_INV",
        "variantCalls/pindel/{sample}/{sample}_LI",
        "variantCalls/pindel/{sample}/{sample}_RP",
        "variantCalls/pindel/{sample}/{sample}_SI",
        "variantCalls/pindel/{sample}/{sample}_TD"
    params:
        x = 2,
        B = 60
    log:
        "logs/pindel/{sample}.pindel.log"
    singularity:
        "pindel-0.2.5b8.simg"
    threads:    4
    shell:
        " (pindel -f {input.ref} -i {input.bamconfig} -T {threads} -x {params.x} -B {params.B} -j {input.bed} -o variantCalls/pindel/{wildcards.sample}/{wildcards.sample} ) &> {log}"

rule pindel2vcf:
    input:
        ref = "/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta",
        bp = "variantCalls/pindel/{sample}/{sample}_BP",
        closeend = "variantCalls/pindel/{sample}/{sample}_CloseEndMapped",
        d = "variantCalls/pindel/{sample}/{sample}_D",
        final = "variantCalls/pindel/{sample}/{sample}_INT_final",
        inv = "variantCalls/pindel/{sample}/{sample}_INV",
        li = "variantCalls/pindel/{sample}/{sample}_LI",
        rp = "variantCalls/pindel/{sample}/{sample}_RP",
        si = "variantCalls/pindel/{sample}/{sample}_SI",
        td = "variantCalls/pindel/{sample}/{sample}_TD"
    output:
        "variantCalls/pindel/{sample}.pindel.vcf"
    params:
        e = 35, #min supporting reads
        mc = 10, #min coverage
        minsize = 5, #min size of reported
        refname = "hg19",
        refdate = 000000
    log:
        "logs/pindel/{sample}.pindel2vcf.log"
    singularity:
        "pindel-0.2.5b8.simg"
    threads:    1
    shell:
        "(pindel2vcf -P variantCalls/pindel/{wildcards.sample}/{wildcards.sample} -r {input.ref} -R {params.refname} -d {params.refdate} -v {output} -e {params.e} -mc {params.mc} -G -is {params.minsize} ) &> {log}"

rule pindelIndex:
    input:
        "variantCalls/pindel/{sample}.pindel.vcf"
    output:
        gz = "variantCalls/pindel/{sample}.pindel.vcf.gz",
        tbi = "variantCalls/pindel/{sample}.pindel.vcf.gz.tbi"
    log:
        "logs/pindel/{sample}.index.log"
    singularity:
        "bcftools-1.9--8.simg"
    shell:
        "( bgzip {input} && tabix {input}.gz ) &> {log}"

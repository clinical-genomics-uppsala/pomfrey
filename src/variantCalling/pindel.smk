localrules:
    fixContigPindel,
    pindelConf,
    fixPindelDPoAF,
    filterPindel,
    bgzipPindel,


rule pindelConf:
    input:
        bam="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        bai="Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai",
    output:
        "variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}-config.txt",
    log:
        "logs/variantCalling/pindel/{sample}_{seqID}.config.log",
    shell:
        "( echo -e '{input.bam}\t250\t{wildcards.sample}'>{output} ) &> {log}"


rule pindel:
    input:
        bed=config["bed"]["pindel"],
        ref=config["reference"]["ref"],
        bamconfig="variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}-config.txt",
    output:
        "variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_BP",
        "variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_CloseEndMapped",
        "variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_D",
        "variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_INT_final",
        "variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_INV",
        "variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_LI",
        "variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_RP",
        "variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_SI",
        "variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_TD",
    params:
        x=2,
        B=60,
    log:
        "logs/variantCalling/pindel/{sample}_{seqID}.pindel.log",
    singularity:
        config["singularitys"]["pindel"]
    threads: 5
    shell:
        " (pindel -f {input.ref} -i {input.bamconfig} -T {threads} -x {params.x} -B {params.B} -j {input.bed} "
        "-o variantCalls/pindel/{wildcards.sample}_{wildcards.seqID}/{wildcards.sample}_{wildcards.seqID} ) &> {log}"


rule pindel2vcf:
    input:
        ref=config["reference"]["ref"],
        bp="variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_BP",
        closeend="variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_CloseEndMapped",
        d="variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_D",
        final="variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_INT_final",
        inv="variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_INV",
        li="variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_LI",
        rp="variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_RP",
        si="variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_SI",
        td="variantCalls/pindel/{sample}_{seqID}/{sample}_{seqID}_TD",
    output:
        temp("variantCalls/pindel/{sample}_{seqID}.pindel.noDP.noContig.vcf"),
    params:
        e=10,  # min supporting reads 35
        mc=10,  # min coverage
        minsize=5,  # min size of reported 5
        refname="hg19",
        refdate=000000,  # Can I add seqID instead? config["seqID"]["sequencerun"]
        he=0.01,  # Hetrozygot call to be included in QCI
    log:
        "logs/variantCalling/pindel/{sample}_{seqID}.pindel2vcf.log",
    singularity:
        config["singularitys"]["pindel"]
    threads: 1
    shell:
        "(pindel2vcf -P variantCalls/pindel/{wildcards.sample}_{wildcards.seqID}/{wildcards.sample}_{wildcards.seqID} "
        "-r {input.ref} -R {params.refname} -d {params.refdate} -v {output} -he {params.he} -e {params.e} -mc {params.mc} "
        "-G -is {params.minsize} ) &> {log}"


rule fixContigPindel:
    input:
        vcf="variantCalls/pindel/{sample}_{seqID}.pindel.noDP.noContig.vcf",
        fasta=config["reference"]["ref"],
    output:
        temp("variantCalls/pindel/{sample}_{seqID}.pindel.noDP.vcf"),
    log:
        "logs/variantCalling/pindel/{sample}_{seqID}.fixContig.log",
    singularity:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk UpdateVcfSequenceDictionary -I {input.vcf} -SD {input.fasta} -O {output}) &> {log}"


rule fixPindelDPoAF:
    input:
        "variantCalls/pindel/{sample}_{seqID}.pindel.noDP.vcf",
    output:
        "variantCalls/pindel/{sample}_{seqID}.pindel.vcf",
    params:
        config["programdir"]["dir"],
    log:
        "logs/variantCalling/{sample}_{seqID}.fixDP.log",
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params}/src/variantCalling/fix_pindelDPoAF.py {input} {output}) &> {log}"


rule annotatePindel:
    input:
        vcf="variantCalls/pindel/{sample}_{seqID}.pindel.vcf",
        fasta=config["reference"]["ref"],
        cache=config["configCache"]["vep"],
    output:
        temp("variantCalls/pindel/{sample}_{seqID}.pindel.ann.vcf"),
    params:
        "--check_existing --pick --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory "
        "--canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --af --af_1kg --af_gnomad --max_af "
        "--pubmed --variant_class ",
    log:
        "logs/variantCalling/pindel/{sample}_{seqID}.ann.log",
    threads: 8
    singularity:
        config["singularitys"]["vep"]
    shell:
        "(if [[ $(cat {input.vcf} | grep -v '^#' | wc -l) -eq 0 ]]; then mv {input.vcf} {output}; else "
        "vep --vcf --no_stats -o {output} -i {input.vcf} --dir_cache {input.cache} --fork {threads} --cache --refseq "
        "--offline --fasta {input.fasta} {params} ; fi) &> {log}"


rule filterPindel:
    input:
        vcf="variantCalls/pindel/{sample}_{seqID}.pindel.ann.vcf",
    output:
        temp("variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf"),
    params:
        config["programdir"]["dir"],
    log:
        "logs/variantCalling/pindel.{sample}_{seqID}.filt.log",
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params}/src/variantCalling/filter_vcf.py {input.vcf} {output}) &> {log}"


rule bgzipPindel:
    input:
        "variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf",
    output:
        "variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz",
        "variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz.tbi",
    log:
        "logs/variantCalling/pindel/{sample}_{seqID}.bgzip-index.log",
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "( bgzip {input} && tabix {input}.gz ) &> {log}"

rule vep:
    input:
        vcf="variantCalls/recall/{sample}_{seqID}.vcf.gz",
        cache=config["configCache"]["vep"],
        fasta=config["reference"]["ref"],
    output:
        vcf=temp("variantCalls/annotation/raw/{sample}_{seqID}.raw.vcf"),
    params:
        "--check_existing --pick --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory "
        "--canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --af --af_1kg --af_gnomad --max_af "
        "--pubmed --variant_class ",
    log:
        "logs/variantCalling/vep/{sample}_{seqID}.log",
    container:
        config["singularitys"]["vep"]
    threads: 10
    shell:
        "(vep --vcf --no_stats -o {output.vcf} -i {input.vcf} --dir_cache {input.cache} --fork {threads} --cache "
        "--refseq --offline --fasta {input.fasta} {params} ) &> {log}"


rule bgzipVep:
    input:
        "variantCalls/annotation/raw/{sample}_{seqID}.raw.vcf",
    output:
        "variantCalls/annotation/raw/{sample}_{seqID}.raw.vcf.gz",
    log:
        "logs/variantCalling/vep/{sample}_{seqID}.bgzip.log",
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input}) &> {log}"


rule filterVep:
    input:
        vcf="variantCalls/annotation/raw/{sample}_{seqID}.raw.vcf.gz",
    output:
        vcf=temp("variantCalls/annotation/{sample}_{seqID}.filt.vcf"),
    log:
        "logs/variantCalling/vep/filter/{sample}_{seqID}.log",
    container:
        config["singularitys"]["python"]
    script:
        "filter_vcf.py"


rule bgzipSNV:
    input:
        "variantCalls/annotation/{sample}_{seqID}.filt.vcf",
    output:
        "variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz",
        "variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz.tbi",
    log:
        "logs/variantCalling/{sample}_{seqID}.bgzip.log",
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input} && tabix {input}.gz) &> {log}"

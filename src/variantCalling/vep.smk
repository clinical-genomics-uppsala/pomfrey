localrules: bgzipVep, filterVep, bgzipSNV
rule vep:
    input:
        vcf = "variantCalls/recall/{sample}_{seqID}.vcf.gz",
        cache = config["configCache"]["vep"], #"/opt/vep/.vep", ## always remeber the --bind vep-data:/opt/vep/.vep command in singularity args
        fasta = config["reference"]["ref"]
    output:
        vcf = temp("variantCalls/annotation/raw/{sample}_{seqID}.raw.vcf")
    params:
        "--check_existing --pick --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --af --af_1kg --af_gnomad --max_af --pubmed --variant_class "
        # "--everything --check_existing --pick"  #--exclude_null_alleles
    log:
        "logs/variantCalling/vep/{sample}_{seqID}.log"
    container:
        config["singularitys"]["vep"]
    threads:    8
    shell:
        "(vep --vcf --no_stats -o {output.vcf} -i {input.vcf} --dir_cache {input.cache} --fork {threads} --cache --refseq --offline --fasta {input.fasta} {params} ) &> {log}"

rule bgzipVep:
    input:
        "variantCalls/annotation/raw/{sample}_{seqID}.raw.vcf"
    output:
        "variantCalls/annotation/raw/{sample}_{seqID}.raw.vcf.gz"
    log:
        "logs/variantCalling/vep/{sample}_{seqID}.bgzip.log"
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input}) &> {log}"

rule filterVep:
    input:
        vcf="variantCalls/annotation/raw/{sample}_{seqID}.raw.vcf.gz"
    output:
        temp("variantCalls/annotation/{sample}_{seqID}.filt.vcf")
    params:
        config["programdir"]["dir"]
    log:
        "logs/variantCalling/vep/filter/{sample}_{seqID}.log"
    container:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params}/src/variantCalling/filter_vcf.py {input.vcf} {output}) &> {log}"

rule bgzipSNV:
    input:
        "variantCalls/annotation/{sample}_{seqID}.filt.vcf"
    output:
        "variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz",
        "variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz.tbi"
    log:
        "logs/variantCalling/{sample}_{seqID}.bgzip.log"
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input} && tabix {input}.gz) &> {log}"

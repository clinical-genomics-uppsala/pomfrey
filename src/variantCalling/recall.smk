localrules: sort_recall,createMultiVcf

rule recall:
    input:
        vcfs = expand("variantCalls/callers/{method}/{{sample}}_{{seqID}}.{method}.normalized.vcf.gz",  method=config["methods"]) ,  # same order as methods in config!! make sure that is correct
        tabix = expand("variantCalls/callers/{method}/{{sample}}_{{seqID}}.{method}.normalized.vcf.gz.tbi",  method=config["methods"]),
        ref = config["reference"]["ref"]
    output:
        vcf = "variantCalls/recall/{sample}_{seqID}.unsorted.vcf.gz"
    params:
        support =  "1", #"{support}" ,
        order = ",".join([ s for s in config["methods"]]) #"Vardict,Mutect2,Pisces,Freebayes" #,Manta" #Make sure that the order is correct! Order of methods in configfile
    log:
        "logs/variantCalling/recall/{sample}_{seqID}.log"
    container:
        config["singularitys"]["recall"]
        # "/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle/somaticpipeline/src/singularity/bcbio-variation-recall-0.2.6-0.simg" Fungerar inte beh;ver bcftool!!
        # "/gluster-storage-volume/projects/wp4/nobackup/workspace/somatic_dev/bcbio-variation-recall.simg" #Dev
    shell: ##Remove filtered?? if so --nofiltered
        "(bcbio-variation-recall ensemble -n {params.support} --names {params.order} {output.vcf} {input.ref} {input.vcfs}) &> {log}"

rule sort_recall:
    input:
        "variantCalls/recall/{sample}_{seqID}.unsorted.vcf.gz" #multiAllelic.vcf"
    output:
        tbiUnSorted = temp("variantCalls/recall/{sample}_{seqID}.unsorted.vcf.gz.tbi"),
        vcf = temp("variantCalls/recall/{sample}_{seqID}.notMulti.all.vcf.gz"),
        tbi = temp("variantCalls/recall/{sample}_{seqID}.notMulti.all.vcf.gz.tbi")
    log:
        "logs/variantCalling/recall/{sample}_{seqID}.sort.log"
    container:
        config["singularitys"]["bcftools"]
    shell:
        "( tabix -f {input} && \
        bcftools sort -o {output.vcf} -O z {input} && \
        tabix {output.vcf} ) &> {log}"

##Nya regeler
rule filter_recall:
    input:
        "variantCalls/recall/{sample}_{seqID}.notMulti.all.vcf.gz"
    output:
        "variantCalls/recall/{sample}_{seqID}.notMulti.vcf.gz"
    params:
        dir = config["programdir"]["dir"],
        indelArte = config["bed"]["indelartefact"]
    log:
        "logs/variantCalling/recall/{sample}_{seqID}.filter_recall.log"
    container:
        config["singularitys"]["python"]
    shell:
        "(python3 {params.dir}/src/variantCalling/filter_recall.py {input} {output} {params.indelArte}) &> {log}"

rule index_filterRecall:
    input:
        "variantCalls/recall/{sample}_{seqID}.notMulti.vcf.gz"
    output:
        tbi = "variantCalls/recall/{sample}_{seqID}.notMulti.vcf.gz.tbi"
    log:
        "logs/variantCalling/recall/{sample}_{seqID}.index_recallFilter.log"
    container:
        config["singularitys"]["bcftools"]
    shell:
        "( tabix {input} ) &> {log}"

# ##Add in multiallelic Variants
rule createMultiVcf: #Behovs denna?? Eller ar den onodig nu?
     input:
         "variantCalls/callers/pisces/{sample}_{seqID}.pisces.normalized.vcf.gz" ##based on dubbletter i genomeVCF fran pisces!!!
     output:
         "variantCalls/recall/{sample}_{seqID}.multiPASS.vcf"
     log:
         "logs/recall/{sample}_{seqID}.multiPASS.log"
     shell:
         """(zcat {input} | grep '^#' >{output} &&
         for pos in $(zcat {input} | grep -v '^#' | awk '{{print $2}}' | sort |uniq -d);do
         zcat {input}|grep $pos | grep PASS  >> {output} || true ; done) &> {log}"""

rule sort_multiPASS:
    input:
        "variantCalls/recall/{sample}_{seqID}.multiPASS.vcf"
    output:
        vcf = "variantCalls/recall/{sample}_{seqID}.multiPASS.sort.vcf.gz", #temp
        tbi = "variantCalls/recall/{sample}_{seqID}.multiPASS.sort.vcf.gz.tbi" #temp
    log:
        "logs/recall/{sample}_{seqID}.multiPASS.sort.log"
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools sort -o {output.vcf} -O z {input} && tabix {output.vcf}) &> {log}"

rule concatMulti:
    input:
        vcf = "variantCalls/recall/{sample}_{seqID}.notMulti.vcf.gz",
        tbi = "variantCalls/recall/{sample}_{seqID}.notMulti.vcf.gz.tbi",
        multi = "variantCalls/recall/{sample}_{seqID}.multiPASS.sort.vcf.gz"
    output:
        "variantCalls/recall/{sample}_{seqID}.vcf.gz"
    params:
        "--allow-overlaps -d all -O z"
    log:
        "logs/recall/{sample}_{seqID}.concat.log"
    container:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools concat {params} -o {output} {input.vcf} {input.multi}) &> {log}"

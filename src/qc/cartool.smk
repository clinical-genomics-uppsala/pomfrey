# localrules: fixoutput

rule cartool:
    input:
        bam = "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        bed = config["bed"]["cartool"],
        bai = "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai"
    output:
        statstable = "qc/{sample}_{seqID}/{sample}_{seqID}_Stat_table.csv",
        cartoollog =  "qc/{sample}_{seqID}/{sample}_{seqID}_Log.csv",
        coverage = "qc/{sample}_{seqID}/{sample}_{seqID}_coverage.tsv",
        full = "qc/{sample}_{seqID}/{sample}_{seqID}_MeanCoverageFullList.csv",
        short = "qc/{sample}_{seqID}/{sample}_{seqID}_MeanCoverageShortList.csv"
    params:
        user = "arielle",
        coverage = config["cartool"]["cov"],
        extra = "-k" #k = combine, p= mapQ
    log:
        "logs/qc/CARTool/{sample}_{seqID}.cartool.log",
    container:
        config["singularitys"]["cartool"]
    shell: #Need to fix -o so no space is needed.
        "( python3.6 /opt/CARtool/ProgramLancher.py -a {input.bed} -b {input.bam} -c {params.coverage} -e {params.user} -o qc/{wildcards.sample}_{wildcards.seqID}/ {wildcards.sample}_{wildcards.seqID} {params.extra} )&> {log}"

#
# rule fixoutput:
#     input:
#         header = "/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/report/multiqc-header.txt", #Always change to correct full path
#         stat = "qc/{sample}/{sample}_Stat_table.csv"
#     output:
#         "qc/{sample}/{sample}_cartool_mqc.csv"
#     log:
#         "logs/qc/CARTool/{sample}.fix.log"
#     shell:
#         """ (cat {input.header} > {output} && cut -d ',' -f2- {input.stat} | tr -d "\15\32" >>{output} ) &> {log}"""

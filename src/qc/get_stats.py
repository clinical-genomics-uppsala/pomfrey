#!/bin/python3
import sys
import subprocess
import csv
import json
from collections import OrderedDict


def extractMatchingLines(expressionMatch, artefactFile, grepvarible):
    if grepvarible == '':
        grepvarible = '-wE '
    cmdArt = 'grep '+grepvarible+' '+str(expressionMatch)+' '+artefactFile
    matchLines = subprocess.run(cmdArt, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
    return matchLines


min_cov = int(snakemake.params.thresholds.split(',')[0])
med_cov = int(snakemake.params.thresholds.split(',')[1])
max_cov = int(snakemake.params.thresholds.split(',')[2])

data_json = {}
for picard_dup in snakemake.input.picard_dup:
    # percent duplicateLevel
    duplicateLevel = extractMatchingLines('PERCENT', picard_dup, '-A1 ').split('\n')[-1].split('\t')[8]  # need *100 to be percent
    duplicates_percent = str(round(float(duplicateLevel)*100, 2))
    sample = picard_dup.split("/")[-1].split("_")[0]
    data_json[str(sample)] = {}
    data_json[str(sample)]['Duplicates [%]'] = duplicates_percent

for picard_metrics in snakemake.input.picard_metrics:
    met = extractMatchingLines("BAIT_SET", picard_metrics, "-A1 ")
    metrics = met.split('\n')
    zipObject = zip(metrics[0].split('\t'), metrics[1].split('\t'))
    metricsDict = dict(zipObject)
    sample = picard_metrics.split("/")[-1].split("_")[0]
    data_json[str(sample)]['Bases on Target'] = metricsDict['PCT_SELECTED_BASES']
    data_json[str(sample)]['Breadth 100x'] = metricsDict['PCT_TARGET_BASES_100X']
    data_json[str(sample)]['Breadth 50x'] = metricsDict['PCT_TARGET_BASES_50X']
    # metricsDict['PCT_SELECTED_BASES'] ##bases on target
    # metricsDict['PCT_TARGET_BASES_50X'] #
    # metricsDict['PCT_TARGET_BASES_100X']

for samtools_stats in snakemake.input.samtools_stats:
    sam = extractMatchingLines("SN", samtools_stats, "")
    sam = sam.split('\nSN\t')
    listOfList = [i.split('\t') for i in sam[1:]]
    samDict = {item[0].strip(':'): item[1] for item in listOfList}
    sample = samtools_stats.split("/")[-1].split("_")[0]
    data_json[str(sample)]['Average Quality'] = samDict['average quality']
    data_json[str(sample)]['Insert Size s.d.'] = samDict['insert size standard deviation']
    data_json[str(sample)]['Insert Size'] = samDict['insert size average']
    data_json[str(sample)]['Reads Paired [%]'] = samDict['percentage of properly paired reads (%)']
    data_json[str(sample)]['Reads Mapped'] = samDict['reads mapped']
    data_json[str(sample)]['Tot Seq'] = samDict['raw total sequences']
    # samDict['raw total sequences']
    # samDict['reads mapped']
    # samDict['percentage of properly paired reads (%)']
    # samDict['insert size average']
    # samDict['insert size standard deviation']
    # samDict['average quality']

for mosdepth_summary in snakemake.input.mosdepth_summary:
    avg_cov = extractMatchingLines("total_region", mosdepth_summary, "").split('\t')[3]
    sample = mosdepth_summary.split("/")[-1].split("_")[0]
    data_json[str(sample)]['Avg Coverage'] = avg_cov

for mosdepth_thresh_summary in snakemake.input.mosdepth_thresh_summary:
    with open(mosdepth_thresh_summary, 'r') as threshold_summary:
        thresholds = threshold_summary.read().strip().split("\t")
    sample = mosdepth_thresh_summary.split("/")[-1].split("_")[0]
    data_json[str(sample)]['Breadth '+str(med_cov)+'x'] = thresholds[1]

startReading = 0
samples = []
with open(snakemake.input.samplesheet, 'r') as samplesheet:
    lines = [line.strip() for line in samplesheet]
    for line in lines:
        if startReading == 1:  # Once reached [Data]
            samples.append(line.split(',')[1])
        elif line.startswith("[Data]"):
            startReading = 1
samples = samples[1:]  # Remove header from SampleSheet
samplesheet_order = [string for string in samples if string != ""]  # Remove empty fields
data_json_ordered = {}
data_json_ordered["data"] = dict(sorted(data_json.items(), key=lambda sample: samplesheet_order.index(sample[0])))


header_json = {"headers": {
                    "Avg Coverage": {
                        "title": "Average Coverage",
                        "description": "Avg cov of bedfile from Mosdepth",
                    },
                    "Average Quality": {
                        "title": "Average Quality",
                        "description": "Average mapping quality from Samtools",
                        "min": 0,
                        "max": 60,
                        "scale": "RdYlGn",
                    },
                    "Bases on Target": {
                        "title": "Bases on Target",
                        "description": "Ratio of bases on target from Picard HsMetrics",
                        "format": "{:,.3f}",
                    },
                    "Duplicates [%]": {
                        "title": "Duplicates [%]",
                        "description": "Percent duplicates marked by Picard MarkDuplicates",
                        "min": 0,
                        "max": 50,
                        "scale": "RdYlGn-rev",
                        "suffix": "%",
                    },
                    "Breadth " + str(med_cov) + "x": {
                        "title": "Coverage breadth " + str(med_cov) + "x",
                        "description": "Design covered to " + str(med_cov) + "x from Mosdepth",
                        "min": 0,
                        "max": 1,
                        "scale": "RdYlGn",
                        "format": "{:,.2f}",
                    },
                    "Breadth 100x": {
                        "title": "Coverage Breadth 100x",
                        "description": "Design covered to 100x or over from Picard HsMetrics",
                        "format": "{:,.2f}",
                    },
                    "Breadth 50x": {
                        "title": "Coverage Breadth 50x",
                        "description": "Design covered to 50x or over from Picard HsMetrics",
                        "format": "{:,.2f}",
                    },
                    "Tot Seq": {
                        "title": "Total Sequences",
                        "description": "Number of reads in fastq from Samtools stats",
                        "format": "{:,.0f}",
                    },
                    "Reads Mapped": {
                        "title": "Reads Mapped",
                        "description": "Number of reads mapped from Samtools stats",
                        "format": "{:,.0f}",
                    },
                    "Reads Paired [%]": {
                        "title": "Reads Properly Paired",
                        "description": "Percent properly paired reads from Samtools stats",
                        "min": 0,
                        "max": 100,
                        "scale": "RdYlGn",
                        "suffix": "%",
                    },
                    "Insert Size": {
                        "title": "Insert Size",
                        "description": "Average insert size from Samtools stats",
                    },
                    "Insert Size s.d.": {
                        "title": "Insert Size s.d.",
                        "description": "Insert size standard deviation from Samtools stats",
                    },
            }
    }

id_json = {
            "id": "qc_table",
            "section_name": "QC Stats",
            "description": "QC-values from Picard, Samtools and Mosdepth",
            "plot_type": "table",
            "pconfig": {
                "namespace": "qc-table"
            },

}
json_output = {**id_json, **header_json, **data_json_ordered}
with open(snakemake.output.json, 'w+') as outfile:
    json.dump(json_output, outfile, ensure_ascii=False, indent=4)

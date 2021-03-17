#!/bin/python3.6

import sys
import csv

batchFile = sys.argv[1]
SampleSheet = sys.argv[2]
outFile = sys.argv[3]  # should end with _mqc.json for MultiQC

# Get all lines from get_stats.py output
with open(batchFile, 'r') as bFile:
    linesUnordered = [line.strip().split(',') for line in bFile]
unOrdSamples = [x[0] for x in linesUnordered]

# Get the order from SampleSheetUsed
samples = []
startReading = 0
with open(SampleSheet, 'r') as file:
    lines = [line.strip() for line in file]
    for line in lines:
        if startReading == 1:  # Once reached [Data]
            samples.append(line.split(',')[1])
        if line.startswith("[Data]"):
            startReading = 1
# samples.pop() #Remove any empty are there empty line at end?!
samples = samples[1:]  # Remove header from SampleSheet
sampleSheetSamples = [string for string in samples if string != ""]  # Remove empty fields
# Remove any HD829 because other pipeline
# HDindices = [i for i, x in enumerate(sampleSheetSamples) if x.startswith("HD829")]
# if len(HDindices) != 0 :
#     for index in HDindices:
#         sampleSheetSamples.pop(index)

header = ['Sample', 'Tot seq', 'Reads mapped', 'Avg Coverage', 'Breadth 500x', 'Reads paired [%]', 'Insert size',
          'Insert size s.d.', 'Average quality', 'Duplicates [%]', 'Breadth 50x', 'Breadth 100x', 'Bases on target']

with open(outFile, 'w') as file:
    # Write all config for custom table
    file.write("{\n")
    file.write("  \"id\": \"qc_table\",\n")
    file.write("  \"section_name\": \"QC stats\",\n")
    file.write("  \"description\": \"QC-values from Picard, Samtools and CARTool\",\n")
    file.write("  \"plot_type\": \"table\",\n")
    file.write("  \"pconfig\": {\n")
    file.write("    \"namespace\": \"qc-table\"\n")
    file.write("  },\n")

    file.write("  \"headers\": {\n")  # All header configs
    file.write('    \"Avg Coverage\": {\n')
    file.write("      \"title\": \"Average Coverage\",\n")
    file.write("      \"description\": \"Avg cov of bedfile from CARTool\"\n")
    file.write("    },\n")
    file.write("    \"Average quality\": {\n")
    file.write("      \"title\": \"Average quality\",\n")
    file.write("      \"description\": \"Average mapping quality from Samtools\",\n")
    file.write("      \"min\": 0,\n")
    file.write("      \"max\": 60,\n")
    file.write("      \"scale\": \"RdYlGn\"\n")
    file.write("    },\n")
    file.write("    \"Bases on target\": {\n")
    file.write("      \"title\": \"Bases on target\",\n")
    file.write("      \"description\": \"Ratio of bases on target from Picard HsMetrics\",\n")
    file.write("      \"format\": \"{:,.3f}\"\n")
    file.write("    },\n")
    file.write("    \"Duplicates [%]\": {\n")
    file.write("      \"title\": \"Duplicates\",\n")
    file.write("      \"description\": \"Percent duplicates marked by MarkDuplicates in Picard\",\n")
    file.write("      \"min\": 0,\n")
    file.write("      \"max\": 50,\n")
    file.write("      \"scale\": \"RdYlGn-rev\",\n")
    file.write("      \"suffix\": \"%\"\n")
    file.write("    },\n")
    file.write("    \"Breadth 500x\": {\n")
    file.write("      \"title\": \"Coverage breadth 500x\",\n")
    file.write("      \"description\": \"Design covered to 500x from CARTool\",\n")
    file.write("      \"min\": 0,\n")
    file.write("      \"max\": 1,\n")
    file.write("      \"scale\": \"RdYlGn\",\n")
    file.write("      \"format\": \"{:,.2f}\"\n")
    file.write("    },\n")
    file.write("    \"Breadth 100x\": {\n")
    file.write("      \"title\": \"Coverage breadth 100x\",\n")
    file.write("      \"description\": \"Design covered to 100x or over from Picard HsMetrics\",\n")
    file.write("      \"format\": \"{:,.2f}\"\n")
    file.write("    },\n")
    file.write("    \"Breadth 50x\": {\n")
    file.write("      \"title\": \"Coverage breadth 50x\",\n")
    file.write("      \"description\": \"Design covered to 50x or over from Picard HsMetrics\",\n")
    file.write("      \"format\": \"{:,.2f}\"\n")
    file.write("    },\n")
    file.write('    \"Tot seq\": {\n')
    file.write("      \"title\": \"Total sequences\",\n")
    file.write("      \"description\": \"Number of reads in fastq from Samtools\",\n")
    file.write("      \"format\": \"{:,.0f}\"\n")
    file.write('    },\n')
    file.write('    \"Reads mapped\": {\n')
    file.write("      \"title\": \"Reads mapped\",\n")
    file.write("      \"description\": \"Number of reads mapped from Samtools\",\n")
    file.write("      \"format\": \"{:,.0f}\"\n")
    file.write("    },\n")
    file.write("    \"Reads paired [%]\": {\n")
    file.write("      \"title\": \"Reads properly paired\",\n")
    file.write("      \"description\": \"Percent properly paired reads from Samtools\",\n")
    file.write("      \"min\": 0,\n")
    file.write("      \"max\": 100,\n")
    file.write("      \"scale\": \"RdYlGn\",\n")
    file.write("      \"suffix\": \"%\"\n")
    file.write("    },\n")
    file.write("    \"Insert size\": {\n")
    file.write("      \"title\": \"Insert size\",\n")
    file.write("      \"description\": \"Averge insert size from Samtools\"\n")
    file.write("    },\n")
    file.write("    \"Insert size s.d.\": {\n")
    file.write("      \"title\": \"Insert size s.d.\",\n")
    file.write("      \"description\": \"Insert size standard deviation from Samtools\"\n")
    file.write("    }\n")
    file.write("  },\n")  # Close headers

    # The table in samplesheet order
    file.write("  \"data\": {\n")
    for sample in sampleSheetSamples:
        line = linesUnordered[unOrdSamples.index(sample)]
        # 'sample1' :{ 4spaces
        file.write("    \""+sample+"\": {\n")
        # 'col1': value1, 6spaces
        for i in range(1, len(header)):
            if i == len(header)-1:
                file.write("      \""+header[i]+"\": " + line[i] + "\n")
            else:
                file.write("      \""+header[i]+"\": " + line[i] + ",\n")
        if sample == sampleSheetSamples[-1]:
            file.write('    }\n')
        else:
            file.write('    },\n')
    file.write('  }\n')  # close data

    file.write('}')

#!/bin/python3.6
import sys
import subprocess
import csv

picardDup = sys.argv[1]
picardMet = sys.argv[2]
samtools = sys.argv[3]
# multiQCheader = sys.argv[4]
cartoolLog = sys.argv[4]
sample = sys.argv[5]
# outFile = sys.argv[7]
batchFile = sys.argv[6]

# sample = samtools.split('/')[1]

# picardDup
# percent duplicateLevel
duplCmd = 'grep -A1 PERCENT '+picardDup+' | tail -1 | cut -f9'
duplicateLevel = subprocess.run(duplCmd, stdout=subprocess.PIPE, shell='TRUE').stdout.decode(
    'utf-8').strip()  # need *100 to be percent
dupliPerc = str(round(float(duplicateLevel)*100, 2))
# picardMet
# insert?, bases on target, coverage 50,100,500x,
metCmd = 'grep -A1 BAIT_SET '+picardMet
met = subprocess.run(metCmd, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8')
metrics = met.split('\n')
zipObject = zip(metrics[0].split('\t'), metrics[1].split('\t'))
metricsDict = dict(zipObject)
# metricsDict['PCT_SELECTED_BASES'] ##bases on target
# metricsDict['PCT_TARGET_BASES_50X'] #
# metricsDict['PCT_TARGET_BASES_100X']
# samtools stats
# readsmapped,raw tot seq, reads paired %, insert +sd, avg qual
samCmd = 'grep SN '+samtools
sam = subprocess.run(samCmd, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
sam = sam.split('\nSN\t')
listOfList = [i.split('\t') for i in sam[1:]]
samDict = {item[0].strip(':'): item[1] for item in listOfList}
# samDict['raw total sequences']
# samDict['reads mapped']
# samDict['percentage of properly paired reads (%)']
# samDict['insert size average']
# samDict['insert size standard deviation']
# samDict['average quality']

avgCovCmd = 'grep "Mean Coverage Depth:" '+cartoolLog + ' |cut -f2 -d"," | cut -f1 -d" " '
avgCov = subprocess.run(avgCovCmd, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8')

breadth500Cmd = 'grep "Mean Coverage Breadth:" '+cartoolLog + ' | cut -f3 -d"," '
breadth500 = subprocess.run(breadth500Cmd, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8')

header = ['Sample', 'Tot seq', 'Reads mapped', 'Avg Coverage', 'Breadth 500x', 'Reads paired [%]', 'Insert size',
          'Insert size s.d.', 'Average Quality', 'Duplicates [%]', 'Target bases 50x', 'Target bases 100x', 'Bases on target']
line = [sample, samDict['raw total sequences'], samDict['reads mapped'], avgCov.strip(), breadth500.strip(),
        samDict['percentage of properly paired reads (%)'], samDict['insert size average'],
        samDict['insert size standard deviation'], samDict['average quality'], dupliPerc, metricsDict['PCT_TARGET_BASES_50X'],
        metricsDict['PCT_TARGET_BASES_100X'], metricsDict['PCT_SELECTED_BASES']]

# append to Batch file and write sampleFile
with open(batchFile, 'a') as file:
    writer = csv.writer(file, delimiter=',', lineterminator='\n')
    writer.writerow(line)

# Print multiQCheader
# with open(multiQCheader, 'r') as f:
#     with open(outFile, "w") as file:
#         for mqcline in f:
#             file.write(mqcline)
#         writer = csv.writer(file, delimiter=',', lineterminator='\n')
#         writer.writerow(header)
#         writer.writerow(line)

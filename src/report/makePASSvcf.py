#!/bin/python3.6
import sys
import subprocess
from pysam import VariantFile

# Define sys.argvs
vcf_in = VariantFile(sys.argv[1])  # variatnts
pindel_in = VariantFile(sys.argv[2])
artefactFile = sys.argv[3]
germlineFile = sys.argv[4]
hematoCountFile = sys.argv[5]

new_header = vcf_in.header

vcf_out = VariantFile(sys.argv[6], 'w', header=new_header)

for x in vcf_in.header.records:
    if 'CSQ' in str(x):
        csqIndex = str(x).split('Format: ')[1].strip().strip('">').split('|')

# SNVs
for record in vcf_in.fetch():
    synoCosmicN = 0
    spliceVariant = False
    csq = record.info["CSQ"][0].split("|")
    consequence = csq[csqIndex.index('Consequence')]
    if record.filter.keys() == ["Syno"]:  # If only syno! No popAF.    any(x in "Syno" for x in record.filter.keys()):
        csq = record.info["CSQ"][0]
        synoCosmicVepList = [cosmic for cosmic in csq.split("|")[17].split(
            "&") if cosmic.startswith('CO')]  # Get all cosmicID in list
        if len(synoCosmicVepList) != 0:
            for synoCosmicId in synoCosmicVepList:
                cmdCosmic = 'grep -w '+synoCosmicId+' '+hematoCountFile+' | cut -f 16 '
                synoCosmicNew = subprocess.run(cmdCosmic, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
                if len(synoCosmicNew) == 0:
                    synoCosmicNew = 0
                synoCosmicN += int(synoCosmicNew)
        if 'splice' in consequence:
            spliceVariant = True

    if record.filter.keys() == ["PASS"] or synoCosmicN != 0 or spliceVariant:
        # if len(record.ref) > len(record.alts[0]): Deletion but need to remove first base as well
        #     record.pos = record.pos-1
        if record.info["AF"][0] >= 0.03:  # Change for pindel, get all pindels
            cmdArt = 'grep -w '+str(record.pos)+' '+artefactFile
            artLines = subprocess.run(cmdArt, stdout=subprocess.PIPE, shell='TRUE').stdout.decode(
                'utf-8').strip()  # What happens if two hits?
            artefact_variant = 0

            for artLine in artLines.split("\n"):
                # if pos exists and match in artefact file.
                if artLine and record.ref == artLine.split()[2] and record.alts[0] == artLine.split()[3]:
                    # Artefact do not print
                    artefact_variant = 1
                    continue

            if artefact_variant == 0:
                cmdGerm = 'grep -w '+str(record.pos)+' '+germlineFile
                germLines = subprocess.run(cmdGerm, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
                germline_variant = 0
                for germLine in germLines.split("\n"):
                    # if exists in germline file
                    if germLine and record.ref == germLine.split("\t")[2] and record.alts[0] == germLine.split("\t")[3]:
                        # Germline match, do nothing
                        germline_variant = 1
                        continue
                if germline_variant == 0:
                    vcf_out.write(record)

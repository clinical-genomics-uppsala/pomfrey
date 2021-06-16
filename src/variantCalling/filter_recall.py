#!/bin/python3.6
import sys
from pysam import VariantFile
import subprocess

vcf_in = VariantFile(sys.argv[1])
new_header = new_header = vcf_in.header
vcf_out = VariantFile(sys.argv[2], 'w', header=new_header)
sv_out = sys.argv[2]+'.svtypeDEL.txt'
indelArteFile = sys.argv[3]

for record in vcf_in.fetch():
    # Remove strange vardict DEL variants
    try:
        if record.info["SVTYPE"] == 'DEL':
            with open(sv_out, 'a+') as svtype_out:
                svtype_out.write(str(record))
    except KeyError:
        if len(record.ref) != len(record.alts[0]):  # if InDel
            # Support by either Vardict or Manta, ok.
            if ("Mutect2".lower() in [i.lower() for i in record.info["CALLERS"]] or
               "Vardict".lower() in [i.lower() for i in record.info["CALLERS"]]):
                # Check if indel artefact
                write = 1
                cmdIndelArte = 'grep -w '+str(record.pos)+' '+indelArteFile
                artefactLines = subprocess.run(cmdIndelArte, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
                for artefactLine in artefactLines.split("\n"):
                    # if ref and alt is same as artefacts
                    if artefactLine and record.ref == artefactLine.split()[2] and record.alts[0] == artefactLine.split()[3]:
                        write = 0  # Could exit for loop?
                if write == 1:
                    vcf_out.write(record)
        elif len(record.info["CALLERS"]) > 2:  # if substitution demand 3/4 support
            vcf_out.write(record)

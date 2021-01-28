#!/bin/python3.6
import sys
import csv
from pysam import VariantFile
import subprocess

vcf_in = VariantFile(sys.argv[1])
multiVcf = VariantFile(sys.argv[2])
new_header = vcf_in.header
# new_header.generic.add("Multi allelic variants added from Pisces.")

vcf_out = VariantFile(sys.argv[3], 'w', header=new_header)

for record in vcf_in.fetch():
    vcf_out.write(record)
    for mRecord in multiVcf.fetch():
        if record.contig == mRecord.contig and record.pos == mRecord.pos:
            # import pdb; pdb.set_trace()
            if record.alts[0] != mRecord.alts[0]:
                vcf_out.write(mRecord)

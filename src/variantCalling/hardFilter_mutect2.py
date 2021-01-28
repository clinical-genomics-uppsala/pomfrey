#!/bin/python3.6
import sys
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1])
vcf_out = VariantFile(sys.argv[2], 'w', header=vcf_in.header)

for record in vcf_in.fetch():
    if record.filter.keys() == ["PASS"]:
        vcf_out.write(record)

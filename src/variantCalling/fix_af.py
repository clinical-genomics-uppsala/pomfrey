#!/bin/python3.6
import sys
import re
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1])  # dosen't matter if bgziped or not. Automatically recognizes
method = re.search('callers/(.+?)/', sys.argv[1]).group(1)  # The folder after callers/
# Add new filter descriptions to new header
new_header = vcf_in.header
# import pdb; pdb.set_trace()

if method == "freebayes" or method == "pisces":  # Byta description on AF for freebayes?
    sample = vcf_in.header.samples[0]

if method == "pisces" or method == "mutect2":
    new_header.info.add("AF", "A", "Float", "DescriptionDescription")

if method == "snver":
    new_header.info.add("AF", "A", "Float", "Allel count divided on depth, crude")

# start new vcf with the new_header
vcf_out = VariantFile(sys.argv[2], 'w', header=new_header)

for record in vcf_in.fetch():
    if method == "freebayes":
        ad = record.samples[sample].get("AD")
        af = []
        for item in ad:
            af.append(item/sum(ad))
    if method == "pisces":
        af = record.samples[sample].get("VF")
    if method == "mutect2" or method == "vardict":
        af = record.samples[0].get("AF")
    if method == "snver":
        dp = record.info["DP"]
        ac = record.info["AC"]
        af = ac/dp

    record.info["AF"] = tuple(af)

    vcf_out.write(record)

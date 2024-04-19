#!/bin/python3
import re
from pysam import VariantFile
import gzip

vcf_in = VariantFile(snakemake.input.vcf)
caller = re.search('callers/(.+?)/', snakemake.input.vcf).group(1)  # The folder after callers/
# Add new filter descriptions to new header
new_header = vcf_in.header
if caller == "pisces" or caller == "mutect2":
    new_header.info.add("AF", "A", "Float", "DescriptionDescription")


# start new vcf with the new_header
vcf_out = VariantFile(snakemake.output.vcf, 'w', header=new_header)

for record in vcf_in.fetch():
    if caller == "freebayes":
        ads = record.samples[0].get("AD")
        ad_afs = []
        for ad in ads:
            ad_afs.append(ad/sum(ads))
        af = tuple(ad_afs[1:])
    elif caller == "pisces":
        af = record.samples[0].get("VF")
    elif caller == "mutect2" or caller == "vardict":
        af = record.samples[0].get("AF")
    else:
        raise ValueError(
            "{} is not a valid caller for this script. Choose between:"
            "freebayes, mutect2, pisces, vardict.".format(caller)
        )
    record.info["AF"] = af

    vcf_out.write(record)

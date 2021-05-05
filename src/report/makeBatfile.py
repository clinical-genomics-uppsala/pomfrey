#!/bin/python3.6
import sys
from pysam import VariantFile

batfile = open(sys.argv[1], 'w')
vcf_in = VariantFile(sys.argv[2], 'r')
bamfile = sys.argv[3]
reffile = sys.argv[4]
bedfile = sys.argv[5]
outfolder = sys.argv[6]
padding = int(sys.argv[7])
sort = sys.argv[8]
view = sys.argv[9]
format = sys.argv[10]

batfile.write("new")
batfile.write("\ngenome "+reffile)
batfile.write("\nload "+bedfile+"\nload "+bamfile)
batfile.write("\nsnapshotDirectory "+outfolder)
batfile.write("\npreference SAM.SHOW_SOFT_CLIPPED false")
batfile.write("\npreference SAM.DOWNSAMPLE false")  # should we keep? goes outside picture usually


for record in vcf_in.fetch():
    # vcfRow = line.split()
    chr = record.contig
    pos = record.pos-1
    low = pos-padding
    high = pos+len(record.alts[0])+padding
    csq = record.info["CSQ"][0]
    gene = csq.split("|")[3]

    batfile.write("\ngoto "+chr+":"+str(low)+"-"+str(high))
    if len(record.ref) > len(record.alts[0]):  # Deletions
        sortPos = record.pos+1
    else:
        sortPos = record.pos
    batfile.write("\nsort "+sort+" "+str(sortPos)+"\n"+view+" "+bamfile.split("/")[-1])
    batfile.write("\nsnapshot "+gene+"-"+chr+"_"+str(pos)+"_"+str(pos+len(record.alts[0]))+"."+format)

batfile.write("\nexit")
batfile.close()

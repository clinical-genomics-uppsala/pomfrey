import subprocess
import gzip

threshold_table = []
total_min_breadth = 0
total_med_breadth = 0
total_max_breadth = 0
total_length = 0
with open(snakemake.output.fixed, 'w+') as outfile:
    with gzip.open(snakemake.input.mosdepth_thresholds, 'rt') as thresfile:
        next(thresfile)
        for lline in thresfile:
            line = lline.strip().split("\t")
            length = int(line[2])-int(line[1])
            total_length += length
            total_min_breadth += int(line[4])
            total_med_breadth += int(line[5])
            total_max_breadth += int(line[6])

            min_bredth = round(int(line[4])/length, 4)
            med_breadth = round(int(line[5])/length, 4)
            max_breadth = round(int(line[6])/length, 4)

            outline = line[0:4]+[str(min_bredth), str(med_breadth), str(max_breadth)]
            outfile.write("\t".join(outline)+"\n")

with open(snakemake.output.summary, "w+") as summaryfile:
    outline = [str(round(total_min_breadth/total_length, 4)),
               str(round(total_med_breadth/total_length, 4)), str(round(total_max_breadth/total_length, 4))]
    summaryfile.write("\t".join(outline)+"\n")

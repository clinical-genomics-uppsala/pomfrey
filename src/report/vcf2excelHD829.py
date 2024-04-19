#!/bin/python3
import csv
from pysam import VariantFile
import xlsxwriter
from datetime import date
import subprocess
from operator import itemgetter
import yaml
import gzip

# Known variation in HD829: gene,variat, pos, type, cosmic, af
known = [
    ["ABL1", "T315I", "133748283", "SNP", "COSM12560", "0.050"],
    ["ASXL1", "G646fs*12", "31022441", "INS", "COSM1411076", "0.400"],
    ["ASXL1", "W796C", "31022903", "SNP", "COSM1681610", "0.050"],
    ["BCOR", "Q1174fs*8", "39923086", "INS", "COSM3732385", "0.700"],
    ["CBL", "S403F", "119148988", "SNP", "COSM1676499", "0.050"],
    ["DNMT3A", "R882C", "25457243", "SNP", "COSM53042", "0.050"],
    ["EZH2", "R418Q", "148514471", "SNP", "COSM3259655", "0.050"],
    ["FLT3", "D835Y", "28592642", "SNP", "COSM783", "0.050"],
    ["FLT3", "ITD300", "28608047", "300bp INS", "N/A", "0.050"],
    ["GATA1", "Q119*", "48650385", "SNP", "N/A", "0.100"],
    ["GATA2", "G200fs*18", "128204841", "DEL", "COSM1418772", "0.350"],
    ["IDH1", "R132C", "209113113", "SNP", "COSM28747", "0.050"],
    ["IDH2", "R172K", "90631838", "SNP", "COSM33733", "0.050"],
    ["JAK2", "F537-K539>L", "5070021", "DEL", "COSM24437", "0.050"],
    ["JAK2", "V617F", "5073770", "SNP", "COSM12600", "0.050"],
    ["KRAS", "G13D", "25398281", "SNP", "COSM532", "0.400"],
    ["NPM1", "W288fs*12", "170837543", "INS", "COSM158604", "0.050"],
    ["NRAS", "Q61L", "115256529", "SNP", "COSM583", "0.100"],
    ["RUNX1", "M267I", "36206711", "SNP", "COSM1681955", "0.350"],
    ["SF3B1", "G740E", "198266713", "SNP", "COSM133120", "0.050"],
    ["TET2", "R1261H", "106164914", "SNP", "COSM211643", "0.050"],
    ["TP53", "S241F", "7577559", "SNP", "COSM10812", "0.050"],
]

knownPos = [x[2] for x in known]
known_found_temp = []
known_found = []

vcf_snv = VariantFile(snakemake.input.vcf_snv)
vcf_indel = VariantFile(snakemake.input.vcf_pindel)

runid = snakemake.params.seqid
sample = list(vcf_snv.header.samples)[0]
today = date.today()
emptyList = ["", "", "", "", "", ""]


# VEP fileds in list to get index
def index_vep(variantfile):
    csqIndex = []
    for x in variantfile.header.records:
        if "CSQ" in str(x):
            csqIndex = str(x).split("Format: ")[1].strip().strip('">').split("|")
    return csqIndex


# Return matching lines in file
def extractMatchingLines(expressionMatch, artefactFile, grepvarible):
    if grepvarible == "":
        grepvarible = "-wE "
    cmdArt = "grep " + grepvarible + " " + str(expressionMatch) + " " + artefactFile
    matchLines = subprocess.run(cmdArt, stdout=subprocess.PIPE, shell="TRUE").stdout.decode("utf-8").strip()
    return matchLines


def file_length(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


""" Loop through vcf_snv to find known (and unknown) variants """
for x in vcf_snv.header.records:
    if x.key == "reference":
        refV = x.value
    if x.key == "VEP":
        vepline = x.value

white = []
green = []
orange = []
underFive = []  # put after green and orange but still white
csqIndex = index_vep(vcf_snv)


for record in vcf_snv.fetch():
    # Get annotation from vep field
    csq = record.info["CSQ"][0].split("|")
    gene = csq[csqIndex.index("SYMBOL")]
    consequence = csq[csqIndex.index("Consequence")]
    clinical = csq[csqIndex.index("CLIN_SIG")]
    existing = csq[csqIndex.index("Existing_variation")].split("&")

    if len([rs for rs in existing if rs.startswith("rs")]) == 0:
        rs = ""
    else:
        rs = ", ".join([rs for rs in existing if rs.startswith("rs")])

    try:
        if record.info["CALLERS"]:
            callers = " & ".join(record.info["CALLERS"])
    except KeyError:
        callers = "Pisces-multi"

    transcript = csq[csqIndex.index("HGVSc")].split(":")[0]
    if len(csq[csqIndex.index("HGVSc")].split(":")) > 1:
        codingName = csq[csqIndex.index("HGVSc")].split(":")[1]
    else:
        codingName = ""
    ensp = csq[csqIndex.index("HGVSp")]

    popFreqsPop = csqIndex[csqIndex.index("AF") : csqIndex.index("gnomAD_SAS_AF") + 1]
    popFreqAllRaw = csq[csqIndex.index("AF") : csqIndex.index("gnomAD_SAS_AF") + 1]
    if any(popFreqAllRaw) and max([float(x) if x else 0 for x in popFreqAllRaw]) != 0:  # if all not empty
        popFreqAll = [float(x) if x else 0 for x in popFreqAllRaw]
        maxPopAf = max(popFreqAll)
        maxIndex = [i for i, j in enumerate(popFreqAll) if j == maxPopAf]
        if len(maxIndex) == 1:
            maxPop = popFreqsPop[maxIndex[0]]
        else:
            popFreqPops = [popFreqsPop[x] for x in maxIndex]
            maxPop = "&".join(popFreqPops)
    else:
        maxPopAf = ""
        maxPop = ""

    if len(record.info["AF"]) == 1:
        af = record.info["AF"][0]
    else:
        print(record.info["AF"])
        sys.exit()  # Fails if vt decompose didn't work

    if len(record.alts) == 1:
        alt = record.alts[0]
    else:
        print(record.alts)
        sys.exit()  # Fails if vt decompose didn't work

    # Add variants if syno and in cosmic hemto list
    synoCosmicN = 0
    spliceVariant = False
    if record.filter.keys() == ["Syno"]:  # Only if Syno not and popAF.
        synoCosmicVepList = [
            cosmic for cosmic in csq[csqIndex.index("Existing_variation")].split("&") if cosmic.startswith("CO")
        ]  # Get all cosmicID in list
        # COSMIC Hemato hit
        if len(synoCosmicVepList) != 0:
            for synoCosmicId in synoCosmicVepList:
                try:
                    synoCosmicNew = extractMatchingLines(synoCosmicId, snakemake.input.hemato_count, "-wE").split("\t")[15]
                except IndexError:
                    synoCosmicNew = 0
                synoCosmicN += int(synoCosmicNew)
        if "splice" in consequence:
            spliceVariant = True

    if record.filter.keys() == ["PASS"] or synoCosmicN != 0 or spliceVariant:
        # Test if variant in known list
        for known_line in known:
            if gene == known_line[0] and str(record.pos) == known_line[2]:
                known_found_temp.append(known_line + [af, record.info["DP"], record.ref, alt])

        # Total number of cosmic hemato hits on the position. Vep reports all cosmicId for that position.
        cosmicVepList = [cosmic for cosmic in existing if cosmic.startswith("CO")]
        if len(cosmicVepList) == 0:
            cosmicVep = ""
        else:
            cosmicVep = ", ".join(cosmicVepList)

        if len(cosmicVepList) == 0:
            cosmicN = ""
        else:
            cosmicN = 0
            for cosmicId in cosmicVepList:
                try:
                    cosmicNew = extractMatchingLines(cosmicId, snakemake.input.hemato_count, "-wE").split("\t")[15]
                except IndexError:
                    cosmicNew = 0
                cosmicN += int(cosmicNew)

        if float(af) >= 0.01:
            snv = [
                runid,
                sample,
                gene,
                record.contig,
                record.pos,
                record.ref,
                alt,
                af,
                record.info["DP"],
                transcript,
                codingName,
                ensp,
                consequence,
                cosmicVep,
                cosmicN,
                clinical,
                rs,
                maxPopAf,
                maxPop,
                callers,
            ]

            # Append line with sample and rundate to rolling list of artefacts..
            with open(snakemake.input.variantslog, "a") as appendfile:
                variants = snv + ["\n"]
                appendfile.write("\t".join(str(e) for e in variants))

            # Check if variant in artefact or germline file
            artLines = extractMatchingLines(str(record.pos), snakemake.input.artefact_snv, "-wE")
            artefact_variant = False

            for artLine in artLines.split("\n"):
                # if pos exists and match in artefact file.
                if artLine and record.ref == artLine.split()[2] and alt == artLine.split()[3]:
                    orange.append(snv)
                    artefact_variant = True
                    break
            if not artefact_variant:
                germLines = extractMatchingLines(str(record.pos), snakemake.input.germline, "-wE")
                germline_variant = False
                for germLine in germLines.split("\n"):
                    # if exists in germline file
                    if germLine and record.ref == germLine.split()[2] and alt == germLine.split()[3]:
                        green.append(snv)
                        germline_variant = True
                        break
                if not germline_variant:
                    if float(af) < 0.05:
                        underFive.append(snv)
                    else:
                        white.append(snv)


""" Loop through vcf_indel for unknown variants and FLT3 ITD """
itd_lines = []
orangeIndel = []
whiteIndel = []
underFiveIndel = []
csqIndex = index_vep(vcf_indel)

for indel in vcf_indel.fetch():
    if indel.filter.keys() == ["PASS"]:
        svlen = indel.info["SVLEN"]
        af = indel.info["AF"]

        if len(indel.alts) == 1:
            alt = indel.alts[0]
        else:
            print(indel.alts)
            sys.exit()

        csqIndel = indel.info["CSQ"][0].split("|")  # VEP annotation
        indelGene = csqIndel[csqIndex.index("SYMBOL")]

        # Not using ExAC pop
        popFreqsPop = csqIndex[csqIndex.index("AF") : csqIndex.index("gnomAD_SAS_AF") + 1]
        popFreqAllRawIndel = csqIndel[csqIndex.index("AF") : csqIndex.index("gnomAD_SAS_AF") + 1]
        if any(popFreqAllRawIndel) and max([float(x) if x else 0 for x in popFreqAllRawIndel]) != 0:  # if all not empty
            popFreqAllIndel = [float(x) if x else 0 for x in popFreqAllRawIndel]
            maxPopAfIndel = max(popFreqAllIndel)
            maxIndexIndel = [i for i, j in enumerate(popFreqAllIndel) if j == maxPopAfIndel]
            if len(maxIndexIndel) == 1:
                maxPopIndel = popFreqsPop[maxIndexIndel[0]]
            else:
                popFreqPopsIndel = [popFreqsPop[x] for x in maxIndexIndel]
                maxPopIndel = "&".join(popFreqPopsIndel)
        else:
            maxPopAfIndel = ""
            maxPopIndel = ""

        indelTranscript = csqIndel[csqIndex.index("HGVSc")].split(":")[0]
        if len(csqIndel[csqIndex.index("HGVSc")].split(":")) > 1:
            indelCodingName = csqIndel[csqIndex.index("HGVSc")].split(":")[1]
        else:
            indelCodingName = ""
        indelEnsp = csqIndel[csqIndex.index("HGVSp")]

        indelRow = [
            runid,
            sample,
            indelGene,
            indel.contig,
            indel.pos,
            indel.stop,
            svlen,
            af,
            indel.ref,
            alt,
            indel.info["DP"],
            indelTranscript,
            indelCodingName,
            indelEnsp,
            maxPopAfIndel,
            maxPopIndel,
        ]

        # If FLT3 ITD from known list
        if indel.pos > 28608040 and indel.pos < 28608054:
            itd_lines.append([af, indel.info["DP"], indel.ref, alt, svlen])

        # Mark artefact based on artefactfile
        artLines = extractMatchingLines(str(indel.contig) + ".*" + str(indel.pos), snakemake.input.artefact_pindel, "-wE")

        if len(artLines) > 0:  # if pos exists and match in artefact file.
            orangeIndel.append(indelRow)
        else:
            if float(af) < 0.05:
                underFiveIndel.append(indelRow)
            else:
                whiteIndel.append(indelRow)


""" Organize known dictionary to match known list for known sheet """
# Order known found if several variants in same pos and order based on odre of known.
known_found_temp = sorted(known_found_temp, key=itemgetter(0))
if len(known_found_temp) < len(known):
    i = 0
    for known_variant in known:
        if known_found_temp[i] and known_found_temp[i][2] != known_variant[2]:
            known_found.append(known_variant)
        else:
            known_found.append(known_found_temp[i])
            i += 1
else:
    known_found = known_found_temp

# Add FLT3 ITD from pindel
num_known = 0
last_line = ["", "", "", "", "", ""]
known_lines_final = []
first_found_itd = False
for known_line in known_found:
    line = [runid, sample] + known_line
    if known_line[1] == "ITD300":
        numITD = 1
        if len(itd_lines) < 1:
            known_lines_final.append(line)
        else:
            for itd in itd_lines:
                if not first_found_itd:
                    known_lines_final.append(line + itd)
                    first_found_itd = True
                    num_known += 1
                else:
                    known_lines_final.append(["", "", "", "", "", "", "", ""] + itd)
    elif len(known_line) < 7:  # Saknas
        known_lines_final.append(line)
    elif last_line[4] == line[4]:  # Forra raden ifall fler varianter pa samma pos
        known_lines_final.append(["", "", "", "", "", "", "", ""] + line[8:])
    else:
        known_lines_final.append(line)
        num_known += 1
    last_line = line


""" Coverage of all regions in bedfile """
min_cov = int(snakemake.params.thresholds.split(",")[0])
med_cov = int(snakemake.params.thresholds.split(",")[1])
max_cov = int(snakemake.params.thresholds.split(",")[2])

bed_table = []
with open(snakemake.input.bedfile, "r") as bedfile:
    next(bedfile)
    # Skip header
    for line in bedfile:
        bed_table.append(line.strip().split("\t"))
# Coverage per region for Coverage sheet
cov_table_lines = []
with gzip.open(snakemake.input.mosdepth_regions, "rt") as regionsfile:
    for lline in regionsfile:
        line = lline.strip().split("\t")
        cov_row = [line[3], line[0], line[1], line[2], line[4]]
        cov_table_lines.append(cov_row)


""" Low cov file """
low_cov_lines = []
condensed_line = ["", "", "", "", ""]
with open(snakemake.input.mosdepth_lowcov, "r") as lowfile:
    for lline in lowfile:
        line = lline.strip().split("\t")
        if int(line[3]) <= int(min_cov):
            if condensed_line[0] == line[0] and condensed_line[2] == line[1]:
                condensed_length = int(condensed_line[2]) - int(condensed_line[1])
                line_length = int(line[2]) - int(line[1])
                new_avgcov = (float(condensed_line[3]) * condensed_length + float(line[3]) * line_length) / (
                    line_length + condensed_length
                )
                condensed_line[3] = new_avgcov
                condensed_line[2] = line[2]
            else:
                low_cov_lines.append(condensed_line)
                condensed_line = line

low_cov_lines.pop(0)
for line in low_cov_lines:
    for bedline in bed_table:  # blir det en rad eller element?
        if line[0] == bedline[0] and int(line[1]) >= int(bedline[1]) and int(line[2]) <= int(bedline[2]):
            line[:] = [bedline[3]] + line[:] + [int(line[2]) - int(line[1])]
            break

# Sort based on coverage
low_cov_lines.sort(key=lambda x: float(x[4]))

# Number of low cov regions for the Overview sheet.
low_regions = len(low_cov_lines)


""" Hotspot file """
low_pos = 0
hotspotTable = []
with open(snakemake.input.mosdepth_hotspot, "r") as hotFile:
    for dpLine in hotFile:
        line = dpLine.strip().split("\t")
        chr = line[0]
        pos = line[2]
        dp = line[3]
        for bedline in bed_table:  # blir det en rad eller element?
            if chr == bedline[0] and int(pos) >= int(bedline[1]) and int(pos) <= int(bedline[2]):
                gene = bedline[3]
                break

        hotspotTable.append([str(gene), chr, str(pos), str(dp)])
        if int(dp) <= med_cov:
            low_pos += 1
hotspotTable.sort(key=lambda x: float(x[2]))


""" Create xlsx file and sheets. """
workbook = xlsxwriter.Workbook(snakemake.output[0])
worksheetOver = workbook.add_worksheet("Overview")
worksheetKnown = workbook.add_worksheet("Known")
worksheetSNV = workbook.add_worksheet("SNVs")
worksheetIndel = workbook.add_worksheet("InDel")
worksheetLowCov = workbook.add_worksheet("Low Coverage")
worksheetHotspot = workbook.add_worksheet("Hotspot")
worksheetCov = workbook.add_worksheet("Coverage")
worksheetQCI = workbook.add_worksheet("QCI")
worksheetVersions = workbook.add_worksheet("Version")
# Define formats to be used.
headingFormat = workbook.add_format({"bold": True, "font_size": 18})
lineFormat = workbook.add_format({"top": 1})
tableHeadFormat = workbook.add_format({"bold": True, "text_wrap": True})
textwrapFormat = workbook.add_format({"text_wrap": True})
italicFormat = workbook.add_format({"italic": True})
redFormat = workbook.add_format({"font_color": "red"})

greenFormat = workbook.add_format({"bg_color": "#85e085"})
orangeFormat = workbook.add_format({"bg_color": "#ffd280"})
green_italicFormat = workbook.add_format({"bg_color": "#85e085", "italic": "True"})
orange_italicFormat = workbook.add_format({"bg_color": "#ffd280", "italic": "True"})

# Overview
worksheetOver.write(0, 0, sample, headingFormat)
worksheetOver.write(1, 0, "RunID: " + runid)
worksheetOver.write(2, 0, "Processing date: " + today.strftime("%B %d, %Y"))
worksheetOver.write_row(3, 0, emptyList, lineFormat)

worksheetOver.write(4, 0, "Created by: ")
worksheetOver.write(4, 4, "Valid from: ")
worksheetOver.write(5, 0, "Signed by: ")
worksheetOver.write(5, 4, "Document nr: ")
worksheetOver.write_row(6, 0, emptyList, lineFormat)

worksheetOver.write(7, 0, "Sheets:", tableHeadFormat)
worksheetOver.write_url(8, 0, "internal:'Known'!A1", string="Known Variants")
worksheetOver.write_url(9, 0, "internal:'SNVs'!A1", string="Variants analysis")
worksheetOver.write_url(10, 0, "internal:'Indel'!A1", string="Indel variants")
worksheetOver.write_url(11, 0, "internal:'Low Coverage'!A1", string="Positions with coverage lower than " + str(min_cov) + "x")
worksheetOver.write_url(12, 0, "internal:'Hotspot'!A1", string="Coverage of hotspot positions")
worksheetOver.write_url(13, 0, "internal:'Coverage'!A1", string="Average coverage of all regions in bed")
worksheetOver.write_url(14, 0, "internal:'Version'!A1", string="Version Log")
worksheetOver.write_row(16, 0, emptyList, lineFormat)

worksheetOver.write_row(17, 0, ["Percent known variants found:"])
if num_known < len(known):
    worksheetOver.write_row(18, 0, [str(num_known / len(known) * 100) + " %"], redFormat)
else:
    worksheetOver.write_row(18, 0, [str(num_known / len(known) * 100) + " %"])

avgCov = extractMatchingLines("total_region", snakemake.input.mosdepth_summary, "-wE").split("\t")[3]
duplicateLevel = extractMatchingLines("PERCENT", snakemake.input.picard_dup, "-A1").split("\n")[-1].split("\t")[8]
with open(snakemake.input.mosdepth_thresh_summary) as threshold_summary:
    thresholds = threshold_summary.read().strip().split("\t")

worksheetOver.write_row(
    20,
    0,
    ["RunID", "DNAnr", "Avg. coverage [x]", "Duplicationlevel [%]", str(min_cov) + "x", str(med_cov) + "x", str(max_cov) + "x"],
    tableHeadFormat,
)
worksheetOver.write_row(
    21,
    0,
    [
        runid,
        sample,
        avgCov,
        str(round(float(duplicateLevel) * 100, 2)),
        str(thresholds[0]),
        str(thresholds[1]),
        str(thresholds[2]),
    ],
)

if low_pos == 0:  # From Hotspot sheet
    worksheetOver.write(24, 0, "Number of positions from the hotspot list not covered by at least " + str(med_cov) + "x: ")
    worksheetOver.write(25, 0, str(low_pos))
else:
    worksheetOver.write(24, 0, "Number of positions from the hotspot list not covered by at least " + str(med_cov) + "x: ")
    worksheetOver.write(25, 0, str(low_pos), redFormat)
    worksheetOver.write_url(26, 0, "internal:'Hotspot'!A1", string="For more detailed list see hotspotsheet ")


worksheetOver.write(27, 0, "Number of regions not covered by at least " + str(min_cov) + "x: ")  # From Cov sheet
worksheetOver.write(28, 0, str(low_regions))  # From Cov sheet

worksheetOver.write(30, 0, "Bedfile: " + snakemake.input.bedfile)
worksheetOver.write(31, 0, "Hotspotlist: " + snakemake.input.hotspot)
worksheetOver.write(32, 0, "Artefact file: " + snakemake.input.artefact_snv)
worksheetOver.write(33, 0, "Germline file: " + snakemake.input.germline)
worksheetOver.write(37, 0, "Bedfile for pindel: " + snakemake.input.bedfile_pindel)
worksheetOver.write(38, 0, "Pindel artefact file: " + snakemake.input.artefact_pindel)


# Known variants
worksheetKnown.set_column("E:E", 10)

worksheetKnown.write("A1", "Known variants ", headingFormat)
worksheetKnown.write("A3", "Sample: " + str(sample))
worksheetKnown.write("A4", "Reference used: " + str(refV))
worksheetKnown.write("A6", "VEP: " + vepline)  # , textwrapFormat)
worksheetKnown.write("A8", "The following filters were applied: ")
worksheetKnown.write("B9", "Coverage >= " + str(min_cov) + "x")
worksheetKnown.write("B10", "Population freq (KGP, gnomAD, NHLBI_ESP ) <= 2%")
worksheetKnown.write("B11", "Biotype is protein coding")
worksheetKnown.write("B12", "Consequence not deemed relevant")

worksheetKnown.write("A14", "For all variants see: ")
worksheetKnown.write_url("B14", "internal:'SNVs'!A1", string="SNVs")

worksheetKnown.write("A16", "Several variants found at same position.", italicFormat)
tableheading = [
    "RunID",
    "DNAnr",
    "Gene",
    "Variant",
    "Pos",
    "Type of Variant",
    "COSMIC ID",
    "Known AF",
    "Found AF",
    "DP",
    "Ref",
    "Alt",
    "SV length",
]
worksheetKnown.write_row("A18", tableheading, tableHeadFormat)  # 1 index
row = 18  # 0 index
col = 0

for line in known_lines_final:
    if len(line) < 9:
        worksheetKnown.write_row(row, col, line, redFormat)
    elif line[0] == "":
        worksheetKnown.write_row(row, col, line, italicFormat)
    else:
        worksheetKnown.write_row(row, col, line)
    row += 1

# SNV
worksheetSNV.set_column("E:E", 10)  # Width for position column

worksheetSNV.write("A1", "Variants found", headingFormat)
worksheetSNV.write("A3", "Sample: " + str(sample))
worksheetSNV.write("A4", "Reference used: " + str(refV))
worksheetSNV.write("A6", "VEP: " + vepline)
worksheetSNV.write("A8", "The following filters were applied: ")
worksheetSNV.write("B9", "Coverage >= " + str(min_cov) + "x")
worksheetSNV.write("B10", "Population freq (KGP, gnomAD, NHLBI_ESP ) <= 2%")
worksheetSNV.write("B11", "Biotype is protein coding")
worksheetSNV.write("B12", "Consequence not deemed relevant")

worksheetSNV.write("A14", "Coverage below " + str(med_cov) + "x", italicFormat)
worksheetSNV.write("A15", "Variant in artefact list ", orangeFormat)
worksheetSNV.write("A16", "Variant likely germline", greenFormat)
worksheetSNV.write("A17", "Variants with frequency 0.01 <= AF < 0.05 are located below artefact and germline variants.")

# Variant table
tableheading = [
    "RunID",
    "DNAnr",
    "Gene",
    "Chr",
    "Pos",
    "Ref",
    "Alt",
    "AF",
    "DP",
    "Transcript",
    "Mutation cds",
    "ENSP",
    "Consequence",
    "COSMIC ids on position",
    "N COSMIC Hemato hits on position",
    "Clinical significance",
    "dbSNP",
    "Max popAF",
    "Max Pop",
    "Callers",
]
worksheetSNV.write_row("A19", tableheading, tableHeadFormat)  # 1 index
row = 19  # 0 index
col = 0

for line in white:
    if line[8] < med_cov:
        worksheetSNV.write_row(row, col, line, italicFormat)
    else:
        worksheetSNV.write_row(row, col, line)
    row += 1

for line in green:
    if line[8] < med_cov:
        worksheetSNV.write_row(row, col, line, green_italicFormat)
    else:
        worksheetSNV.write_row(row, col, line, greenFormat)
    row += 1

for line in orange:
    if line[8] < med_cov:
        worksheetSNV.write_row(row, col, line, orange_italicFormat)
    else:
        worksheetSNV.write_row(row, col, line, orangeFormat)
    row += 1

for line in underFive:
    if line[8] < med_cov:
        worksheetSNV.write_row(row, col, line, italicFormat)
    else:
        worksheetSNV.write_row(row, col, line)
    row += 1

# (P)Indel sheet
col = 0
worksheetIndel.set_column("E:F", 10)  # pos
worksheetIndel.write("A1", "Pindel results", headingFormat)
worksheetIndel.write_row(1, col, emptyList, lineFormat)

for x in vcf_indel.header.records:
    if x.key == "reference":
        refI = x.value
worksheetIndel.write("A3", "Sample: " + str(sample))
worksheetIndel.write("A4", "Reference used: " + str(refI))
worksheetIndel.write("A5", "Genes included: ")
worksheetIndel.write("A6", "Genes included: ")
with open(snakemake.input.bedfile_pindel) as bed:
    genesDup = [line.split("\t")[3].strip() for line in bed]
    genes = set(genesDup)
row = 6
for gene in genes:
    worksheetIndel.write("B" + str(row), gene)
    row += 1
worksheetIndel.write(row, col, "Coverage below " + str(med_cov) + "x", italicFormat)
worksheetIndel.write(row + 1, col, "Variant in pindel artefact list.", orangeFormat)
worksheetIndel.write(row + 2, col, "Variants with frequency 0.01 <= AF < 0.05 are located below artefact variants.")
row += 5
tableheading = [
    "RunID",
    "DNAnr",
    "Gene",
    "Chr",
    "Start",
    "End",
    "SV length",
    "Af",
    "Ref",
    "Alt",
    "Dp",
    "Transcript",
    "Mutation cds",
    "ENSP",
    "Max popAF",
    "Max Pop",
]
worksheetIndel.write_row("A" + str(row), tableheading, tableHeadFormat)  # 1 index

for line in whiteIndel:
    if line[10] < med_cov:
        worksheetIndel.write_row(row, col, line, italicFormat)
    else:
        worksheetIndel.write_row(row, col, line)
    row += 1

for line in orangeIndel:
    if line[10] < med_cov:
        worksheetIndel.write_row(row, col, line, orange_italicFormat)
    else:
        worksheetIndel.write_row(row, col, line, orangeFormat)
    row += 1

for line in underFiveIndel:
    if line[10] < med_cov:
        worksheetIndel.write_row(row, col, line, italicFormat)
    else:
        worksheetIndel.write_row(row, col, line)
    row += 1

# Low cvergage sheet
worksheetLowCov.set_column(1, 3, 10)
worksheetLowCov.set_column(1, 4, 10)
worksheetLowCov.write("A1", "Mosdepth coverage analysis", headingFormat)
worksheetLowCov.write_row("A2", emptyList, lineFormat)
worksheetLowCov.write("A3", "Sample: " + str(sample))
description = "Gene Regions with coverage lower than " + str(min_cov) + "x."
worksheetLowCov.write("A4", description)
covHeadings = ["Region Name", "Chr", "Start", "Stop", "Mean Coverage", "Length of Region"]
worksheetLowCov.write_row("A6", covHeadings, tableHeadFormat)  # 1 index
row = 6  # 0 index
for line in low_cov_lines:
    worksheetLowCov.write_row(row, col, line)
    row += 1

# Hotspot
# should we have gene as well not just pos?
worksheetHotspot.set_column(1, 2, 10)
worksheetHotspot.write("A1", "Hotspot Coverage", headingFormat)
worksheetHotspot.write("A3", "Sample: " + str(sample))
worksheetHotspot.write_row("A5", ["Gene", "Chr", "Pos", "Depth"], tableHeadFormat)
row = 5
for hotLine in hotspotTable:
    if int(hotLine[2]) <= med_cov:
        worksheetHotspot.write_row(row, 0, hotLine, redFormat)
    else:
        worksheetHotspot.write_row(row, 0, hotLine)
    row += 1

# Coverage
worksheetCov.write("A1", "Average Coverage", headingFormat)
worksheetCov.write_row("A2", emptyList, lineFormat)
worksheetCov.write("A3", "Sample: " + str(sample))
worksheetCov.write("A4", "Averge coverage of each region in bedfile from Mosdepth")

tableArea = "A6:F" + str(len(cov_table_lines) + 6)  # rows of full list
headerListDict = [
    {"header": "Region Name"},
    {"header": "Chr"},
    {"header": "Start"},
    {"header": "End"},
    {"header": "Avg Coverage"},
]
worksheetCov.add_table(tableArea, {"data": cov_table_lines, "columns": headerListDict, "style": "Table Style Light 1"})

# QCI
worksheetQCI.set_column("C:C", 10)
worksheetQCI.write("A1", "Results from QCI ", headingFormat)
worksheetQCI.write_row("A2", emptyList, lineFormat)

worksheetQCI.write("A5", "Analysen utfördes i enlighet med dokumentationen.")
worksheetQCI.write("A6", "Eventuella avikelser:")
iva = [
    "DNA nr",
    "Chromosome",
    "Position",
    "Gene Region",
    "Gene Symbol",
    "Transcript ID",
    "Transcript Variant",
    "Protein Variant",
    "Variant Findings",
    "Sample Genotype Quality",
    "Read Depth",
    "Allele Fraction",
    "Translation Impact",
    "dbSNP ID",
    "1000 Genomes Frequency",
    "ExAC Frequency",
    "HGMD",
    "COSMIC ID",
    "Artefacts_without_ASXL1",
    "ASXL1_variant_filter",
]
worksheetQCI.write_row(9, 0, iva, tableHeadFormat)

# Versions
worksheetVersions.write("A1", "Version Log", headingFormat)
worksheetVersions.write_row(1, 0, emptyList, lineFormat)
worksheetVersions.write("A3", "Sample: " + str(sample))
worksheetVersions.write("A5", "Variant calling reference used: " + str(refV))
worksheetVersions.write("A6", "Pindel reference used: " + str(refI))
worksheetVersions.write("A7", "Containers used: ", tableHeadFormat)
containers = [clist for clist in snakemake.params.singularitys.items()]
row = 8
col = 0
for containerTuple in containers:
    container = list(containerTuple)
    worksheetVersions.write_row("A" + str(row), container)
    row += 1

workbook.close()

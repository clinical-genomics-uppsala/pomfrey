#!/bin/python3.6
import sys
import csv
from pysam import VariantFile
import xlsxwriter
from datetime import date
import subprocess
import yaml
import gzip

# Define sys.argvs
vcf_snv = VariantFile(snakemake.input.vcf_snv)
vcf_indel = VariantFile(snakemake.input.vcf_pindel)

runid = snakemake.params.seqid
sample = list(vcf_snv.header.samples)[0]
today = date.today()
emptyList = ['', '', '', '', '', '']

sample_purity = 0.8


# VEP fileds in list to get index
def index_vep(variantfile):
    csqIndex = []
    for x in variantfile.header.records:
        if 'CSQ' in str(x):
            csqIndex = str(x).split('Format: ')[1].strip().strip('">').split('|')
    return csqIndex


# Return matching lines in file
def extractMatchingLines(expressionMatch, artefactFile, grepvarible):
    if grepvarible == '':
        grepvarible = '-wE '
    cmdArt = 'grep '+grepvarible+' '+str(expressionMatch)+' '+artefactFile
    matchLines = subprocess.run(cmdArt, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
    return matchLines


def file_length(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def float_or_na(value):
    return float(value) if value != '' else None


def int_or_na(value):
    return int(value) if value != '' else None


shortListGenes = ['ABL1', 'ANKRD26', 'ASXL1', 'ATRX', 'BCOR', 'BCORL1', 'BRAF', 'CALR', 'CBL', 'CBLB', 'CDKN2A', 'CEBPA', 'CSF3R',
                  'CUX1', 'DDX41', 'DNMT3A', 'ETV6', 'ETNK1', 'TEL', 'EZH2', 'FBXW7', 'FLT3', 'GATA1', 'GATA2', 'GNAS', 'HRAS',
                  'IDH1', 'IDH2', 'IKZF1', 'JAK2', 'JAK3', 'KDM6A', 'KIT', 'KRAS', 'KMT2A', 'MPL', 'MYD88', 'NF1', 'NOTCH1',
                  'NPM1', 'NRAS', 'PDGFRA', 'PHF6', 'PPM1D', 'PTEN', 'PTPN11', 'RAD21', 'RUNX1', 'SAMD9', 'SAMD9L', 'SETBP1',
                  'SF3B1', 'SMC1A', 'SMC3', 'SRSF2', 'STAG2', 'STAT3', 'STAT5B', 'TET2', 'TP53', 'U2AF1', 'WT1', 'ZRSR2']


intronDict = {'GATA2': ['chr3', 128201827,  128202419],
              'TERC': ['chr3', 169482182, 169483654],
              'NOTCH1': ['chr9', 139388885, 139390523],
              'ANKRD26': ['chr10', 27389007, 27389433],
              'TP53': ['chr17', 7590690, 7590874]}

introns = {}
for key in intronDict:
    chr = intronDict[key][0]
    if chr in introns:  # If chr in dict already
        introns[chr].append(intronDict[key][1:])
    else:
        introns[chr] = [intronDict[key][1:]]

synoVariants = [['c.1416G>A', 'chr3', '128199889', 'C', 'T'], ['c.1023C>T', 'chr3', '128200782', 'G', 'A'],
                ['c.981G>A', 'chr3', '128202739', 'C', 'T'], ['c.649C>T', 'chr3', '128204792', 'G', 'A'],
                ['c.351C>G', 'chr3', '128205090', 'G', 'C'], ['c.375G>A', 'chr17', '7579312', 'C', 'T'],
                ['c.375G>T', 'chr17', '7579312', 'C', 'A'], ['c.375G>C', 'chr17', '7579312', 'C', 'G'],
                ['c.672G>A', 'chr17', '7578177', 'C', 'T'], ['c.993G>A', 'chr17', '7576853', 'C', 'T']]

min_cov = int(snakemake.params.thresholds.split(',')[0])
med_cov = int(snakemake.params.thresholds.split(',')[1])
max_cov = int(snakemake.params.thresholds.split(',')[2])

''' Coverage of all regions in bedfile '''
bed_table = []
with open(snakemake.input.bedfile, 'r') as bedfile:
    next(bedfile)
    # Skip header
    for line in bedfile:
        bed_table.append(line.strip().split('\t'))
# Coverage per region for Coverage sheet
cov_table_lines = []
with gzip.open(snakemake.input.mosdepth_regions, 'rt') as regionsfile:
    for lline in regionsfile:
        line = lline.strip().split('\t')
        cov_row = [line[3], line[0], line[1], line[2], line[4]]
        cov_table_lines.append(cov_row)

''' Low cov file '''
low_cov_lines = []
condensed_line = ['', '', '', '', '']
with open(snakemake.input.mosdepth_lowcov, 'r') as lowfile:
    for lline in lowfile:
        line = lline.strip().split("\t")
        if int(line[3]) <= int(min_cov):
            if condensed_line[0] == line[0] and condensed_line[2] == line[1]:
                condensed_length = int(condensed_line[2]) - int(condensed_line[1])
                line_length = int(line[2]) - int(line[1])
                new_avgcov = (float(condensed_line[3])*condensed_length+float(line[3])*line_length)/(line_length+condensed_length)
                condensed_line[3] = new_avgcov
                condensed_line[2] = line[2]
            else:
                low_cov_lines.append(condensed_line)
                condensed_line = line

low_cov_lines.pop(0)
for line in low_cov_lines:
    for bedline in bed_table:  # blir det en rad eller element?
        if line[0] == bedline[0] and int(line[1]) >= int(bedline[1]) and int(line[2]) <= int(bedline[2]):
            line[:] = [bedline[3]]+line[:]+[int(line[2])-int(line[1])]
            break

# Sort based on coverage
low_cov_lines.sort(key=lambda x: float(x[4]))

# Number of low cov regions for the Overview sheet.
low_regions = len(low_cov_lines)

''' Hotspot file '''
low_pos = 0
hotspotTable = []
with open(snakemake.input.mosdepth_hotspot, 'r') as hotFile:
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


''' Loop through vcf_snv to get snvs and intron '''
# SNV
white = []
green = []
orange = []
whiteIGV = []
underFive = []  # put after green and orange but still white
underFiveIGV = []  # put after green and orange but still white
# Reported
shortListSNV = []  # Reported genes only
shortListSNVigv = []  # Reported genes only
greenShortList = []  # Germline in reported genes
# Intron
intron_variants = []
synoFound = []

for x in vcf_snv.header.records:
    if (x.key == 'reference'):
        refV = x.value
    if (x.key == 'VEP'):
        vepline = x.value

for record in vcf_snv.fetch():
    # Get info from vep annotation
    csqIndex = index_vep(vcf_snv)
    csq = record.info["CSQ"][0].split("|")
    consequence = csq[csqIndex.index('Consequence')]
    gene = csq[csqIndex.index('SYMBOL')]
    clinical = csq[csqIndex.index('CLIN_SIG')]
    existing = csq[csqIndex.index('Existing_variation')].split('&')  # Vad hander om ej existerar?

    if len([rs for rs in existing if rs.startswith('rs')]) == 0:
        rs = ''
    else:
        rs = ', '.join([rs for rs in existing if rs.startswith('rs')])

    try:
        if record.info["CALLERS"]:
            callers = ' & '.join(record.info["CALLERS"])
    except KeyError:
        callers = 'Pisces-multi'

    transcript = csq[csqIndex.index('HGVSc')].split(":")[0]
    if len(csq[csqIndex.index('HGVSc')].split(":")) > 1:
        codingName = csq[csqIndex.index('HGVSc')].split(":")[1]
    else:
        codingName = ''
    ensp = csq[csqIndex.index('HGVSp')]

    popFreqsPop = csqIndex[csqIndex.index('AF'):csqIndex.index('gnomAD_SAS_AF')+1]
    popFreqAllRaw = csq[csqIndex.index('AF'):csqIndex.index('gnomAD_SAS_AF')+1]
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
        maxPopAf = ''
        maxPop = ''

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

    # SNV variant
    synoCosmicN = 0
    spliceVariant = False
    if record.filter.keys() == ["Syno"]:  # Only if Syno not and popAF.   any(x in "Syno" for x in record.filter.keys()):
        synoCosmicVepList = [cosmic for cosmic in csq[csqIndex.index('Existing_variation')].split("&")
                             if cosmic.startswith('CO')]  # Get all cosmicID in list
        # COSMIC Hemato
        if len(synoCosmicVepList) != 0:
            for synoCosmicId in synoCosmicVepList:
                try:
                    synoCosmicNew = extractMatchingLines(synoCosmicId, snakemake.input.hemato_count, '-wE').split("\t")[15]
                except IndexError:
                    synoCosmicNew = 0
                synoCosmicN += int(synoCosmicNew)
        if 'splice' in consequence:
            spliceVariant = True

    if (record.filter.keys() == ["PASS"] or synoCosmicN != 0 or spliceVariant) and af >= 0.01:
        # Total number of cosmic hemato hits on the position. Vep reports all cosmicId for that position.
        cosmicVepList = [cosmic for cosmic in existing if cosmic.startswith('CO')]
        if len(cosmicVepList) == 0:
            cosmicVep = ''
        else:
            cosmicVep = ', '.join(cosmicVepList)

        if len(cosmicVepList) == 0:
            cosmicN = ''
        else:
            cosmicN = 0
            for cosmicId in cosmicVepList:
                try:
                    cosmicNew = extractMatchingLines(cosmicId, snakemake.input.hemato_count, '-wE').split("\t")[15]
                except IndexError:
                    cosmicNew = 0
                cosmicN += int(cosmicNew)

        # IGV image path for each SNV
        igv = "external:IGV/"+gene+"-"+record.contig+"_"+str(int(record.pos)-1)+"_"+str(int(record.pos)-1+len(alt))+".svg"

        snv = [runid, sample, gene, record.contig, record.pos, record.ref, alt, af, record.info["DP"], transcript,
               codingName, ensp, consequence, cosmicVep, cosmicN, clinical, rs, maxPopAf, maxPop, '', callers]

        # Append line with sample and rundate to rolling list of artefacts..
        with open(snakemake.input.variantslog, "a") as appendfile:
            variants = snv+["\n"]
            appendfile.write('\t'.join(str(e) for e in variants))

        # Check if variant in artefact or germline file
        artLines = extractMatchingLines(str(record.pos), snakemake.input.artefact_snv, '-wE')
        artefact_variant = False

        for artLine in artLines.split("\n"):
            # if pos exists and match in artefact file.
            if artLine and record.ref == artLine.split()[2] and alt == artLine.split()[3]:
                orange.append(snv)
                artefact_variant = True
                break
        if not artefact_variant:
            germLines = extractMatchingLines(str(record.pos), snakemake.input.germline, '-wE')
            germline_variant = False
            for germLine in germLines.split("\n"):
                # if exists in germline file
                if germLine and record.ref == germLine.split()[2] and alt == germLine.split()[3]:
                    green.append(snv)
                    germline_variant = True
                    if gene in shortListGenes:
                        greenShortList.append(snv)
                    break
            if not germline_variant:
                if float(af) < 0.05:
                    underFive.append(snv)
                    underFiveIGV.append(igv)
                else:
                    white.append(snv)
                    whiteIGV.append(igv)
                if gene in shortListGenes:
                    shortListSNV.append(snv)
                    shortListSNVigv.append(igv)

    # Intron variants with allel freq over 20 % in specified regions (intronDict)
    if "PopAF" not in record.filter.keys() and record.contig in introns:
        for pair in introns[record.contig]:
            if record.pos >= pair[0] and record.pos <= pair[1] and af >= 0.2 and int(record.info["DP"]) >= 100:
                intron_line = [runid, sample, gene, record.contig, str(record.pos), record.ref, alt, af, str(record.info["DP"]),
                               transcript, codingName, ensp, consequence, maxPopAf, maxPop, callers]
                intron_variants.append(intron_line)

    # Synonymous variants in GATA2 and TP53 (synoVariants)
    for synoVariant in synoVariants:
        if record.contig == synoVariant[1] and record.pos == int(synoVariant[2]) and alt == synoVariant[4]:
            syno_line = [runid, sample, gene, record.contig, str(record.pos), record.ref, alt, af, str(record.info["DP"]),
                         transcript, codingName, ensp, consequence, maxPopAf, maxPop, callers]
            synoFound.append(syno_line)


''' Loop through pindel vcf '''
orangeIndel = []
whiteIndel = []
whiteIGVIndel = []
underFiveIndel = []
underFiveIGVIndel = []
csqIndex = index_vep(vcf_indel)

for indel in vcf_indel.fetch():
    # Borde man ta med alla och istallet lagga till en filterkolumn? Hur blir det med icke proteincoding och
    # konsekvens som kanske inte blir samma sak.
    if indel.filter.keys() == ["PASS"] and float(indel.info["AF"]) >= 0.01:
        svlen = indel.info["SVLEN"]
        af = indel.info["AF"]

        if len(indel.alts) == 1:
            alt = indel.alts[0]
        else:
            print(indel.alts)
            sys.exit()

        csqIndel = indel.info["CSQ"][0].split("|")  # VEP annotation
        indelGene = csqIndel[csqIndex.index('SYMBOL')]

        # Not using ExAC pop
        popFreqsPop = csqIndex[csqIndex.index('AF'):csqIndex.index('gnomAD_SAS_AF')+1]
        popFreqAllRawIndel = csqIndel[csqIndex.index('AF'):csqIndex.index('gnomAD_SAS_AF')+1]
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
            maxPopAfIndel = ''
            maxPopIndel = ''

        indelTranscript = csqIndel[csqIndex.index('HGVSc')].split(":")[0]
        if len(csqIndel[csqIndex.index('HGVSc')].split(":")) > 1:
            indelCodingName = csqIndel[csqIndex.index('HGVSc')].split(":")[1]
        else:
            indelCodingName = ''
        indelEnsp = csqIndel[csqIndex.index('HGVSp')]

        indelRow = [runid, sample, indelGene, indel.contig, indel.pos, indel.stop, svlen, af, indel.ref,
                    alt, indel.info["DP"], indelTranscript, indelCodingName, indelEnsp, maxPopAfIndel, maxPopIndel]

        # Mark artefact based on artefactfile
        artLines = extractMatchingLines(str(indel.contig)+'.*'+str(indel.pos), snakemake.input.artefact_pindel, '-wE')

        if len(artLines) > 0:  # if pos exists and match in artefact file.
            orangeIndel.append(indelRow)
        else:
            indelIgv = "external:IGV/"+indelGene+"-"+indel.contig+"_" + \
                str(int(indel.pos)-1)+"_"+str(int(indel.pos)-1+len(alt))+".svg"
            if float(af) < 0.05:
                underFiveIndel.append(indelRow)
                underFiveIGVIndel.append(indelIgv)
            else:
                whiteIndel.append(indelRow)
                whiteIGVIndel.append(indelIgv)


''' Loop through GATKs CNV vcf '''
# Import genomic pos to cytoCoord translation to list
chrBands = []
with open(snakemake.input.cyto_coord_convert, 'r') as chrBandFile:
    for line in chrBandFile:
        chrBands.append(line.split("\t"))

# Load in bedfile with generegions condensed
gene_regions = {}
with open(snakemake.input.bedfile_cnv, 'r') as cnv_bed_file:
    for line in cnv_bed_file:
        lline = line.strip().split("\t")
        chrom = lline[0]
        start = lline[1]
        end = lline[2]
        name = lline[3]
        gene = name.split("_")[0]

        if gene in gene_regions:  # gene lengths
            gene_regions[gene][2] = end
        else:
            gene_regions[gene] = [chrom, start, end]

# Process GATK4 CNV seg-file
cnv_lines = []
with open(snakemake.input.gatk_seg, 'r') as GATK_file:
    header = True
    for line in GATK_file:  # Skip ahead until data
        genes = []
        if header:
            if line[0:3] == 'chr':
                header = False
            else:
                continue
        lline = line.strip().split("\t")
        start_pos = int(lline[1])
        end_pos = int(lline[2])
        chrom = lline[0]
        # Translate genomic coordinate to cytogen coordinates
        cytoCoord = ['', '']
        for chrBand in chrBands:
            if chrBand[0] == chrom:
                if (start_pos >= int(chrBand[1]) and start_pos <= int(chrBand[2])):
                    cytoCoord[0] = chrBand[3]
                if (end_pos >= int(chrBand[1]) and end_pos <= int(chrBand[2])):
                    cytoCoord[1] = chrBand[3]
        if cytoCoord[0] == cytoCoord[1]:
            cytoCoordString = chrom[3:]+cytoCoord[0]
        else:
            cytoCoordString = chrom[3:]+cytoCoord[0]+'-'+cytoCoord[1]
        # Only look at cnv that gatk marked as not neutral
        if lline[5] != '0':
            for gene, coordinates in gene_regions.items():
                if coordinates[0] == chrom:
                    if (start_pos >= int(coordinates[1]) and start_pos <= int(coordinates[2])) or (end_pos >= int(coordinates[1])
                       and end_pos <= int(coordinates[2])) or (start_pos < int(coordinates[1]) and end_pos > int(coordinates[2])):
                        genes.append(gene)
            geneString = ', '.join([x for x in genes if not x.startswith('CNV')])
            logRatio = float(lline[4])
            copyNumberTumor = round(2*pow(2, logRatio), 2)

            cnv_lines.append([sample, geneString, chrom, str(start_pos)+'-'+str(end_pos), cytoCoordString, str(sample_purity),
                             '', '', str(round(logRatio, 4)), str(copyNumberTumor)])

# Process CNVkit cns file

chromosomes = ['chr'+str(i) for i in range(1, 23)]+['chrX', 'chrY']
relevant_cnvs = {i: [] for i in chromosomes}
relevant_cnvs_header = ['Sample', 'Chromosome', 'Start', 'End', 'CytoCoordinates', 'Log2',
                        'CI high', 'CI low', 'BAF', 'Copy Number',
                        'Copies Allele 1', 'Copies Allele 2', 'Depth', 'Probes', 'Weight', 'Genes']
with open(snakemake.input.cnvkit_calls, 'r+') as cnsfile:
    cns_header = next(cnsfile).rstrip().split("\t")
    for cnv_line in cnsfile:
        cnv = cnv_line.strip().split("\t")
        if not (cnv[cns_header.index('cn')] == '2' and cnv[cns_header.index('cn1')] == '1'):
            cnv_chr = cnv[cns_header.index('chromosome')]
            cnv_start = int(cnv[cns_header.index('start')])
            cnv_end = int(cnv[cns_header.index('end')])
    #        import pdb; pdb.set_trace()
            cnv_baf = float_or_na(cnv[cns_header.index('baf')])
            cytoCoord = ['', '']
            for chrBand in chrBands:
                if chrBand[0] == cnv_chr:
                    if (cnv_start >= int(chrBand[1]) and cnv_start <= int(chrBand[2])):
                        cytoCoord[0] = chrBand[3]
                    if (cnv_end >= int(chrBand[1]) and cnv_end <= int(chrBand[2])):
                        cytoCoord[1] = chrBand[3]
            if cytoCoord[0] == cytoCoord[1]:
                cytoCoordString = cnv_chr[3:]+cytoCoord[0]
            else:
                cytoCoordString = cnv_chr[3:]+cytoCoord[0]+'-'+cytoCoord[1]
            outline = [sample, cnv_chr, cnv_start, cnv_end, cytoCoordString, float(cnv[cns_header.index('log2')]),
                       float(cnv[cns_header.index('ci_hi')]), float(cnv[cns_header.index('ci_lo')]), cnv_baf,
                       cnv[cns_header.index('cn')], int_or_na(cnv[cns_header.index('cn1')]),
                       int_or_na(cnv[cns_header.index('cn2')]), cnv[cns_header.index('depth')],
                       cnv[cns_header.index('probes')], cnv[cns_header.index('weight')], str(cnv[cns_header.index('gene')])]
            relevant_cnvs[cnv_chr].append(outline)


''' Xlsx sheets '''
workbook = xlsxwriter.Workbook(snakemake.output[0])
worksheetOver = workbook.add_worksheet('Overview')
worksheetShortList = workbook.add_worksheet('ShortList')
worksheetSNV = workbook.add_worksheet('SNVs')
worksheetIndel = workbook.add_worksheet('InDel')
worksheetIntron = workbook.add_worksheet('Intron & Synonymous')
worksheetCNV = workbook.add_worksheet('CNV GATK')
worksheetCNVkit = workbook.add_worksheet('CNVkit')
worksheetLowCov = workbook.add_worksheet('Low Coverage')
worksheetHotspot = workbook.add_worksheet('Hotspot')
worksheetCov = workbook.add_worksheet('Coverage')
worksheetQCI = workbook.add_worksheet('QCI')
worksheetVersions = workbook.add_worksheet('Version')
# Define formats to be used.
headingFormat = workbook.add_format({'bold': True, 'font_size': 18})
lineFormat = workbook.add_format({'top': 1})
tableHeadFormat = workbook.add_format({'bold': True, 'text_wrap': True})
textwrapFormat = workbook.add_format({'text_wrap': True})
italicFormat = workbook.add_format({'italic': True})
redFormat = workbook.add_format({'font_color': 'red'})

greenFormat = workbook.add_format({'bg_color': '#85e085'})
orangeFormat = workbook.add_format({'bg_color': '#ffd280'})
green_italicFormat = workbook.add_format({'bg_color': '#85e085', 'italic': 'True'})
orange_italicFormat = workbook.add_format({'bg_color': '#ffd280', 'italic': 'True'})

# Overview
worksheetOver.write(0, 0, sample, headingFormat)
worksheetOver.write(1, 0, "RunID: "+runid)
worksheetOver.write(2, 0, "Processing date: "+today.strftime("%B %d, %Y"))
worksheetOver.write_row(3, 0, emptyList, lineFormat)

worksheetOver.write(4, 0, "Created by: ")
worksheetOver.write(4, 4, "Valid from: ")
worksheetOver.write(5, 0, "Signed by: ")
worksheetOver.write(5, 4, "Document nr: ")
worksheetOver.write_row(6, 0, emptyList, lineFormat)

worksheetOver.write(7, 0, "Sheets:", tableHeadFormat)
worksheetOver.write_url(8, 0, "internal:'ShortList'!A1", string='Variants in genes to report')
worksheetOver.write_url(9, 0, "internal:'SNVs'!A1", string='Variants analysis')
worksheetOver.write_url(10, 0, "internal:'Indel'!A1", string='Indel variants')
worksheetOver.write_url(11, 0, "internal:'Intron & Synonymous'!A1", string='Intron & synonymous variants')
worksheetOver.write_url(12, 0, "internal:'CNV GATK'!A1", string='CNVs found with GATK4')
worksheetOver.write_url(13, 0, "internal:'CNVkit'!A1", string='CNVs found with CNVkit')
worksheetOver.write_url(14, 0, "internal:'Low Coverage'!A1", string='Positions with coverage lower than '+str(min_cov)+'x')
worksheetOver.write_url(15, 0, "internal:'Hotspot'!A1", string='Coverage of hotspot positions')
worksheetOver.write_url(16, 0, "internal:'Coverage'!A1", string='Average coverage of all regions in bed')
worksheetOver.write_url(17, 0, "internal:'Version'!A1", string='Version Log')
worksheetOver.write_row(18, 0, emptyList, lineFormat)

avgCov = extractMatchingLines("total_region", snakemake.input.mosdepth_summary, '-wE').split('\t')[3]
duplicateLevel = extractMatchingLines('PERCENT', snakemake.input.picard_dup, '-A1').split('\n')[-1].split('\t')[8]
with open(snakemake.input.mosdepth_thresh_summary) as threshold_summary:
    thresholds = threshold_summary.read().strip().split("\t")

worksheetOver.write_row(20, 0, ['RunID', 'DNAnr', 'Avg. coverage [x]', 'Duplicationlevel [%]',
                                str(min_cov)+'x', str(med_cov)+'x', str(max_cov)+'x'], tableHeadFormat)
worksheetOver.write_row(21, 0, [runid, sample, avgCov, str(round(float(duplicateLevel)*100, 2)),
                                str(thresholds[0]), str(thresholds[1]), str(thresholds[2])])

if low_pos == 0:  # From Hotspot sheet
    worksheetOver.write(24, 0, 'Number of positions from the hotspot list not covered by at least '+str(med_cov)+'x: ')
    worksheetOver.write(25, 0, str(low_pos))
else:
    worksheetOver.write(24, 0, 'Number of positions from the hotspot list not covered by at least '+str(med_cov)+'x: ')
    worksheetOver.write(25, 0, str(low_pos), redFormat)
    worksheetOver.write_url(26, 0, "internal:'Hotspot'!A1", string='For more detailed list see hotspotsheet ')

worksheetOver.write(27, 0, 'Number of regions not covered by at least '+str(min_cov)+'x: ')
worksheetOver.write(28, 0, str(low_regions))  # From low cov sheet

cov_chrX = extractMatchingLines('chrX_region', snakemake.input.mosdepth_summary, '-wE').split('\t')[3]
cov_chrY = extractMatchingLines('chrY_region', snakemake.input.mosdepth_summary, '-wE').split('\t')[3]

worksheetOver.write(30, 0, 'Average coverage of region in bedfile:', tableHeadFormat)
worksheetOver.write_row(31, 0, ['chrX', cov_chrX])
worksheetOver.write_row(32, 0, ['chrY', cov_chrY])


worksheetOver.write(34, 0, 'Bedfile: ' + snakemake.input.bedfile)
worksheetOver.write(35, 0, 'Hotspotlist: ' + snakemake.input.hotspot)
worksheetOver.write(36, 0, 'Artefact file: ' + snakemake.input.artefact_snv)
worksheetOver.write(37, 0, 'Germline file: ' + snakemake.input.germline)
worksheetOver.write(38, 0, 'Bedfile for pindel: ' + snakemake.input.bedfile_pindel)
worksheetOver.write(39, 0, 'Pindel artefact file: ' + snakemake.input.artefact_pindel)
worksheetOver.write(40, 0, 'CNVkit artefact file: ' + snakemake.input.cnvkit_artefact)

# Reported variants
tableheading = ['RunID', 'DNAnr', 'Gene', 'Chr', 'Pos', 'Ref', 'Alt', 'AF', 'DP', 'Transcript', 'Mutation cds', 'ENSP',
                'Consequence', 'COSMIC ids on position', 'N COSMIC Hemato hits on position', 'Clinical significance', 'dbSNP',
                'Max popAF', 'Max Pop', 'IGV image', 'Callers']
worksheetShortList.set_column('E:E', 10)
worksheetShortList.write('A1', 'Variants in genes to report', headingFormat)
worksheetShortList.write('A3', 'Sample: ' + str(sample))
worksheetShortList.write('A6', 'VEP: ' + vepline)  # , textwrapFormat)
worksheetShortList.write('A8', 'The following filters were applied: ')
worksheetShortList.write('B9', 'Coverage >= ' + str(min_cov) + 'x')
worksheetShortList.write('B10', 'Population freq (KGP, gnomAD, NHLBI_ESP ) <= 2%')
worksheetShortList.write('B11', 'Biotype is protein coding')
worksheetShortList.write('B12', 'Consequence not deemed relevant')


worksheetShortList.write('A14', 'Short list of variants in genes to report:')
worksheetShortList.write_row(14, 0, shortListGenes)
worksheetShortList.write('A16', 'For all variants see: ')
worksheetShortList.write_url('B16', "internal:'SNVs'!A1", string='SNVs')

worksheetShortList.write('A18', 'Coverage below '+str(med_cov)+'x', italicFormat)
worksheetShortList.write('A19', 'Variant likely germline', greenFormat)
worksheetShortList.write_row('A21', tableheading, tableHeadFormat)  # 1 index
row = 21  # 0 index
col = 0
i = 0
for line in shortListSNV:
    if line[8] < med_cov:
        worksheetShortList.write_row(row, col, line, italicFormat)
        worksheetShortList.write_url('T'+str(row+1), shortListSNVigv[i], string="IGV image")
    else:
        worksheetShortList.write_row(row, col, line)
        worksheetShortList.write_url('T'+str(row+1), shortListSNVigv[i], string="IGV image")
    row += 1
    i += 1
for line in greenShortList:
    if line[8] < med_cov:
        worksheetShortList.write_row(row, col, line, green_italicFormat)
    else:
        worksheetShortList.write_row(row, col, line, greenFormat)
    row += 1


# SNV sheet
worksheetSNV.set_column('E:E', 10)  # Width for position column
worksheetSNV.write('A1', 'Variants found', headingFormat)
worksheetSNV.write('A3', 'Sample: '+str(sample))
worksheetSNV.write('A4', 'Reference used: '+str(refV))
worksheetSNV.write('A6', 'VEP: '+vepline)  # , textwrapFormat)
worksheetSNV.write('A8', 'The following filters were applied: ')
worksheetSNV.write('B9', 'Coverage >= '+str(min_cov)+'x')
worksheetSNV.write('B10', 'Population freq (KGP, gnomAD, NHLBI_ESP ) <= 2%')
worksheetSNV.write('B11', 'Biotype is protein coding')
worksheetSNV.write('B12', 'Consequence not deemed relevant')

worksheetSNV.write('A14', 'Coverage below '+str(med_cov)+'x', italicFormat)
worksheetSNV.write('A15', 'Variant in artefact list ', orangeFormat)
worksheetSNV.write('A16', 'Variant likely germline', greenFormat)
worksheetSNV.write('A17', 'Variants with frequency 0.01 <= AF < 0.05 are located below artefact and germline variants.')
worksheetSNV.write_row('A19', tableheading, tableHeadFormat)  # 1 index
row = 19  # 0 index
col = 0
i = 0
for line in white:
    if line[8] < med_cov:
        worksheetSNV.write_row(row, col, line, italicFormat)
        worksheetSNV.write_url('T'+str(row+1), whiteIGV[i], string="IGV image")
    else:
        worksheetSNV.write_row(row, col, line)
        worksheetSNV.write_url('T'+str(row+1), whiteIGV[i], string="IGV image")
    row += 1
    i += 1
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
i = 0
for line in underFive:
    if line[8] < med_cov:
        worksheetSNV.write_row(row, col, line, italicFormat)
        worksheetSNV.write_url('T'+str(row+1), underFiveIGV[i], string="IGV image")
    else:
        worksheetSNV.write_row(row, col, line)
        worksheetSNV.write_url('T'+str(row+1), underFiveIGV[i], string="IGV image")
    row += 1
    i += 1

# (P)Indel sheet
worksheetIndel.set_column('E:F', 10)  # pos
worksheetIndel.write('A1', 'Pindel results', headingFormat)
worksheetIndel.write_row(1, 0, emptyList, lineFormat)

for x in vcf_indel.header.records:
    if (x.key == 'reference'):
        refI = x.value
worksheetIndel.write('A3', 'Sample: '+str(sample))
worksheetIndel.write('A4', 'Reference used: '+str(refI))
worksheetIndel.write('A6', 'Genes included: ')
with open(snakemake.input.bedfile_pindel) as bed:
    genesDup = [line.split("\t")[3].strip() for line in bed]
    genes = set(genesDup)
row = 6
for gene in genes:
    worksheetIndel.write('B'+str(row), gene)
    row += 1
worksheetIndel.write(row, 0, 'Coverage below '+str(med_cov)+'x', italicFormat)
worksheetIndel.write(row+1, 0, 'Variant in artefact list ', orangeFormat)
worksheetIndel.write(row+2, 0, 'Variants with frequency 0.01 <= AF < 0.05 are located below artefact and germline variants.')
row += 5
tableheading = ['RunID', 'DNAnr', 'Gene', 'Chr', 'Start', 'End', 'SV length', 'Af', 'Ref',
                'Alt', 'Dp', 'Transcript', 'Mutation cds', 'ENSP', 'Max popAF', 'Max Pop', 'IGV']
worksheetIndel.write_row('A'+str(row), tableheading, tableHeadFormat)  # 1 index

col = 0
i = 0
for line in whiteIndel:
    if line[10] < med_cov:
        worksheetIndel.write_row(row, col, line, italicFormat)
        worksheetIndel.write_url('Q'+str(row+1), whiteIGVIndel[i], string="IGV image")
    else:
        worksheetIndel.write_row(row, col, line)
        worksheetIndel.write_url('Q'+str(row+1), whiteIGVIndel[i], string="IGV image")
    row += 1
    i += 1

for line in orangeIndel:
    if line[10] < med_cov:
        worksheetIndel.write_row(row, col, line, orange_italicFormat)
    else:
        worksheetIndel.write_row(row, col, line, orangeFormat)
    row += 1

i = 0
for line in underFiveIndel:
    if line[10] < med_cov:
        worksheetIndel.write_row(row, col, line, italicFormat)
        worksheetIndel.write_url('Q'+str(row+1), underFiveIGVIndel[i], string="IGV image")
    else:
        worksheetIndel.write_row(row, col, line)
        worksheetIndel.write_url('Q'+str(row+1), underFiveIGVIndel[i], string="IGV image")
    row += 1
    i += 1


# Intron and synonymous
worksheetIntron.set_column('D:E', 10)
worksheetIntron.write('A1', 'Intron, non-coding and synonymous variants', headingFormat)
worksheetIntron.write_row(1, 0, emptyList, lineFormat)
worksheetIntron.write('A3', 'Sample: '+str(sample))
worksheetIntron.write('A6', 'The following filters for the introns were applied: ')
worksheetIntron.write('B7', 'Coverage >= '+str(min_cov)+'x')
worksheetIntron.write('B8', 'PopAF <= 2 %')
worksheetIntron.write('B9', 'Allele Frequency >= 20%')
worksheetIntron.write('A10', 'For synonymous variants all matching are reported.')

worksheetIntron.write('A12', 'Intron Regions: ')
row = 12
col = 0
for gene in intronDict:
    worksheetIntron.write_row('B'+str(row), [gene]+intronDict[gene])
    row += 1

worksheetIntron.write('A'+str(row+1), 'GATA2 (NM_032638.4) and TP53 synonymous variants: ')
row += 2

for synoVariant in synoVariants:
    worksheetIntron.write_row('B'+str(row), synoVariant)
    row += 1

row += 1
worksheetIntron.write('A'+str(row), 'Coverage below '+str(med_cov)+'x', italicFormat)
row += 2
tableheading = ['RunID', 'DNAnr', 'Gene', 'Chr', 'Pos', 'Ref', 'Alt', 'AF', 'DP',
                'Transcript', 'Mutation cds', 'ENSP', 'Consequence', 'Max popAF', 'Max Pop', 'Callers']
worksheetIntron.write('A'+str(row), 'Intron variants', tableHeadFormat)
worksheetIntron.write_row('A'+str(row+1), tableheading, tableHeadFormat)  # 1 index
row += 1

for line in intron_variants:
    if int(line[8]) < med_cov:
        worksheetIntron.write_row(row, col, line, italicFormat)
    else:
        worksheetIntron.write_row(row, col, line)
    row += 1

row += 1
worksheetIntron.write(row, col, 'Synonymous variants', tableHeadFormat)
worksheetIntron.write_row(row+1, col, tableheading, tableHeadFormat)
row += 2
for line in synoFound:
    if int(line[8]) < med_cov:
        worksheetIntron.write_row(row, col, line, italicFormat)
    else:
        worksheetIntron.write_row(row, col, line)
        row += 1

# CNV GATK
worksheetCNV.conditional_format('G43:G70', {'type': 'cell', 'criteria': 'between',
                                            'minimum': -0.25, 'maximum':  0.2, 'format':   redFormat})
worksheetCNV.conditional_format('I43:I70', {'type': 'cell', 'criteria': 'between',
                                            'minimum': -0.25, 'maximum':  0.2, 'format':   redFormat})
worksheetCNV.set_column('D:D', 23)
worksheetCNV.set_column('E:E', 15)

worksheetCNV.write('A1', 'CNVs found', headingFormat)
worksheetCNV.write_row(1, 0, emptyList, lineFormat)
worksheetCNV.write('A3', 'Sample: '+str(sample))
worksheetCNV.write('A5',
                   'Log2 ratio between -0.25<=x<=0.2 are marked red since they are very weak signals, '
                   + 'and should be interpret with care. ', redFormat)
# Insert png picture to sheets
worksheetCNV.insert_image('A7', snakemake.input.gatk_png)

header = ['Sample', 'Genes', 'Chr', 'Region', 'CytoCoordinates', 'Purity', 'Adapted log2CopyRatio',
          'Adapted CopyNumber', 'log2CopyRatio', 'CopyNumber']
worksheetCNV.write_row('A42', header, tableHeadFormat)
col = 0
row = 42
for line in cnv_lines:
    cn_formula = '= 2 + (J'+str(row+1)+'-2)*(1/F'+str(row+1)+')'  # To get adapted CN based on purity

    worksheetCNV.write_row(row, col, line[0:5])
    worksheetCNV.write_number(row, 5, float(line[5]))  # purity
    worksheetCNV.write_formula(row, 6, '= LOG(J'+str(row+1)+'/2, 2)')  # Adapted log2CR
    worksheetCNV.write_formula(row, 7, cn_formula)  # Adapted CN
    worksheetCNV.write_number(row, 8, float(line[8]))  # log2CR
    worksheetCNV.write_number(row, 9, float(line[9]))  # CN
    row += 1

# CNVkit
worksheetCNVkit.set_column('C:D', 10)
worksheetCNVkit.set_column('B:B', 12)
worksheetCNVkit.set_column('E:E', 15)

worksheetCNVkit.write('A1', 'CNVkit calls', headingFormat)
worksheetCNVkit.write('A3', 'Sample: '+str(sample))
worksheetCNVkit.write('A5', 'Only non-diploid calls or calls with allelic imbalance included')
worksheetCNVkit.write('A7', 'Variant in artefact list ', orangeFormat)

worksheetCNVkit.insert_image('A9', snakemake.input.cnvkit_scatter)

worksheetCNVkit.write_row('A31', relevant_cnvs_header, tableHeadFormat)
row = 31
col = 0
for chromosome in chromosomes:
    for line in relevant_cnvs[chromosome]:
        if len(extractMatchingLines('"' + str(line[1]) + ' ' + str(line[2]) + ' ' +
                                    str(line[3]) + ' ' + str(line[9]) + ' ' + str(line[10]) +
                                    ' ' + str(line[11]) + '"',
                                    snakemake.input.cnvkit_artefact, '-wE')) > 0:
            worksheetCNVkit.write_row(row, col, line, orangeFormat)
            row += 1
        else:
            worksheetCNVkit.write_row(row, col, line)
            row += 1


relevant_chroms = [key for key, value in relevant_cnvs.items() if value != []]
row = row+2
worksheetCNVkit.write(row, col, 'Results per chromosome with aberrant calls',
                      workbook.add_format({'bold': True, 'font_size': 14}))
row = row+1

for i in relevant_chroms:
    if i == 'chrX':
        chr_int = 22
    elif i == 'chrY':
        chr_int = 23
    else:
        chr_int = int(i.replace('chr', ''))-1
    worksheetCNVkit.write(row, col, str(i),  workbook.add_format({'bold': True, 'font_size': 14}))
    row += 1
    worksheetCNVkit.insert_image(row, col, snakemake.input.cnvkit_scatter_perchr[chr_int])
    row += 22
    worksheetCNVkit.write_row(row, col, relevant_cnvs_header, tableHeadFormat)
    row += 1
    for line in relevant_cnvs[i]:
        if len(extractMatchingLines('"' + str(line[1]) + ' ' + str(line[2]) + ' ' +
                                    str(line[3]) + ' ' + str(line[9]) + ' ' + str(line[10]) +
                                    ' ' + str(line[11]) + '"',
                                    snakemake.input.cnvkit_artefact, '-wE')) > 0:
            worksheetCNVkit.write_row(row, col, line, orangeFormat)
            row += 1
        else:
            worksheetCNVkit.write_row(row, col, line)
            row += 1
    row = row+2

# Low Coverage
worksheetLowCov.set_column(1, 3, 10)
worksheetLowCov.set_column(1, 4, 10)
worksheetLowCov.write('A1', 'Mosdepth coverage analysis', headingFormat)
worksheetLowCov.write_row('A2', emptyList, lineFormat)
worksheetLowCov.write('A3', 'Sample: '+str(sample))
description = 'Gene Regions with coverage lower than '+str(min_cov)+'x.'
worksheetLowCov.write('A4', description)
covHeadings = ['Region Name', 'Chr', 'Start', 'Stop', 'Mean Coverage']
worksheetLowCov.write_row('A6', covHeadings, tableHeadFormat)  # 1 index
row = 6  # 0 index
for line in low_cov_lines:
    worksheetLowCov.write_row(row, col, line)
    row += 1

# Hotspot
worksheetHotspot.set_column(1, 2, 10)
worksheetHotspot.write('A1', 'Hotspot Coverage', headingFormat)
worksheetHotspot.write('A3', 'Sample: '+str(sample))
worksheetHotspot.write_row('A5', ['Gene', 'Chr', 'Pos', 'Depth'], tableHeadFormat)
row = 5
for hotLine in hotspotTable:
    if int(hotLine[2]) <= med_cov:
        worksheetHotspot.write_row(row, 0, hotLine, redFormat)
    else:
        worksheetHotspot.write_row(row, 0, hotLine)
    row += 1

# Coverage
worksheetCov.write('A1', 'Average Coverage', headingFormat)
worksheetCov.write_row('A2', emptyList, lineFormat)
worksheetCov.write('A3', 'Sample: '+str(sample))
worksheetCov.write('A4', 'Averge coverage of each region in bedfile from Mosdepth')

tableArea = 'A6:F'+str(len(cov_table_lines)+6)  # rows of full list
headerListDict = [{'header': 'Region Name'}, {'header': 'Chr'}, {'header': 'Start'}, {'header': 'End'},
                  {'header': 'Avg Coverage'}]
worksheetCov.add_table(tableArea, {'data': cov_table_lines, 'columns': headerListDict, 'style': 'Table Style Light 1'})

# QCI
worksheetQCI.set_column('C:C', 10)
worksheetQCI.write('A1', 'Results from QCI ', headingFormat)
worksheetQCI.write_row('A2', emptyList, lineFormat)

worksheetQCI.write('A5', "Analysen utfördes i enlighet med dokumentationen.")
worksheetQCI.write('A6', "Eventuella avikelser:")
iva = ['DNA nr', 'Chromosome', 'Position', 'Gene Region', 'Gene Symbol', 'Transcript ID', 'Transcript Variant',
       'Protein Variant', 'Variant Findings', 'Sample Genotype Quality', 'Read Depth', 'Allele Fraction', 'Translation Impact',
       'dbSNP ID', '1000 Genomes Frequency', 'ExAC Frequency', 'HGMD', 'COSMIC ID', 'Artefacts_without_ASXL1',
       'ASXL1_variant_filter']
worksheetQCI.write_row(9, 0, iva, tableHeadFormat)

# Versions
worksheetVersions.write('A1', 'Version Log', headingFormat)
worksheetVersions.write_row(1, 0, emptyList, lineFormat)
worksheetVersions.write('A3', 'Sample: '+str(sample))
worksheetVersions.write('A5', 'Variant calling reference used: '+str(refV))
worksheetVersions.write('A6', 'Pindel reference used: '+str(refI))
worksheetVersions.write('A7', 'Containers used: ', tableHeadFormat)
containers = [clist for clist in snakemake.params.singularitys.items()]
row = 8
col = 0
for containerTuple in containers:
    container = list(containerTuple)
    worksheetVersions.write_row('A'+str(row), container)
    row += 1


workbook.close()

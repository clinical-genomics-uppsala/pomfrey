#!/bin/python3.6
import sys
import csv
from pysam import VariantFile
import xlsxwriter
from datetime import date
import subprocess
import yaml

# Define sys.argvs
vcf_snv = VariantFile(sys.argv[1])
vcf_indel = VariantFile(sys.argv[2])
cnv_file = sys.argv[3]
cnv_image_path = sys.argv[4]
cartool = sys.argv[5]  # sys.argv[6]
output = sys.argv[6]  # sys.argv[19]
configfile = sys.argv[7]

with open(configfile, 'r') as file:
    config_list = yaml.load(file, Loader=yaml.FullLoader)

runID = config_list['seqID']['sequencerun']  # sys.argv[5]
minCov = int(config_list['cartool']['cov'].split(' ')[0])  # int(sys.argv[7])
medCov = int(config_list['cartool']['cov'].split(' ')[1])  # int(sys.argv[8])
maxCov = int(config_list['cartool']['cov'].split(' ')[2])  # int(sys.argv[9])
bedfile = config_list["bed"]["bedfile"]
pindelBedfile = config_list["bed"]["pindel"]  # sys.argv[10]
# sys.argv[11] #bedfile with - annotated as CNV and exon nummers but exon number not really used right now
cnv_bed_file_path = config_list["CNV"]["bedPoN"]
chrBandFilePath = config_list["CNV"]["cyto"]  # sys.argv[12]
hotspotFile = config_list["bed"]["hotspot"]  # sys.argv[13]
artefactFile = config_list["bed"]["artefact"]  # sys.argv[14]
pindelArtefactFile = config_list["bed"]["pindelArtefact"]  # sys.argv[15]
germlineFile = config_list["bed"]["germline"]  # sys.argv[16]
hematoCountFile = config_list["configCache"]["hemato"]  # sys.argv[17]
variantLog = config_list["configCache"]["variantlist"]  # sys.argv[18]

sample_purity = 0.8

''' Create execl file and sheets. '''
workbook = xlsxwriter.Workbook(output)
worksheetOver = workbook.add_worksheet('Overview')
worksheetShortList = workbook.add_worksheet('ShortList')
worksheetSNV = workbook.add_worksheet('SNVs')
worksheetIndel = workbook.add_worksheet('InDel')  # .... sys.argv[2]
worksheetIntron = workbook.add_worksheet('Intron & Synonymous')
worksheetCNV = workbook.add_worksheet('CNV')
worksheetLowCov = workbook.add_worksheet('Low Coverage')  # ... sys.argv[3]
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

greenFormat = workbook.add_format({'bg_color': '#85e085'})  # , 'font_color': '#b3b3b3'
orangeFormat = workbook.add_format({'bg_color': '#ffd280'})  # font color gray 70%
green_italicFormat = workbook.add_format({'bg_color': '#85e085', 'italic': 'True'})  # 'font_color': '#b3b3b3',
orange_italicFormat = workbook.add_format({'bg_color': '#ffd280', 'italic': 'True'})  # 'font_color': '#b3b3b3',

# Define sample based on annotated vcf
sample = list(vcf_snv.header.samples)[0]
today = date.today()
emptyList = ['', '', '', '', '', '']

shortListGenes = ['ABL1', 'ANKRD26', 'ASXL1', 'ATRX', 'BCOR', 'BCORL1', 'BRAF', 'CALR', 'CBL', 'CBLB', 'CBLC', 'CDKN2A', 'CEBPA',
                  'CSF3R', 'CUX1', 'DDX41', 'DNMT3A', 'ETV6', 'ETNK1', 'TEL', 'EZH2', 'FBXW7', 'FLT3', 'GATA1', 'GATA2', 'GNAS',
                  'HRAS', 'IDH1', 'IDH2', 'IKZF1', 'JAK2', 'JAK3', 'KDM6A', 'KIT', 'KRAS', 'KMT2A', 'MPL', 'MYD88', 'NF1',
                  'NOTCH1', 'NPM1', 'NRAS', 'PDGFRA', 'PHF6', 'PPM1D', 'PTEN', 'PTPN11', 'RAD21', 'RUNX1', 'SAMD9', 'SAMD9L',
                  'SETBP1', 'SF3B1', 'SMC1A', 'SMC3', 'SRP72', 'SRSF2', 'STAG2', 'TET2', 'TP53', 'U2AF1', 'WT1', 'ZRSR2']

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

# VEP fields in list to get index
for x in vcf_snv.header.records:
    if 'CSQ' in str(x) :
        csqIndex = str(x).split('Format: ')[1].strip().strip('">').split('|')


''' Prog Version sheet (9) '''
worksheetVersions.write('A1', 'Version Log', headingFormat)
worksheetVersions.write_row(1, 0, emptyList, lineFormat)
worksheetVersions.write('A3', 'Sample: '+str(sample))
# worksheetVersions.write('A5', 'Variant calling reference used: '+str(refV))
# worksheetVersions.write('A6', 'Pindel reference used: '+str(refI))
worksheetVersions.write('A7', 'Containers used: ', tableHeadFormat)
containers = [clist for clist in config_list['singularitys'].items()]
row = 8
col = 0
for containerTuple in containers:
    container = list(containerTuple)
    worksheetVersions.write_row('A'+str(row), container)
    row += 1


''' QCI sheet (8) '''
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


''' Coverage (7) '''
# Lägga till hela MeanFullCoverage... som filtrerbar lista
# Number of lines in MeanFullCoverage
covFullFile = cartool.replace("_MeanCoverageShortList.csv", "_MeanCoverageFullList.csv")


def file_lengthy(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


numRows = file_lengthy(covFullFile)
worksheetCov.write('A1', 'Average Coverage', headingFormat)
worksheetCov.write_row('A2', emptyList, lineFormat)
worksheetCov.write('A3', 'Sample: '+str(sample))
worksheetCov.write('A4', 'Averge coverage of each region in bedfile')

# Fix data to normal table
tableLines = []
with open(covFullFile) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    # Skip header
    for line in readCSV:
        while not line[-1]:  # Remove empty fields at the end of line.
            line.pop()
        for i in range(1, int((len(line)-1)/5+1)):  # For each region/ with lower cov, create a new row.
            start = 1+5*(i-1)
            end = 1+5*i
            covRow = [line[0]]+line[start:end]
            tableLines.append(covRow)
            # worksheetLowCov.write_row(row,col,[line[0]]+line[start:end])

tableArea = 'A6:F'+str(len(tableLines)+6)  # rows of full list
headerListDict = [{'header': 'Region Name'}, {'header': 'Chr'}, {'header': 'Start'},
                  {'header': 'End'}, {'header': 'Mean Coverage'}, {'header': 'Length of Region'}, ]
worksheetCov.add_table(tableArea, {'data': tableLines, 'columns': headerListDict, 'style': 'Table Style Light 1'})
# worksheetCov.write_row('A6', )#, tableHeadFormat) #table header How sokbar??


''' Hotspot sheet (6) '''
worksheetHotspot.write('A1', 'Hotspot Coverage', headingFormat)  # headerFormat)
worksheetHotspot.write('A3', 'Sample: '+str(sample))
worksheetHotspot.set_column(1, 2, 10)
worksheetHotspot.write_row('A5', ['Chr', 'Pos', 'Depth', 'Gene'], tableHeadFormat)

lowPos = 0
row = 5
hotspotTable = []
# Always the same as bedfile just without region
with open(cartool.replace("_MeanCoverageShortList.csv", "_coverageShortHotspot.tsv"), 'r') as hotFile:
    for dpLine in hotFile:
        cov = dpLine.split("\t")
        chrCov = cov[0]
        posCov = cov[1]
        dp = cov[2].rstrip()

        with open(hotspotFile, "r") as hotspotBed:
            for hotspotLine in hotspotBed:
                hotspot = hotspotLine.split("\t")
                chrHS = hotspot[0]
                lowHS = hotspot[1]
                # highHS = hotspot[2]
                # Always the same as bedfile just without region, How does CARTool handle if hotspot is longer than 1bp?
                if chrCov == chrHS and lowHS == posCov:
                    hotspotTable.append([chrCov, posCov, dp, hotspot[3].rstrip()])

hotspotTable.sort(key=lambda x: float(x[2]))
for hotLine in hotspotTable:
    if int(hotLine[2]) <= medCov:  # How to get the number from configfile?!
        worksheetHotspot.write_row(row, 0, hotLine, redFormat)
        row += 1
        lowPos += 1
    else:
        worksheetHotspot.write_row(row, 0, hotLine)
        row += 1


''' Low Coverage sheet (5) '''
worksheetLowCov.set_column(1, 3, 10)
worksheetLowCov.set_column(1, 4, 10)
# Heading in sheet
worksheetLowCov.write('A1', 'CARTools coverage analysis', headingFormat)
worksheetLowCov.write_row('A2', emptyList, lineFormat)
worksheetLowCov.write('A3', 'Sample: '+str(sample))
description = 'Gene Regions with coverage lower than '+str(minCov)+'x.'
worksheetLowCov.write('A4', description)
covHeadings = ['Region Name', 'Chr', 'Start', 'Stop', 'Mean Coverage', 'Length of Region']
worksheetLowCov.write_row('A6', covHeadings, tableHeadFormat)  # 1 index
row = 6  # 0 index

covLines = []
with open(cartool) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    # Skip header

    for line in readCSV:
        while not line[-1]:  # Remove empty fields at the end of line.
            line.pop()
        for i in range(1, int((len(line)-1)/5+1)):  # For each region/ with lower cov, create a new row.
            start = 1+5*(i-1)
            end = 1+5*i
            covRow = [line[0]]+line[start:end]
            covLines.append(covRow)
            # worksheetLowCov.write_row(row,col,[line[0]]+line[start:end])
            # row += 1
    # sort based on Coverage
    covLines.sort(key=lambda x: x[4])
    for line in covLines:
        worksheetLowCov.write_row(row, col, line)
        row += 1
        # covLines
# Number of low cov regions for the Overview sheet.
lowRegions = row - 6


''' CNV (5) '''
# Import genomic pos to cytoCoord translation to list
chrBands = []
with open(chrBandFilePath, 'r') as chrBandFile:
    for line in chrBandFile:
        chrBands.append(line.split("\t"))

# Load in bedfile with generegions condensed
gene_regions = {}
cnv_bed_file = open(cnv_bed_file_path, 'r')
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
cnv_bed_file.close()

# Process GATK4 CNV File
outLines = []
with open(cnv_file, 'r') as GATK_file:
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
        # Only look at cnv that gatk marked as not neutral   # Borde vi lista alla istallet?
        if lline[5] != '0':
            for gene, coordinates in gene_regions.items():
                if coordinates[0] == chrom:
                    if (start_pos >= int(coordinates[1]) and start_pos <= int(coordinates[2])) or (end_pos >= int(coordinates[1])
                       and end_pos <= int(coordinates[2])) or (start_pos < int(coordinates[1]) and end_pos > int(coordinates[2])):
                        genes.append(gene)
            geneString = ', '.join([x for x in genes if not x.startswith('CNV')])
            logRatio = float(lline[4])
            copyNumberTumor = round(2*pow(2, logRatio), 2)

            outLines.append([sample, geneString, chrom, str(start_pos)+'-'+str(end_pos), cytoCoordString, str(sample_purity),
                             '', '', str(round(logRatio, 4)), str(copyNumberTumor)])

# Output results in xlsx worksheet together with the plot
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
worksheetCNV.insert_image('A7', cnv_image_path)

header = ['Sample', 'Genes', 'Chr', 'Region', 'CytoCoordinates', 'Purity', 'Adapted log2CopyRatio',
          'Adapted CopyNumber', 'log2CopyRatio', 'CopyNumber']  # ,'SamplePurity']
worksheetCNV.write_row('A42', header, tableHeadFormat)
col = 0
row = 42
for line in outLines:
    formula = '= 2 + (J'+str(row+1)+'-2)*(1/F'+str(row+1)+')'  # To get adapted CN

    worksheetCNV.write_row(row, col, line[0:5])  # Must divide to get numbers correct.
    worksheetCNV.write_number(row, 5, float(line[5]))  # purity
    worksheetCNV.write_formula(row, 6, '= LOG(J'+str(row+1)+'/2, 2)')  # Adapted log2CR
    worksheetCNV.write_formula(row, 7, formula)  # Adapted CN
    worksheetCNV.write_number(row, 8, float(line[8]))  # log2CR
    worksheetCNV.write_number(row, 9, float(line[9]))  # CN
    row += 1


''' Intron and synonymous (4) '''
worksheetIntron.set_column('D:E', 10)
worksheetIntron.write('A1', 'Intron, non-coding and synonymous variants', headingFormat)
worksheetIntron.write_row(1, 0, emptyList, lineFormat)
worksheetIntron.write('A3', 'Sample: '+str(sample))

worksheetIntron.write('A6', 'The following filters were applied: ')
worksheetIntron.write('B7', 'Coverage >= 100x')
worksheetIntron.write('B8', 'PopAF <= 2 %')
worksheetIntron.write('B9', 'Allele Frequency >= 20%')
# worksheetIntron.write('B9','Biotype is protein coding') #Nej
# worksheetIntron.write('B10','Consequence not deemed relevant') #nej

worksheetIntron.write('A11', 'Intron Regions: ')
row = 11
col = 0
for gene in intronDict:
    worksheetIntron.write_row('B'+str(row), [gene]+intronDict[gene])
    row += 1

# Add GATA2 synonymous Variants
worksheetIntron.write('A'+str(row+1), 'GATA2 (NM_032638.4) synonymous variants: ')
row += 2
gata2Variants = [[ 'c.1416G>A', 'chr3', '128199889', 'C', 'T'], ['c.1023C>T', 'chr3', '128200782', 'G', 'A'],
                ['c.981G>A', 'chr3', '128202739', 'C', 'T'], ['c.649C>T', 'chr3', '128204792', 'G', 'A'],
                ['c.351C>G', 'chr3', '128205090', 'G', 'C']]
for gata2Variant in gata2Variants:
    worksheetIntron.write_row('B'+str(row), gata2Variant)
    row += 1

row += 1
worksheetIntron.write('A'+str(row), 'Coverage below '+str(medCov)+'x', italicFormat)
row += 2
tableheading = ['RunID', 'DNAnr', 'Gene', 'Chr', 'Pos', 'Ref', 'Alt', 'AF', 'DP',
                'Transcript', 'Mutation cds', 'ENSP', 'Consequence', 'Max popAF', 'Max Pop', 'Callers']
worksheetIntron.write('A'+str(row), 'Intron variants', tableHeadFormat)
worksheetIntron.write_row('A'+str(row+1), tableheading, tableHeadFormat)  # 1 index
row += 1

for snv in vcf_snv.fetch():
    if "PopAF" not in snv.filter.keys():
        if snv.contig in introns:
            for pair in introns[snv.contig]:
                if snv.pos >= pair[0] and snv.pos <= pair[1] and snv.info["AF"][0] >= 0.2:
                    csq = snv.info["CSQ"][0]
                    gene = csq.split("|")[3]
                    transcript = csq.split("|")[10].split(":")[0]
                    if len(csq.split("|")[10].split(":")) > 1:
                        codingName = csq.split("|")[10].split(":")[1]
                    else:
                        codingName = ''
                    ensp = csq.split("|")[11]
                    consequence = csq.split("|")[1]
                    popFreqsPop = ['AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomAD_AF', 'gnomAD_AFR_AF',
                                   'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF',
                                   'gnomAD_OTH_AF', 'gnomAD_SAS_AF']
                    popFreqAllRaw = snv.info["CSQ"][0].split("|")[41:56]  # [42:57]
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
                    try:
                        if snv.info["CALLERS"]:
                            callers = ' & '.join(snv.info["CALLERS"])
                    except KeyError:
                        callers = 'Pisces-multi'

                    line = [runID, sample, gene, snv.contig, str(snv.pos), snv.ref, snv.alts[0], str(snv.info["AF"][0]), str(
                        snv.info["DP"]), transcript, codingName, ensp, consequence, maxPopAf, maxPop, callers]
                    if snv.info["DP"] < medCov:
                        worksheetIntron.write_row(row, col, line, italicFormat)
                    else:
                        worksheetIntron.write_row(row, col, line)
                    row += 1

row += 1
worksheetIntron.write(row, col, 'GATA2 variants', tableHeadFormat)
worksheetIntron.write_row(row+1, col, tableheading, tableHeadFormat)
row += 2
for snv in vcf_snv.fetch('chr3'):
    for gata2Variant in gata2Variants:
        if snv.pos == int(gata2Variant[2]) and snv.alts[0] == gata2Variant[4]:
            csq = snv.info["CSQ"][0]
            gene = csq.split("|")[3]
            transcript = csq.split("|")[10].split(":")[0]
            consequence = csq.split("|")[1]
            if len(csq.split("|")[10].split(":")) > 1:
                codingName = csq.split("|")[10].split(":")[1]
            else:
                codingName = ''
            popFreqsPop = ['AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomAD_AF', 'gnomAD_AFR_AF',
                           'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF',
                           'gnomAD_OTH_AF', 'gnomAD_SAS_AF']
            popFreqAllRaw = snv.info["CSQ"][0].split("|")[41:56]  # [42:57]
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
            try:
                if snv.info["CALLERS"]:
                    callers = ' & '.join(snv.info["CALLERS"])
            except KeyError:
                callers = 'Pisces-multi'

            line = [runID, sample, gene, snv.contig, str(snv.pos), snv.ref, snv.alts[0], str(snv.info["AF"][0]), str(
                snv.info["DP"]), transcript, codingName, ensp, consequence, maxPopAf, maxPop, callers]
            if snv.info["DP"] < medCov:
                worksheetIntron.write_row(row, col, line, italicFormat)
            else:
                worksheetIntron.write_row(row, col, line)
            row += 1


''' Indel sheet (3) '''
# Add genes as info before the actual table. Just use bed table as input? Sort uniq
worksheetIndel.set_column('E:F', 10)  # pos
# worksheetIndel.set_column(1,3,10)
worksheetIndel.write('A1', 'Pindel results', headingFormat)
worksheetIndel.write_row(1, 0, emptyList, lineFormat)
with open(pindelBedfile) as bed:
    genesDup = [line.split("\t")[3].strip() for line in bed]
    genes = set(genesDup)

genesString = ['Genes looked at: ']+list(genes)

for x in vcf_indel.header.records:
    if (x.key == 'reference'):
        refI = x.value
worksheetIndel.write('A3', 'Sample: '+str(sample))
worksheetIndel.write('A4', 'Reference used: '+str(refI))
worksheetIndel.write('A6', 'Genes included: ')
row = 6
for gene in genes:
    worksheetIndel.write('B'+str(row), gene)
    row += 1
worksheetIndel.write(row, 0, 'Coverage below '+str(medCov)+'x', italicFormat)
worksheetIndel.write(row+1, 0, 'Variant in artefact list ', orangeFormat)
worksheetIndel.write(row+2, 0, 'Variants with frequency 0.01 <= AF < 0.05 are located below artefact and germline variants.')
row += 5
tableheading = ['RunID', 'DNAnr', 'Gene', 'Chr', 'Start', 'End', 'SV length', 'Af', 'Ref',
                'Alt', 'Dp', 'Transcript', 'Mutation cds', 'ENSP', 'Max popAF', 'Max Pop', 'IGV']
worksheetIndel.write_row('A'+str(row), tableheading, tableHeadFormat)  # 1 index

orangeIndel = []
whiteIndel = []
whiteIGVIndel = []
underFiveIndel = []
underFiveIGVIndel = []
col = 0

for indel in vcf_indel.fetch():
    # Borde man ta med alla och istallet lagga till en filterkolumn? Hur blir det med icke proteincoding och
    # konsekvens som kanske inte blir samma sak.
    if indel.filter.keys() == ["PASS"]:
        svlen = indel.info["SVLEN"]
        ads = indel.samples[sample]["AD"]
        af = indel.info["AF"]  # int(ads[1])/(int(ads[0]) + int(ads[1])) # Finns val redan
        # ad = ', '.join([str(i) for i in ads])

        if len(indel.alts) == 1:
            alt = indel.alts[0]
        else:
            print(indel.alts)
            sys.exit()

        csqIndel = indel.info["CSQ"][0]  # VEP annotation
        indelGene = csqIndel.split("|")[3]
        maxPopAfIndel = csqIndel.split("|")[56]  # [57]

        if len(maxPopAfIndel) > 1:
            maxPopAfIndel = round(float(maxPopAfIndel), 4)
        maxPopIndel = csqIndel.split("|")[57]  # [61]

        indelTranscript = csqIndel.split("|")[10].split(":")[0]
        if len(csqIndel.split("|")[10].split(":")) > 1:
            indelCodingName = csqIndel.split("|")[10].split(":")[1]
        else:
            indelCodingName = ''
        indelEnsp = csqIndel.split("|")[11]

        indelRow = [runID, sample, indelGene, indel.contig, indel.pos, indel.stop, svlen, af, indel.ref,
                    alt, indel.info["DP"], indelTranscript, indelCodingName, indelEnsp, maxPopAfIndel, maxPopIndel]

        # Mark artefact based on artefactfile
        cmdArt = 'grep -w '+str(indel.pos)+' '+pindelArtefactFile + '| grep -w '+str(indel.contig)
        artLines = subprocess.run(cmdArt, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()

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

# Write to xlsx file
i = 0
for line in whiteIndel:
    if line[10] < medCov:
        worksheetIndel.write_row(row, col, line, italicFormat)
        worksheetIndel.write_url('Q'+str(row+1), whiteIGVIndel[i], string="IGV image")
    else:
        worksheetIndel.write_row(row, col, line)
        worksheetIndel.write_url('Q'+str(row+1), whiteIGVIndel[i], string="IGV image")
    row += 1
    i += 1

for line in orangeIndel:
    if line[10] < medCov:
        worksheetIndel.write_row(row, col, line, orange_italicFormat)
    else:
        worksheetIndel.write_row(row, col, line, orangeFormat)
    row += 1

i = 0
for line in underFiveIndel:
    if line[10] < medCov:
        worksheetIndel.write_row(row, col, line, italicFormat)
        worksheetIndel.write_url('Q'+str(row+1), underFiveIGVIndel[i], string="IGV image")
    else:
        worksheetIndel.write_row(row, col, line)
        worksheetIndel.write_url('Q'+str(row+1), underFiveIGVIndel[i], string="IGV image")
    row += 1
    i += 1


''' SNV sheet (2) '''
worksheetSNV.set_column('E:E', 10)  # Width for position column
# worksheetSNV.set_column(1,11,9) #Width for MaxPopAf column.

for x in vcf_snv.header.records:
    if (x.key == 'reference'):
        refV = x.value
    if (x.key == 'VEP'):
        vepline = x.value

# Headers before variants
worksheetSNV.write('A1', 'Variants found', headingFormat)
worksheetSNV.write('A3', 'Sample: '+str(sample))
worksheetSNV.write('A4', 'Reference used: '+str(refV))
worksheetSNV.write('A6', 'VEP: '+vepline)  # , textwrapFormat)
worksheetSNV.write('A8', 'The following filters were applied: ')
worksheetSNV.write('B9', 'Coverage >= 100x')
worksheetSNV.write('B10', 'Population freq (KGP, gnomAD, NHLBI_ESP ) <= 2%')
worksheetSNV.write('B11', 'Biotype is protein coding')
worksheetSNV.write('B12', 'Consequence not deemed relevant')

worksheetSNV.write('A14', 'Coverage below '+str(medCov)+'x', italicFormat)
worksheetSNV.write('A15', 'Variant in artefact list ', orangeFormat)
worksheetSNV.write('A16', 'Variant likely germline', greenFormat)
worksheetSNV.write('A17', 'Variants with frequency 0.03 <= AF < 0.05 are located below artefact and germline variants.')

# Variant table
tableheading = ['RunID', 'DNAnr', 'Gene', 'Chr', 'Pos', 'Ref', 'Alt', 'AF', 'DP', 'Transcript', 'Mutation cds', 'ENSP',
                'Consequence', 'COSMIC ids on position', 'N COSMIC Hemato hits on position', 'Clinical significance', 'dbSNP',
                'Max popAF', 'Max Pop', 'IGV image', 'Callers']
worksheetSNV.write_row('A19', tableheading, tableHeadFormat)  # 1 index
row = 19  # 0 index
col = 0
white = []
green = []
orange = []
whiteIGV = []
underFive = []  # put after green and orange but still white
underFiveIGV = []  # put after green and orange but still white
shortListSNV = []  # Reported genes only
shortListSNVigv = []  # Reported genes only
greenShortList = []  # Germline in reported genes

for record in vcf_snv.fetch():
    # Check if the synonomus variant exists in the COSMIC hemato file or is a splice variant
    synoCosmicN = 0
    spliceVariant = False
    csq = record.info["CSQ"][0]
    consequence = csq.split("|")[1]

    if record.filter.keys() == ["Syno"]:  # Only if Syno not and popAF.   any(x in "Syno" for x in record.filter.keys()):
        synoCosmicVepList = [cosmic for cosmic in csq.split("|")[17].split(
            "&") if cosmic.startswith('CO')]  # Get all cosmicID in list
        # COSMIC Hemato
        if len(synoCosmicVepList) != 0:
            for synoCosmicId in synoCosmicVepList:
                cmdCosmic = 'grep -w '+synoCosmicId+' '+hematoCountFile+' | cut -f 16 '
                synoCosmicNew = subprocess.run(cmdCosmic, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
                if len(synoCosmicNew) == 0:
                    synoCosmicNew = 0
                synoCosmicN += int(synoCosmicNew)
        # Splice variant
        if 'splice' in consequence:
            spliceVariant = True

    if record.filter.keys() == ["PASS"] or synoCosmicN != 0 or spliceVariant:

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

        if af >= 0.03:
            try:
                if record.info["CALLERS"]:
                    callers = ' & '.join(record.info["CALLERS"])
            except KeyError:
                callers = 'Pisces-multi'
#            csq = record.info["CSQ"][0]
            gene = csq.split("|")[3]
            clinical = csq.split("|")[58]  # [59]
            existing = csq.split("|")[17].split("&")

            # rs IDs use more than just the first!
            rsList = [rs for rs in existing if rs.startswith('rs')]
            if len(rsList) == 0:
                rs = ''
            else:
                rs = ', '.join(rsList)
                # rs = rsList[0]

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
                    cmdCosmic = 'grep -w '+cosmicId+' '+hematoCountFile+' | cut -f 16 '
                    cosmicNew = subprocess.run(cmdCosmic, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
                    if len(cosmicNew) == 0:
                        cosmicNew = 0
                    cosmicN += int(cosmicNew)

            transcript = csq.split("|")[10].split(":")[0]
            if len(csq.split("|")[10].split(":")) > 1:
                codingName = csq.split("|")[10].split(":")[1]
            else:
                codingName = ''
            # consequence = csq.split("|")[1]
            ensp = csq.split("|")[11]

            # Population allel freq
            # maxPopAf = record.info["CSQ"][0].split("|")[57] #[60]
            # if len(maxPopAf) > 1:
            #     maxPopAf = round(float(maxPopAf),4)
            # maxPop = record.info["CSQ"][0].split("|")[58] #[61]

            popFreqsPop = ['AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF',
                           'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF']
            popFreqAllRaw = record.info["CSQ"][0].split("|")[41:56]  # [42:57]
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
            # IGV image path for each SNV
            igv = "external:IGV/"+gene+"-"+record.contig+"_"+str(int(record.pos)-1)+"_"+str(int(record.pos)-1+len(alt))+".svg"

            snv = [runID, sample, gene, record.contig, record.pos, record.ref, alt, af, record.info["DP"], transcript,
                   codingName, ensp, consequence, cosmicVep, cosmicN, clinical, rs, maxPopAf, maxPop, '', callers]

            # Append line with sample and rundate to rolling list of artefacts..
            with open(variantLog, "a") as appendfile:
                variants = snv+["\n"]
                appendfile.write('\t'.join(str(e) for e in variants))

            # Artefact_file
            cmdArt = 'grep -w '+str(record.pos)+' '+artefactFile
            artLines = subprocess.run(cmdArt, stdout=subprocess.PIPE, shell='TRUE').stdout.decode(
                'utf-8').strip()  # What happens if two hits?!!!
            artefact_variant = 0

            for artLine in artLines.split("\n"):
                # if pos exists and match in artefact file.
                if artLine and record.ref == artLine.split()[2] and alt == artLine.split()[3]:
                    orange.append(snv)
                    artefact_variant = 1
                    break
            if artefact_variant == 0:
                # Germline /gluster-storage-volume/projects/wp2/nobackup/Twist_Myeloid/Artefact_files/Low_VAF_SNVs.txt
                cmdGerm = 'grep -w '+str(record.pos)+' '+germlineFile
                germLines = subprocess.run(cmdGerm, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
                germline_variant = 0
                for germLine in germLines.split("\n"):
                    # if exists in germline file
                    if germLine and record.ref == germLine.split()[2] and alt == germLine.split()[3]:
                        green.append(snv)
                        germline_variant = 1
                        if gene in shortListGenes:
                            greenShortList.append(snv)
                        break
                if germline_variant == 0:
                    if float(af) < 0.05:
                        underFive.append(snv)
                        underFiveIGV.append(igv)
                    else:
                        white.append(snv)
                        whiteIGV.append(igv)
                    if gene in shortListGenes:
                        shortListSNV.append(snv)
                        shortListSNVigv.append(igv)
# Write to xlsx file
i = 0
for line in white:
    if line[8] < medCov:
        worksheetSNV.write_row(row, col, line, italicFormat)
        worksheetSNV.write_url('T'+str(row+1), whiteIGV[i], string="IGV image")
    else:
        worksheetSNV.write_row(row, col, line)
        worksheetSNV.write_url('T'+str(row+1), whiteIGV[i], string="IGV image")
    row += 1
    i += 1

for line in green:
    if line[8] < medCov:
        worksheetSNV.write_row(row, col, line, green_italicFormat)
    else:
        worksheetSNV.write_row(row, col, line, greenFormat)
    row += 1

for line in orange:
    if line[8] < medCov:
        worksheetSNV.write_row(row, col, line, orange_italicFormat)
    else:
        worksheetSNV.write_row(row, col, line, orangeFormat)
    row += 1

i = 0
for line in underFive:
    if line[8] < medCov:
        worksheetSNV.write_row(row, col, line, italicFormat)
        worksheetSNV.write_url('T'+str(row+1), underFiveIGV[i], string="IGV image")
    else:
        worksheetSNV.write_row(row, col, line)
        worksheetSNV.write_url('T'+str(row+1), underFiveIGV[i], string="IGV image")
    row += 1
    i += 1


''' Reported variants '''
worksheetShortList.set_column('E:E', 10)
# Variants or snv rows from SNV sheet.
worksheetShortList.write('A1', 'Variants in genes to report', headingFormat)
worksheetShortList.write('A3', 'Sample: '+str(sample))
worksheetShortList.write('A6', 'VEP: '+vepline)  # , textwrapFormat)
worksheetShortList.write('A8', 'The following filters were applied: ')
worksheetShortList.write('B9', 'Coverage >= 100x')
worksheetShortList.write('B10', 'Population freq (KGP, gnomAD, NHLBI_ESP ) <= 2%')
worksheetShortList.write('B11', 'Biotype is protein coding')
worksheetShortList.write('B12', 'Consequence not deemed relevant')


worksheetShortList.write('A14', 'Short list of variants in genes to report:')
worksheetShortList.write_row(14, 0, shortListGenes)
worksheetShortList.write('A16', 'For all variants see: ')
worksheetShortList.write_url('B16', "internal:'SNVs'!A1", string='SNVs')

worksheetShortList.write('A18', 'Coverage below '+str(medCov)+'x', italicFormat)
worksheetShortList.write('A19', 'Variant likely germline', greenFormat)
worksheetShortList.write_row('A21', tableheading, tableHeadFormat)  # 1 index
row = 21  # 0 index
i = 0
for line in shortListSNV:
    if line[8] < medCov:
        worksheetShortList.write_row(row, col, line, italicFormat)
        worksheetShortList.write_url('T'+str(row+1), shortListSNVigv[i], string="IGV image")
    else:
        worksheetShortList.write_row(row, col, line)
        worksheetShortList.write_url('T'+str(row+1), shortListSNVigv[i], string="IGV image")
    row += 1
    i += 1
for line in greenShortList:
    if line[8] < medCov:
        worksheetShortList.write_row(row, col, line, green_italicFormat)
    else:
        worksheetShortList.write_row(row, col, line, greenFormat)
    row += 1


''' Overview sheet (1) '''
worksheetOver.write(0, 0, sample, headingFormat)
worksheetOver.write(1, 0, "RunID: "+runID)
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
worksheetOver.write_url(11, 0, "internal:'Intron'!A1", string='Intron & synonymous variants')
worksheetOver.write_url(12, 0, "internal:'CNV'!A1", string='CNVs found with GATK4')
worksheetOver.write_url(13, 0, "internal:'Low Coverage'!A1", string='Positions with coverage lower than '+str(minCov)+'x')
worksheetOver.write_url(14, 0, "internal:'Hotspot'!A1", string='Coverage of hotspot positions')
worksheetOver.write_url(15, 0, "internal: 'Coverage'!A1", string='Average coverage of all regions in bed')
worksheetOver.write_url(16, 0, "internal:'Version'!A1", string='Version Log')
worksheetOver.write_row(17, 0, emptyList, lineFormat)


# Add avg. cov and clonalisy
cartoolLog = cartool.replace("_MeanCoverageShortList.csv", "_Log.csv")

cmdAvgCov = 'grep Depth '+cartoolLog+' | cut -d"," -f2 | cut -f1 -d" "'
avgCov = subprocess.run(cmdAvgCov, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()

duplicationFile = cartool.replace("_MeanCoverageShortList.csv", "_DuplicationMetrics.txt")
cmdDupl = 'grep -A1 PERCENT '+duplicationFile+' | tail -1 | cut -f9'
duplicateLevel = subprocess.run(cmdDupl, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()

breadthCmd = 'grep "Mean Coverage Breadth:" '+cartoolLog + ' | cut -f2- -d"," '
breadth = subprocess.run(breadthCmd, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()


worksheetOver.write_row(19, 0, ['RunID', 'DNAnr', 'Avg. coverage [x]', 'Duplicationlevel [%]',
                                str(minCov)+'x', str(medCov)+'x', str(maxCov)+'x'], tableHeadFormat)
worksheetOver.write_row(20, 0, [runID, sample, avgCov, str(round(float(duplicateLevel)*100, 2))]+breadth.split(','))

if lowPos == 0:  # From Hotspot sheet
    worksheetOver.write(23, 0, 'Number of positions from the hotspot list not covered by at least '+str(medCov)+'x: ')
    worksheetOver.write(24, 0, str(lowPos))
else:
    worksheetOver.write(23, 0, 'Number of positions from the hotspot list not covered by at least '+str(medCov)+'x: ')
    worksheetOver.write(24, 0, str(lowPos), redFormat)
    worksheetOver.write_url(25, 0, "internal:'Hotspot'!A1", string='For more detailed list see hotspotsheet ')

worksheetOver.write(26, 0, 'Number of regions not covered by at least '+str(minCov)+'x: ')  # From Cov sheet
worksheetOver.write(27, 0, str(lowRegions))  # From Cov sheet

worksheetOver.write(29, 0, 'Bedfile: '+bedfile)
worksheetOver.write(30, 0, 'Hotspotlist: '+hotspotFile)
worksheetOver.write(31, 0, 'Artefact file: '+artefactFile)
worksheetOver.write(32, 0, 'Germline file: '+germlineFile)
worksheetOver.write(33, 0, 'Pindel artefact file: '+pindelArtefactFile)


''' Prog Version sheet (8), added last '''
worksheetVersions.write('A5', 'Variant calling reference used: '+str(refV))
worksheetVersions.write('A6', 'Pindel reference used: '+str(refI))


workbook.close()

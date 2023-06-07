#!/bin/python3
import sys
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1])  # dosen't matter if bgziped or not. Automatically recognizes
# Consequences to keep (not to filter out)
consequences = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant',
                'stop_lost', 'start_lost', 'transcript_amplification', 'inframe_insertion', 'inframe_deletion',
                'missense_variant', 'protein_altering_variant', 'splice_region_variant', 'incomplete_terminal_codon_variant',
                'coding_sequence_variant', 'mature_miRNA_variant']
# Add new filter descriptions to new header
new_header = vcf_in.header
new_header.filters.add("PopAF", None, None, "Population AF over two percent")
new_header.filters.add("DP100", None, None, "DP is lower than 100x")
new_header.filters.add("ProtCode", None, None, "Biotype is not protein coding")
new_header.filters.add("Conseq", None, None, "Consequence is not deemed relevant (See XX for more info)")
new_header.filters.add("Syno", None, None, "Consequence is synonymous variant")

# start new vcf with the new_header
vcf_out = VariantFile(sys.argv[2], 'w', header=new_header)


for record in vcf_in.fetch():
    # import pdb; pdb.set_trace()
    # Filter based on coverage in DP field in INFO.
    if record.info["DP"] <= 100:
        record.filter.add("DP100")

    # For every transcript?? Seperated with , not all differs, popfreq is the same on all.
    # Filter on known population freq (KGP phase3 and gnomAD r2.1 exomes only )
    vep = record.info["CSQ"][0].split("|")
    # if not vep[57] == '' and float(vep[57]) >= 0.02 : #vep[60]
    #      record.filter.add("PopAF")

    popFreqAll = vep[41:56]  # [42:57]
    popFreqs = [x for x in popFreqAll if x]
    if any(popFreqs) and any(float(x) >= 0.02 for x in popFreqs):
        record.filter.add("PopAF")

    bioType = 0
    conseq = 0

    for csq in record.info["CSQ"]:
        vep = csq.split("|")
        # Filter out not protein_coding
        if vep[7] == 'protein_coding':  # What if protein_coding&...?
            bioType = 1
            if 'synonymous_variant' in vep[1]:
                record.filter.add("Syno")
            if any(x in vep[1] for x in consequences):
                conseq = 1  # Ok consequnce, no filter

    if bioType == 0:  # if not protein coding
        record.filter.add("ProtCode")
    if conseq == 0:  # if not wanted consequences
        record.filter.add("Conseq")


# Filter based on BIOTYPE? exon, intron splice... splice_region_variant intron_variant
# BIOTYPE 7: protein_coding
# Consequence 1: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
    vcf_out.write(record)


# recod.contig =chr
# record.pos =pos
# record.ref
# record.alts

# # recod.info["CSQ"] = vep annotering!
# Allele 0 |Consequence 1 |IMPACT 2 |SYMBOL 3 |Gene 4 |Feature_type 5 |Feature 6 |BIOTYPE 7 |EXON 8 |INTRON 9 |HGVSc 10 |HGVSp 11
# |cDNA_position 12 |CDS_position 13 |Protein_position 14 |Amino_acids 15 |Codons 16 |Existing_variation 17 |DISTANCE 18
# |STRAND 19 |FLAGS 20 |VARIANT_CLASS 21 |SYMBOL_SOURCE 22 |HGNC_ID 23 |CANONICAL 24 |TSL 25 |APPRIS 26 |CCDS 27 |ENSP 28
# |SWISSPROT 29 |TREMBL 30 |UNIPARC 31 |REFSEQ_MATCH 32 |SOURCE 33 |GIVEN_REF 34 |USED_REF 35 |BAM_EDIT 36 |GENE_PHENO 37
# |SIFT 38 |PolyPhen 39 |DOMAINS 40 |HGVS_OFFSET 41 |AF 42 |AFR_AF 43 |AMR_AF 44 |EAS_AF 45 |EUR_AF 46 |SAS_AF 47 |gnomAD_AF 48
# |gnomAD_AFR_AF 49|gnomAD_AMR_AF 50 |gnomAD_ASJ_AF 51 |gnomAD_EAS_AF 52 |gnomAD_FIN_AF 53 |gnomAD_NFE_AF 54 |gnomAD_OTH_AF 55
# |gnomAD_SAS_AF 56 |MAX_AF 57 |MAX_AF_POPS 58 |CLIN_SIG 59 |SOMATIC 60 |PHENO 61 |PUBMED 62 |MOTIF_NAME 63 |MOTIF_POS 64
# |HIGH_INF_POS 65 |MOTIF_SCORE_CHANGE 66
#################################
# Allele 0 |Consequence 1 |IMPACT 2 |SYMBOL 3 |Gene 4 |Feature_type 5 |Feature 6 |BIOTYPE 7 |EXON 8 |INTRON 9 |HGVSc 10 |HGVSp 11
# |cDNA_position 12 |CDS_position 13 |Protein_position 14 |Amino_acids 15 |Codons 16 |Existing_variation 17 |DISTANCE 18
# |STRAND 19 |FLAGS 20 |VARIANT_CLASS 21 |SYMBOL_SOURCE 22 |HGNC_ID 23 |CANONICAL 24 |TSL 25 |APPRIS 26 |CCDS 27 |ENSP 28
# |SWISSPROT 29 |TREMBL 30 |UNIPARC 31 |REFSEQ_MATCH 32 |SOURCE 33 |GIVEN_REF 34 |USED_REF 35 |BAM_EDIT 36 |GENE_PHENO 37
# |SIFT 38 |PolyPhen 39 |DOMAINS 40 |miRNA 41 |HGVS_OFFSET 42 |AF 43 |AFR_AF 44 |AMR_AF 45 |EAS_AF 46 |EUR_AF 47 |SAS_AF 48
# |AA_AF 49 |EA_AF 50 |gnomAD_AF 51 |gnomAD_AFR_AF 52 |gnomAD_AMR_AF 53 |gnomAD_ASJ_AF 54 |gnomAD_EAS_AF 55 |gnomAD_FIN_AF 56
# |gnomAD_NFE_AF 57 |gnomAD_OTH_AF 58 |gnomAD_SAS_AF 59 |MAX_AF 60 |MAX_AF_POPS 61 |CLIN_SIG 62 |SOMATIC 63 |PHENO 64
# |PUBMED 65 |MOTIF_NAME 66 |MOTIF_POS 67 |HIGH_INF_POS 68 |MOTIF_SCORE_CHANGE 69

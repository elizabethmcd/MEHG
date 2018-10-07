#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description = "Manually change annotations in a gff file based on locus tags and output new genbank files")
parser.add_argument('gff_file', metavar='GFF', help='Annotation file from Prokka in GFF format')
parser.add_argument('locus_list', metavar='LOCUS', help="Locus tag list to replace annotations for")
parser.add_argument('--outputGFF', default='gene_annot.gff', help="Converted GFF file")
parser.add_argument('--outputGBK', default='gene_annot.gbk', help="Converted genbank file")

args = parser.parse_args()

# Input and output files
GFF = args.gff_file
LOCUS = args.locus_list
OUT_GEN = open(args.outputGFF, "w")

# Read in locus tag list
with open(LOCUS) as loci:
    loci_list = loci.read().splitlines()

# Read in the GFF file
annot_in = open(GFF, "r")

for line in annot_in: 
    for loci in loci_list:
        if loci in line:
            print(line)


# Parse Genbank file 
# for record in SeqIO.parse(GBK, "genbank"):
#     contig = record.name
#     for f in record.features:
#         if f.type == 'CDS' and 'locus_tag' in f.qualifiers:
#             locus = f.qualifiers["locus_tag"][0]
#             if locus in loci_list:
#                 product = MutableSeq(f.qualifiers["product"][0])
#                 mehg = "mercury methylation"
#                 new = product == mehg
#                 print(new)
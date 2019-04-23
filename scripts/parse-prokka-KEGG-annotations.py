#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import pandas as pd

# Arguments

parser = argparse.ArgumentParser(description = "Parse Prokka annotations from a genbank file and KEGG annotations to combine in a master annotation file")
parser.add_argument('gbk_file', metavar='GBK', help='Annotation file from Prokka in Genbank format')
parser.add_argument('kegg_file', metavar='KEGG', help='KEGG annotation file in format of protein ID in col1 and KEGG KO in col2')
parser.add_argument('--prokka_annotations', default='prokka_annot.txt', help="Output: Prokka annotation fiile (Default: prokka_annot.txt)")
parser.add_argument('--combo_annotations', default='combined-annotations.txt', help="Output: Combined annotation file (Default: combined-annotations.txt")
args = parser.parse_args()

# Input and output files
GBK = args.gbk_file
KEGG = args.kegg_file
OUT_PROKKA = open(args.prokka_annotations, "w")
IN_PROKKA = open(args.prokka_annotations, "r")

# Headers for Anvi'o output
OUT_PROKKA.write("protein_id\tbeg\tend\tprokka_function\n")

# Parse KO file 
ko_file = pd.read_table(KEGG, sep="\t", names=['tag', 'KO', 'annotation'])

# Parse the GBK file
for record in SeqIO.parse(GBK, "genbank"):
    contig = record.name
    for f in record.features:
        if f.type == 'CDS':
            locus_tag, span, function = (f.qualifiers["protein_id"][0], f.location, f.qualifiers["product"][0])
            protein_id = locus_tag.replace("X:", "")
            beg = span.start # biopython already reads in the genbank file to start counting from 0, so don't have to subtract 1
            end = span.end
            # Write out to gene calls and annotation files
            OUT_PROKKA.write('%s\t%d\t%d\t%s\n' % (protein_id, beg, end, function))

prokka_file = pd.read_table(IN_PROKKA, sep="\t")
combo_annotation = prokka_file.merge(ko_file[['tag', 'KO', 'annotation']], how="left", left_on="protein_id", right_on="tag")
combo_annotation.to_csv("combined-annotations.csv", sep=",", index=False)
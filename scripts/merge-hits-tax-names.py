#! /usr/bin/env python

import pandas as pd
import argparse, os

# Arguments
parser = argparse.ArgumentParser(description = "Merge marker hits with full taxonomy in metadata")
parser.add_argument('hits_file', metavar='HITS', help="List of genome hits of a marker")
parser.add_argument('metadata', metavar='META', help="Location of metadata file")
parser.add_argument('--outfile', default='itol.out.txt', help="Output: Outputs ITOL specific tree-names files")

args = parser.parse_args()

# Input and output files 
HITS = args.hits_file
META = args.metadata
OUTF = open(args.outfile, "w")

# ITOL header 
OUTF.write("LABELS\nSEPARATOR TAB\nDATA\n")

# Hits
marker_hits = pd.read_table(HITS, names=['genome_name'])

# Metadata file
metadata = pd.read_table(META, sep=",")

# Merge
marker_hits_taxonomy = marker_hits.merge(metadata[['genome_name', 'Phylum']], how="left", left_on="genome_name", right_on="genome_name")
marker_hits_names = marker_hits_taxonomy[['genome_name', "Phylum"]]
marker_hits_names.to_csv(OUTF, sep="\t", header=False, index=False)
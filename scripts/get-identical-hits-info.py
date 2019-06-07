#! /usr/bin/env python

import pandas as pd
import argparse, os

# Arguments
parser = argparse.ArgumentParser(description = "Merge marker hits with full taxonomy in metadata")
parser.add_argument('hits_list', metavar='HITS', help="List of genomes with identical marker hits")
parser.add_argument('metadata', metavar='META', help="Location of metadata file")
parser.add_argument('--outfile', default='identical-information.txt', help="Information on identical hits output")
args = parser.parse_args()

# Input and output files 
HITS = args.hits_list
META = args.metadata
OUTF = open(args.outfile, "w")

# Hits
marker_hits = pd.read_table(HITS)

# Metadata file
metadata = pd.read_table(META, sep=",")

# Merge
genome1 = marker_hits.merge(metadata[['genome_name', 'taxonomy', 'Phylum']], how="left", left_on="genome1", right_on="genome_name")
genome1info = genome1[['genome1', 'taxonomy', 'Phylum', 'genome2']]
genome1info.columns = ['genome1', 'taxonomy1', 'phylum1', 'genome2']
genome2 = genome1info.merge(metadata[['genome_name', 'taxonomy', 'Phylum']], how="left", left_on="genome2", right_on="genome_name")
genome2info = genome2[['genome1', 'taxonomy1', 'phylum1', 'genome2', 'taxonomy', 'Phylum']]
genome2info.columns = ['genome1', 'taxonomy1', 'phylum1', 'genome2', 'taxonomy2', 'phylum2']
for col in genome2info.columns:
    print(col)
genome2info.to_csv(OUTF, sep="\t", header=True, index=False)
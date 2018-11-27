#! /usr/bin/env python

import pandas as pd 

# Metadata
metadata = pd.read_table("all-arch-bac-metadata.tsv", sep="\t")

# Hits
hits = pd.read_table("all-hgcA-list-genomes.txt", sep="\t", names=['genome_name', 'locus_tag'])

# Combine hits with metadata
hits_taxonomy = hits.merge(metadata[['new_name', 'taxonomy']], how="left", left_on="genome_name", right_on="new_name")
tree_names=hits_taxonomy[['locus_tag','taxonomy']]
tree_names.to_csv("hgcA-tree-names.tsv", sep="\t", index=False)

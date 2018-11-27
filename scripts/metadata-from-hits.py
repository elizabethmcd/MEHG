#! /usr/bin/env python

import pandas as pd 

# Bacterial and archaeal metadata with new names
arch_metadata = pd.read_table("arch-metadata.tsv", sep="\t")
bac_metadata = pd.read_table("bac-metadata.tsv", sep="\t")

# Hits
arch_hits = pd.read_table("arch-hgcA-list-genomes-hits.txt", sep="\t", names=['genome_name', 'locus_tag'])

# Combine hits with metadata
arch_hits_taxonomy = arch_hits.merge(arch_metadata[['new_name', 'taxonomy']], how="left", left_on="genome_name", right_on="new_name")
arch_tree_names=arch_hits_taxonomy[['locus_tag','taxonomy']]
arch_tree_names.to_csv("arch-hgcA-tree-names.tsv", sep="\t", index=False)

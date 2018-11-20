#! /usr/bin/env python

import pandas as pd 

# All genomes taxonomy
all_tax = pd.read_table("all-mags-refs-genomes.txt", sep="\t", names=['genome_name', 'taxonomy'])

# Archaea names
arch_df = pd.read_table("arch-genomes-new-names.txt", sep="\t", names=['genome_name', 'old_name_fna', 'new_name', 'new_name_fna'])
# Bacteria names
bac_df = pd.read_table("bac-genomes-new-names.txt", sep="\t", names=['genome_name', 'old_name_fna', 'new_name', 'new_name_fna'])

# Create metadata files
arch_metadata=arch_df.merge(all_tax[['genome_name','taxonomy']],how="left", left_on="genome_name", right_on="genome_name")
arch_metadata=arch_metadata[['genome_name', 'new_name', 'taxonomy']]
bac_metadata=bac_df.merge(all_tax[['genome_name', 'taxonomy']],how="left", left_on="genome_name",right_on="genome_name")
bac_metadata=bac_metadata[['genome_name', 'new_name', 'taxonomy']]

arch_metadata.to_csv("arch-metadata.tsv", sep="\t", index=False)
bac_metadata.to_csv("bac-metadata.tsv",sep="\t",index=False)
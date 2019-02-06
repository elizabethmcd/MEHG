#! /usr/bin/env python 

import pandas as pd 

gtdb_classf = pd.read_table("all-gtdbk-classifications.txt", sep="\t", names=['genome_name', 'gtdb_taxonomy'])
all_metadata = pd.read_table("all-arch-bac-metadata.tsv", sep="\t")

all_classifications = gtdb_classf.merge(all_metadata[['genome_name','new_name', 'taxonomy']], how="left", left_on="genome_name", right_on="new_name")
all_classifications.to_csv("all-metadata-gtdb-names.tsv", sep="\t", index=False)
#! /Users/emcdaniel/anaconda3/bin/python

import pandas as pd 

# UBA dataset
locus = pd.read_table("UBA-list-genome-tags.txt", sep="\t", header=None, names=["locus_tag"])
met = pd.read_table("UBA-methylating-bins-metadata.tsv")

locus["genome_name"]=locus["locus_tag"].str.split("_").str.get(0)
tax_names = locus.merge(met[["UBA-Bin", "Classification"]], how="left", left_on="genome_name", right_on="UBA-Bin")
tax_names[["locus_tag", "Classification"]].to_csv("UBA-meth-locus-names.tsv", sep="\t", index=False)

# Peat dataset 
peat_locus = pd.read_table("peat-locus-tags.txt", sep="\t", header=None, names=["locus_tag"])
peat_tax = pd.read_table("peat-methylators-taxonomy.txt", sep="\t", header=None, names=["genome_name", "taxonomy"])
peat_locus["genome_name"] = peat_locus["locus_tag"].str.slice(0,15)
peat_tax = peat_locus.merge(peat_tax, how="left", left_on="genome_name", right_on="genome_name")
peat_tax[["locus_tag", "taxonomy"]].to_csv("Peat-meth-locus-names.tsv", sep="\t", index=False)




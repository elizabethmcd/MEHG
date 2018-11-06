#! /usr/bin/env python 

# Read in loci list of hits to get genome names, organize a master file and rename directory/file names based off of consistent naming code 
import pandas as pd 

loci_list = open("hgcA-locus-list.txt", "r")
with open("genomelist.txt", "a") as f:
    for line in loci_list: 
        if line.startswith("K_DeepCast"):
            new="_".join(line.split("_", 5)[:5])
            f.write(new+"\n")
        elif line.startswith("M_DeepCast"):
            new="_".join(line.split("_", 5)[:5])
            f.write(new+"\n")
        elif line.startswith("UBA"):
            new="_".join(line.split("_", 1)[:1])
            f.write(new+"\n")
        elif line.startswith("GCF"):
            new="_".join(line.split("_", 2)[:2])
            f.write(new+"\n")
        elif line.startswith("GCA"):
            new="_".join(line.split("_", 2)[:2])
            f.write(new+"\n")
        elif line.startswith("KMB"):
            new="_".join(line.split("_", 1)[:1])
            f.write(new+"\n")
        elif ".genes" in line:
            new=line.split(".genes")[0]
            f.write(new+"\n")
loci_list.close()

# Metadata files 
tang=pd.read_table("tang-classf.txt", sep="\t", names=['genome_name', 'taxonomy'], index_col=None)
mendota=pd.read_table("mendota-classf.txt", sep="\t", names=['genome_name', 'taxonomy'])
trout=pd.read_table("trout-bog-classf.txt", sep="\t", header=0)
references=pd.read_table("ref-classf.txt", sep="\t", header=0)
UBA=pd.read_table("UBA-methylating-bins-metadata.tsv", sep="\t", header=0)
peat=pd.read_table("Peat-taxonomy.txt", sep="\t", names=['genome_name', 'taxonomy'])
aquifer=pd.read_table("2500-taxonomy.txt", sep="\t", names=['genome_name', 'taxonomy'])

# Genome list 
# generated after getting rid of duplicate lines and non-hits, and subsets to make matches 
genome_list = pd.read_table("meth-genomes-list.txt", names=['genome_name'])
UBA_list = pd.read_table("UBA-genome-list.txt", names=['genome_name'])
genbank_list=pd.read_table("genbank-meths.txt", names=['genome_name'])
mendota_list=pd.read_table("mendota-genome-list.txt", names=['genome_name'])
tang_list=pd.read_table("tang-genome-list.txt", names=['genome_name'])
trout_list=pd.read_table("trout-bog-genome-list.txt", names=['genome_name'])
ref_list=pd.read_table("ref-genome-list.txt", names=['genome_name'])

# UBA matches 
UBA_table = UBA[['UBA-Bin', 'Classification']]
UBA_tax= UBA_list.merge(UBA_table[['UBA-Bin', 'Classification']], how="left", left_on="genome_name", right_on="UBA-Bin")
UBA_final=UBA_tax[['genome_name', 'Classification']]
UBA_final.rename(columns={'Classification':'taxonomy'}, inplace=True)

# peat and aquifer matches 
frames = (peat,aquifer)
genbank_full = pd.concat(frames)
genbank_tax = genbank_list.merge(genbank_full[['genome_name', 'taxonomy']], how="left", left_on="genome_name", right_on="genome_name")

# mendota
mendota_tax = mendota_list.merge(mendota[['genome_name', 'taxonomy']], how="left", left_on="genome_name", right_on="genome_name")

# tanganyika
tang_tax = tang_list.merge(tang[['genome_name', 'taxonomy']], how="left", left_on='genome_name', right_on='genome_name')

# trout bog
trout_tax=trout_list.merge(trout[['genome_name', 'taxonomy']], how="left", left_on='genome_name', right_on='genome_name')

# references
ref_tax = ref_list.merge(references[['genome_name', 'taxonomy']], how='left', left_on='genome_name', right_on='genome_name')

# write out separate dataframes 
UBA_final.to_csv("UBA-tax.txt", sep="\t", index=False)
genbank_tax.to_csv("genbank-tax.txt", sep="\t", index=False)
mendota_tax.to_csv("mendota-tax.txt", sep="\t", index=False)
tang_tax.to_csv("tang-tax.txt", sep="\t", index=False)
trout_tax.to_csv("trout-tax.txt", sep="\t", index=False)
ref_tax.to_csv("ref-tax.txt", sep="\t", index=False)

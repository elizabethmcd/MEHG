#! /home/emcd/anaconda3/bin/python3

import re, glob, os
import pandas as pd

# Create master file with dictionary of dictionaries
all_reports = glob.glob("*_assembly_report.txt")
all_dicts = {}
for report in all_reports:
    test_assembly_report=open(report, "r")
    report_dict = {}
    for line in test_assembly_report:
        if line.startswith("#\n"):
            break
        try: 
            descrip = re.match('#\s(.*):', line).group(1).lstrip(' ').replace(' ', '_')
        except: 
            descrip = None
        item = line.split(":", 1)[1].lstrip(' ').replace(' ', '_').rstrip()
        report_dict[descrip] = item
        all_dicts[os.path.basename(report).split("_assembly")[0]] = report_dict
df = pd.DataFrame.from_dict(all_dicts, orient='index', dtype=None)
df.to_csv("reference-genbank-genomes-metadata.csv")

# Filename with organism names
masterfile = pd.read_csv("reference-genbank-genomes-metadata.csv", header=0, index_col=0)
key = masterfile[['GenBank_assembly_accession', 'Organism_name']]
key.rename(columns={'GenBank_assembly_accession':'genome_name', 'Organism_name':'taxonomy'},inplace=True)
key['taxonomy'].replace(regex=True, inplace=True, to_replace=r'\([^)]*\)',  value=r'')
key['taxonomy'].replace(regex=True, inplace=True, to_replace=r'_', value=r' ')
key.to_csv('reference-genbank-genomes-taxonomy.txt', sep='\t', index=False)



    
    
    
    
    

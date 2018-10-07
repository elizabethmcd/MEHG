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
df.to_csv("Masterfile.csv")

# Filename with strain name for tree making
masterfile = pd.read_csv("Masterfile.csv", header=0, index_col=0)
key = masterfile[['Organism_name', 'RefSeq_assembly_accession']]
key.rename(columns={'Organism_name':'name','RefSeq_assembly_accession':'contigs_db_path'},inplace=True)
key['name'].replace(regex=True, inplace=True, to_replace=r'\([^)]*\)',  value=r'')
key['name'].replace(regex=True, inplace=True, to_replace=r'_', value=r' ')
key['contigs_db_path']=key['contigs_db_path'].astype(str)
key['contigs_db_path']=key['contigs_db_path'].str.replace('.','_')
key['contigs_db_path']=key['contigs_db_path'].astype(str) + '_contigs.db'
# strain_names['Organism_name'].replace(regex=True, inplace=True, to_replace=r'\([^)]*\)',  value=r'')
# strain_names['Organism_name'].replace(regex=True, inplace=True, to_replace=r'_', value=r' ')
# strain_names['Filenames']=strain_names['Filenames'].astype(str)
# strain_names['Filenames']=strain_names['Filenames'].str.replace('.', '_')
# strain_names.to_csv('files-and-strains.csv', index=False)
key.to_csv('meth-genomes.txt', sep='\t', index=False)



    
    
    
    
    

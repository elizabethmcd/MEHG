#! /usr/bin/env python

import pandas as pd 

combined_tax = pd.read_table("all-methylators-taxonomy-original-names.txt", sep="\t", names=[['genome_name', 'taxonomy']])

